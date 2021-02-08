# Find watch off times from heart rate data and interpolate missing steps counts
# Question: How does the way we handle missing actigraphy data affect rest-activity rhythm measures?

## TO DO:
# 1. decide how to interpolate steps data best (mean, or other function)
# 2. write a loop to perform this for all subjects 

library(lubridate)
library(ggplot2)
library(zoo)

# list all heart rate data files and select the first as an example
hr_files <- list.files("~/Box/CogNeuroLab/Wearables/data/fitbit/", pattern = "heartrate_1min", full.names = TRUE)
hr <- read.csv(hr_files[1])

# find fitbit data start and end times and create a sequence in 1 minute intervals
df <- c()
df$Time <- seq(mdy_hms(head(hr[,1], 1)), mdy_hms(tail(hr[,1], 1)), by = "1 min")
hr$Time <- mdy_hms(hr$Time)

# merge sequence with actual fitbit data to identify nans in hr data and associated timestamps
df2 <- merge(df, hr, by = "Time", all = T)
watch_off <- df2$Time[is.na(df2$Value)]

# calculate total time watch was off in minutes
total_time_watch_off <- sum(is.na(df2$Value))

# load in subject's activity data file
act_files <- list.files("~/Box/CogNeuroLab/Wearables/data/fitbit/", pattern = "minuteStepsNarrow", full.names = TRUE)
act_fit <- read.csv(act_files[1])
head(act_fit)
act_fit$ActivityMinute <- mdy_hms(act_fit$ActivityMinute)

# merge heart rate data and steps data
df3 <- merge(df2, act_fit, by.x = "Time", by.y = "ActivityMinute")
head(df3)
sum(is.na(df3$Value))

# step counts at times when watch was off should all be 0
df3$Steps[is.na(df3$Value)]

# set NA values for steps where HR is NA
df3$Steps <- ifelse(is.na(df3$Value), NA, df3$Steps)

# now we want to interpolate step counts at minutes when watch was off
plot_interpolation <- function(df, starttime, endtime, method = "linear", f = NA, maxgap = "none"){
  
  df3$`Steps Interpolated` <- na.approx(df$Steps,  maxgap = maxgap, method = method, f = f)
  
  df3 %>%
    filter(Time > starttime) %>%
    filter(Time < endtime) %>%
    pivot_longer(cols = c(Steps, `Steps Interpolated`), names_to = "Key") %>%
    ggplot() + 
    geom_line(aes(x = Time, y = value, color = Key), size = 2) + 
    facet_wrap(. ~ Key) + theme_classic() + theme(legend.position = "none") +
    ylab("Steps")
}

# let's try out different interpolation methods
df = df3
starttime = ymd_hms("2019-10-25 18:00:00")
endtime = ymd_hms("2019-10-25 19:00:00")

# constant, replacing NAs with median value
method = "constant"
f = 0.5
plot_interpolation(df, starttime, endtime, method, f)

# constant, replacing NAs with upper value
method = "constant"
f = 1
plot_interpolation(df, starttime, endtime, method, f)

# linear interpolation
method = "linear"
f = NA
plot_interpolation(df, starttime, endtime, method, f)

# linear interpolation with a maximum number of interpolated values
# If exceeds maximum gap, data left unchanged.
method = "linear"
f = NA
maxgap = 30 #minutes
plot_interpolation(df, starttime, endtime, method, f, maxgap)

maxgap = 60 #minutes
plot_interpolation(df, starttime, endtime, method, f, maxgap)

# interpolate method #2: find average steps at time when watch was off from other days
library(imputeTS)
plot(df$Steps, type = "l", xlab = "Time", ylab = "Steps", main = "Raw")

# summary of missing periods for which the watch was off
statsNA(df$Steps)

# mean interpolation
plot(na_mean(df$Steps, option = "mean"), 
     type = "l", xlab = "Time", ylab = "Steps", xlim = c(1500,1700), 
     main = "Mean")

# last observation carried forward
plot(na_locf(df$Steps, option = "locf"), 
     type = "l", xlab = "Time", ylab = "Steps", xlim = c(1500,1700), 
     main = "LOCF")

# next observation carried backward
plot(na_locf(df$Steps, option = "nocb"), 
     type = "l", xlab = "Time", ylab = "Steps", xlim = c(1500,1700), 
     main = "NOCB")

# linear interpolation
plot(na.interpolation(df$Steps, option = "linear"), 
     type = "l", xlab = "Time", ylab = "Steps", xlim = c(1500,1700), 
     main = "Linear")

# spline interpolation
plot(na.interpolation(df$Steps, option = "spline"), 
     type = "l", xlab = "Time", ylab = "Steps", xlim = c(1500,1700), 
     main = "Spline")


# method #2

# create new clocktime variable
df$clocktime=lubridate::hour(df$Time) + lubridate::minute(df$Time)/60
df$StepsInt <- df$Steps

# find missing data chunks
find_missing <- function(stepsdata){
  x <- df %>%
    mutate(missing = ifelse(is.na(stepsdata), missing, 0)) %>%
    group_by(group = cumsum(c(0, diff(missing) != 0))) %>%
    filter(missing == 1 & n() > 1) %>%
    summarize("start_missing"=min(as.character(Time)),
              "end_missing"=max(as.character(Time)),
              "length_missing"=n()) %>%
    ungroup() %>%
    select(-matches("group"))
  
  # create new cloktime variable
  x$startclock=lubridate::hour(x$start_missing) + lubridate::minute(x$start_missing)/60
  x$endclock=lubridate::hour(x$end_missing) + lubridate::minute(x$end_missing)/60
  x <- x[order(x$length_missing),] 
  return(x)
}

# this shows us the data during the missing period of interest
i = 1 # looking at first missing period

# view missing data period
df %>%
  filter(Time > ymd_hms(x$start_missing[i])) %>%
  filter(Time < ymd_hms(x$end_missing[i]))

# now we want to loop through each missing period, interpolate, update what periods are missing, 
# and interpolate some more, using the average value from the same time period on other days 
# of recording

# create new dataframe
df$StepsInt <- df$Steps
x <- find_missing(df$StepsInt)

while (dim(x)[1] > 0) {
  # convert start and end times to clocktimes
  startclock=lubridate::hour(x$start_missing[1]) + lubridate::minute(x$start_missing[1])/60
  endclock=lubridate::hour(x$end_missing[1]) + lubridate::minute(x$end_missing[1])/60
  
  # get index values for missing period
  replaceindex <- which((df$Time >= ymd_hms(x$start_missing[1])) & (df$Time <= ymd_hms(x$end_missing[1])))
  
  # now we want to see the steps data during the time period of interest on all other days
  mean_steps <- df %>%
    filter(clocktime >= startclock) %>%
    filter(clocktime <= endclock) %>%
    summarise(mean_steps = mean(Steps, na.rm = T)) %>%
    unlist()
  
  # replace values with the mean steps value from all other time periods
  df$StepsInt[replaceindex] <- mean_steps
  
  # checking that we replaced them all
  missing <- df %>%
    filter(Time > ymd_hms(x$start_missing[1])) %>%
    filter(Time < ymd_hms(x$end_missing[1])) %>%
    summarise(missing = sum(is.na(StepsInt))) %>%
    unlist()
  
  if (missing > 0){
    print("error - detected missing values!")
  } else {
    print(paste0("all missing values replaced, ", sum(is.na(df$StepsInt)), " remaining"))
  }
  
  # check for missing periods and update on each round of interpolation
  x <- find_missing(df$StepsInt)
  print(dim(x)[1])
}

df %>%
  mutate(Date = lubridate::date(Time)) %>%
  mutate(Hour = lubridate::hour(Time)) %>%
  select(Hour, Date, Steps, StepsInt) %>%
  melt(id.vars = c("Hour", "Date")) %>%
  ggplot() + 
  geom_point(aes(x = Hour, y = value, color = variable, alpha = 0.85)) + 
  facet_wrap(Date ~ variable) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() + 
  xlab("Time") + ylab("Steps") +
  ggsave("~/Box/CogNeuroLab/Wearables/results/figures/interpolated_ts_fitbit.png", dpi = 300)

# to do: highlight all missing periods originally found in find_missing function
x0 <- find_missing(df$Steps)

