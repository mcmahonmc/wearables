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
