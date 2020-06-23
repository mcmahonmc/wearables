#fitbit study
act_raw <- "~/Box/CogNeuroLab/Wearables/data/circadian_measures/raw/actiwatch/"
fit_raw <- "~/Box/CogNeuroLab/Wearables/data/circadian_measures/raw/fitbit/"
out <- "~/Box/CogNeuroLab/Wearables/data/circadian_measures/7_days/"


truncate_all <- function(in_dir, out_dir, ndays){
  library(lubridate)
  # truncate to number of days desired for recording period
  for (f in list.files(in_dir)){
    d <- read_delim(paste0(in_dir, f), delim = " ", col_names = F)
    subject <- substr(f, 1, 5)
    start = ymd_hms(tail(d$X1, 1), tz="UTC") - days(ndays)
    d_truncated <- d[d$X1 >= start,]
    write.table(d_truncated, paste0(out_dir, subject, "_", ndays, "days.txt"), sep = " ", row.names = F, col.names = F)
  }
}

calc_nparact <- function(in_dir, out_dir, SR, device){
  #device = "act" or "fit"
  #SR = 1/60 for validation study if actiwatch has already been resampled from python script
  library(nparACT)
  actall <- nparACT_base_loop(in_dir, SR = SR, fulldays = F)
  actall$record_id <- substr(list.files(in_dir), 1, 5)
  write.csv(actall, paste0(out_dir, "nparact_7days_", device, ".csv"), row.names = F)
}

run_cosinor <- function(d, filename, proc_dir = out_dir, print = FALSE) {
  # modified from Stephanie Sherman
  
  results <- data.frame(stringsAsFactors = FALSE)
  d$record_id <- substr(filename, 1, 5)
  d$cloktime=lubridate::hour(d$time) + lubridate::minute(d$time)/60
  
  if (sum(d$cloktime) != 0) {
    d$twopio24 = (2*3.14159)/24 
    d$xcos = cos(d$twopio24*d$cloktime) 
    d$xsin = sin(d$twopio24*d$cloktime)
    
    #d$activity=as.character(d$ZCM)
    #d$activity=as.numeric(d$PIM, 'NA')
    d$lactivity = log((d$activity +1),10)
    
    allwatch=d[,c('record_id','cloktime','lactivity','xcos','xsin','twopio24')]
    allwatch=na.omit(allwatch)
    
    model=lm(allwatch$lactivity ~ allwatch$xcos + allwatch$xsin)
    allwatch$linactxb=coef(model)['(Intercept)']
    allwatch$linactcos=coef(model)['allwatch$xcos']
    allwatch$linactsin=coef(model)['allwatch$xsin']
    #need column for residuals called linract
    allwatch$linract=model$residuals
    
    # filename = paste0(work_dir, '/residuals/', subject, '_residuals.csv')
    # write.csv(allwatch, file = filename, row.names = FALSE)
    
    actres1 <- allwatch
    
    actres1$linactamp = sqrt(actres1$linactcos^2 + actres1$linactsin^2)
    actres1$linactmin = actres1$linactxb-actres1$linactamp 
    
    for (p in 1:length(actres1$lactivity[1])){
      if (actres1$linactsin[1] > 0 & actres1$linactcos[1] > 0) {
        actres1$phase = atan(actres1$linactsin/actres1$linactcos)}
      else if (actres1$linactsin[1] > 0 & actres1$linactcos[1] < 0) {
        actres1$phase = 3.14159 - atan(actres1$linactsin/abs(actres1$linactcos))}
      else if (actres1$linactsin[1] < 0 & actres1$linactcos[1] < 0) {
        actres1$phase = 3.14159 + atan(abs(actres1$linactsin)/abs(actres1$linactcos))}
      else {(actres1$linactsin[1] < 0 & actres1$linactcos[1] > 0)
        actres1$phase = 2*3.14159 - atan(abs(actres1$linactsin)/(actres1$linactcos))} 
    }
    
    actres1$linactacro = actres1$phase*24/(2*3.14159) 
    
    #get sum of squares (uss variable)
    linractuss=(sum((actres1$linract)^2))-((sum(actres1$linract))^2/(length(actres1$linract))) 
    
    #num_nonmissingvalues
    nlinract=dim(actres1)[1]
    
    #nonlinear regression
    carhythm = function(actphi,actbeta,actalph,actmin,actamp,cloktime) {
      twopio24 = (2*3.14159)/24 
      rhythm = cos(twopio24*(cloktime - actphi ))
      lexpt=actbeta*(rhythm - actalph)
      expt = exp(lexpt)
      er = expt/(1 + expt)
      actmin + actamp*er
      
    }
    
    #if want it to print out iterations change trace=TRUE
    error = try(b <- nls(actres1$lactivity ~carhythm(actphi,actbeta,actalph,actmin,actamp,cloktime),
                         data=actres1, algorithm='port',
                         start=list(actphi = 12,actbeta = 2.00,actalph = 0.0,actmin =0,actamp=1),
                         lower=list(actphi = -3,actbeta = 0,actalph = -1,actmin =0,actamp=1),
                         upper=list(actphi = 27,actbeta = Inf,actalph = 1,actmin =Inf,actamp=5),
                         control=list(maxiter=200), #warnOnly=TRUE
                         trace=FALSE))
    print(error)
    
    if(class(error)!="try-error"){
      actres1$rnlact=resid(b)
      actres1$pnlact=fitted(b)	
      
      
      # take estimates from model and add to actres (in SAS all5) changes parameter names
      ## x beginning variables are the same as the e beginning variables
      actres1$xactphi=coef(b)['actphi']
      actres1$xactbeta=coef(b)['actbeta']
      actres1$xactalph=coef(b)['actalph']
      actres1$xactmin=coef(b)['actmin']
      actres1$xactamp=coef(b)['actamp']
      
      actres1$coact = actres1$linactxb + actres1$linactcos*actres1$xcos + actres1$linactsin*actres1$xsin
      
      ncssrnlact=(sum((actres1$rnlact)^2))-((sum(actres1$rnlact))^2/(length(actres1$rnlact)))
      cssact=(sum((actres1$lactivity)^2))-((sum(actres1$lactivity))^2/(length(actres1$lactivity)))
      nact=length(actres1$lactivity)
      nlinract=length(actres1$lactivity) 
      
      
      actacos=acos(actres1$xactalph[1])/actres1$twopio24[1]
      acthalftimel=-actacos + actres1$xactphi[1]
      acthalftimer=actacos + actres1$xactphi[1]
      actwidthratio = 2*actacos/24
      
      
      if(actres1$xactalph[1] < -0.99 |actres1$xactalph[1] > 0.99){
        actwidthratio = 0.5
        acthalftimel = (actres1$xactphi[1] - 6)
        acthalftimer = actres1$xactphi[1] + 6
      }
      
      actdervl = -sin((acthalftimel - actres1$xactphi[1])*actres1$twopio24[1])
      actdervr = -sin((acthalftimer - actres1$xactphi[1])*actres1$twopio24[1])	
      
      #sd is standard error I can get that from nls output 
      sdactphi=summary(b)$coefficients['actphi',2]
      sdactbeta=summary(b)$coefficients['actbeta',2]
      sdactalph=summary(b)$coefficients['actalph',2]
      sdactmin=summary(b)$coefficients['actmin',2]
      sdactamp=summary(b)$coefficients['actamp',2]
      
      #t is t value from model
      tactphi=summary(b)$coefficients['actphi',3]
      tactbeta=summary(b)$coefficients['actbeta',3]
      tactalph=summary(b)$coefficients['actalph',3]
      tactmin=summary(b)$coefficients['actmin',3]
      tactamp=summary(b)$coefficients['actamp',3]
      
      rsqact = (cssact - ncssrnlact)/cssact  
      fact = ((cssact - ncssrnlact)/4)/(ncssrnlact/(nlinract - 5))
      ndf = 4
      ddfact = nlinract - 5
      efact = ddfact/(ddfact - 2)
      varfact = ( 2/ndf )*( efact**2 )*( (ndf + ddfact -2)/(ddfact - 4) )  #wilks p. 187 */;
      tfact = (fact - efact)/sqrt(varfact)
      varact = cssact/(nlinract - 1)
      mselinact = linractuss/(nlinract - 3)
      msenlinact = (ncssrnlact/(nlinract - 5))
      fnlrgact = ((linractuss - ncssrnlact)/2)/(ncssrnlact/(nlinract - 5)) 
      flinact = ((cssact - linractuss)/2)/(linractuss/(nlinract - 3)) 
      
      actmesor = actres1$xactmin[1] + (actres1$xactamp[1]/2) 
      actupmesor = acthalftimel
      actdownmesor = acthalftimer 
      actamp=actres1$xactamp[1]
      actbeta=actres1$xactbeta[1]
      actphi=actres1$xactphi[1]
      actmin=actres1$xactmin[1]
      actalph=actres1$xactalph[1]
      session=actres1$session[1]
      record_id=actres1$record_id[1]
      rhythm=as.character(c(record_id, actamp,actbeta,actphi,actmin,actmesor,actupmesor,actdownmesor,actalph,actwidthratio,rsqact,fact,fnlrgact))
      newline <- data.frame(t(rhythm), stringsAsFactors = FALSE)
      #results <- rbind(results, newline)
      return(newline)
    }else{
      print(paste0("Unable to obtain cosinor model for subject ", d$record_id[1]))
      newline <- c(d$record_id[1], rep(NA, 12))
      return(newline)
    }
    
  }
}

calc_cosinor <- function(in_dir, out_dir, device){
  library(readr)
  library(lubridate)
  results <- c()
  
  for (filename in list.files(in_dir)){
    print(filename)
    d <- read_delim(paste0(in_dir, filename), delim = " ")
    colnames(d) <- c("date", "time", "activity")
    d$time <- with(d, ymd(date) + hms(time))
    rhythm <- run_cosinor(d, filename, proc_dir)
    results <- rbind(results, rhythm)
    
  }
  colnames(results)=c('record_id','actamp','actbeta','actphi','actmin','actmesor','actupmesor','actdownmesor','actalph','actwidthratio','rsqact','fact','fnlrgact')
  write.csv(results, file = paste0(out_dir, "cosinor_7days_", device, ".csv"), row.names = F)
  return(results)
}



truncate_all(in_dir=act_raw, out_dir=paste0(out_dir, "actiwatch/"), 7)
calc_nparact(in_dir=paste0(out, "actiwatch/"), out_dir=out, SR=1/60, device = "act")
calc_cosinor(in_dir=paste0(out, "actiwatch/"), out_dir=out, device = "act")

truncate_all(in_dir=fit_raw, out_dir=paste0(out_dir, "fitbit/"), 7)
calc_nparact(in_dir=paste0(out, "fitbit/"), out_dir=out, SR=1/60, device = "fit")
calc_cosinor(in_dir=paste0(out, "fitbit/"), out_dir=out, device = "fit")
