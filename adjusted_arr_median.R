# Depends
library(data.table)
library(plyr)
library(stringr)
library(survival)

# Index of significant diseases
sig_either <- fread(file='sig_either_median.csv')

# Index of file names
list <- paste0('/newer_phecode_tables/',sig_either$disease,'.tab.tsv')

# Load exposure/covariate data
ww <- fread('cox_data_ww.csv')
ww[,activity_group_median := factor(activity_group_median,levels=c('Active - Regular','Active - WW','Inactive'))]
setkey(ww,sample_id)

# Format dates
for (j in (c('accel_date','phenotype_censor_date'))){set(ww,j=j,value=as.Date(ww[[j]],format='%Y-%m-%d'))}

# Init vars
outer <- data.table(); n <- 1

# Looping cox model
for (i in list){
  # Merge
  phecode <- NULL; analysis_set <- NULL
  phecode <- read.table(i,sep='\t',header = TRUE); setDT(phecode)
  setkey(phecode,sample_id)
  analysis_set <- ww[phecode,nomatch=0]
  # Format variables
  analysis_set[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]
  # Create analysis variables
  analysis_set[,time_to_event := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - accel_date)/365.25,as.numeric(phenotype_censor_date - accel_date)/365.25),
                                        as.numeric(censor_date - accel_date)/365.25)]
  # Remove prevalent disease or no follow-up
  disease <- analysis_set$disease[1]
  analysis_set <- analysis_set[!is.na(time_to_event) & time_to_event > 0]
  # Fit cox model
  analysis_set[,activity_group_median := factor(activity_group_median,levels=c("Inactive","Active - Regular","Active - WW"))]
  
  model <- coxph(Surv(time_to_event,has_disease) ~ activity_group_median + accel_age + sex + race +
                   tob + etoh + tdi + employment_status + self_health + diet + qual_ea, data=analysis_set)
  
  weights <- data.frame(activity_group_median = levels(analysis_set$activity_group_median),
                        accel_age = rep(mean(analysis_set$accel_age),3),
                        sex = rep('Female',3),
                        race = rep('white',3),
                        tob = rep('Never',3),
                        etoh = rep(mean(analysis_set$etoh),3),
                        tdi = rep(mean(analysis_set$tdi),3),
                        employment_status = rep('Employed',3),
                        self_health = rep('Good',3),
                        diet = rep('intermediate',3),
                        qual_ea = rep(mean(analysis_set$qual_ea),3))
  
  km <- survfit(model,newdata=weights)
  time_index <- km$time - 5
  end_time <- which(time_index == max(time_index[time_index <= 0]))
  obv_est_inactive <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time,1]))(5)
  obv_est2_active_reg <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time,2]))(5)
  obv_est3_active_ww <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time,3]))(5)
  
  ## Bootstrap to obtain standard errors
  out <- data.table(diff_reg=NULL,diff_ww=NULL)
  for (i in 1:200){
    resample <- analysis_set[sample(1:nrow(analysis_set),size=nrow(analysis_set),replace=TRUE)]
    model <- coxph(Surv(time_to_event,has_disease) ~ activity_group_median + accel_age + sex + race +
                     tob + etoh + tdi + employment_status + self_health + diet + qual_ea, data=resample)
    km <- survfit(model,newdata=weights)
    time_index <- km$time - 5
    end_time <- which(time_index == max(time_index[time_index <= 0]))
    est_inactive <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time,1]))(5)
    est2_active_reg <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time,2]))(5)
    est3_active_ww <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time,3]))(5)
    out_iter <- data.table(diff_reg = est2_active_reg - est_inactive,
                           diff_ww = est3_active_ww - est_inactive)
    out <- rbind(out,out_iter)
  }
  obv_diff_reg <- obv_est_inactive - obv_est2_active_reg
  obv_diff_ww <- obv_est_inactive - obv_est3_active_ww
  se_diff_reg <- sd(out$diff_reg)
  se_diff_ww <- sd(out$diff_ww)
  result <- data.table(disease,obv_diff_reg,obv_diff_ww,se_diff_reg,se_diff_ww)
  outer <- rbind(outer,result)
  if (n %% 5 == 0){print(paste0("Just finished model ",n," out of ",length(list),"!"))}
  n <- n+1
}

# Save out
write.csv(outer,file='arr_median.csv',row.names=F)
