# Script to process phecode tables

# Depends
library(data.table)
library(stringr)
library(survival)

# Index of file names
list <- paste0('/newer_phecode_tables/',list.files('/newer_phecode_tables'))

# Load exposure/covariate data
ww <- fread('cox_data_ww.csv')
setkey(ww,sample_id)

# Format dates
for (j in (c('accel_date','phenotype_censor_date'))){set(ww,j=j,value=as.Date(ww[[j]],format='%Y-%m-%d'))}

# Init vars
out <- data.table(); n <- 1

# Set levels
ww[,activity_group_median := factor(activity_group_median,levels=c('Active - Regular','Active - WW','Inactive'))]

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
  analysis_set <- analysis_set[!is.na(time_to_event) & time_to_event > 0]
  # Define events and follow-up
  disease <- analysis_set$disease[1]
  n_events <- nrow(analysis_set[has_disease==1])
  fu_median <- quantile(analysis_set$time_to_event,0.50); fu_q1 <- quantile(analysis_set$time_to_event,0.25); fu_q3 <- quantile(analysis_set$time_to_event,0.75)
  # If less than 10 cases, abort
  if (n_events < 10){
    hr_ww <- NA; lower_ww <- NA; upper_ww <- NA; z_ww <- NA; p_ww <- NA
    result <- data.table(disease,n_events,fu_median,fu_q1,fu_q3,hr_ww,lower_ww,upper_ww,z_ww,p_ww)
    out <- rbind(out,result)
    print(paste0("Skipping phenotype ",analysis_set$disease[1]," since < 10 cases"))
    if (n %% 50 == 0){print(paste0("Just finished model ",n," out of ",length(list),"!"))}
    n <- n+1; next}
  # Fit cox model
  model <- coxph(Surv(time_to_event,has_disease) ~ activity_group_median + accel_age + sex + race +
                   tob + etoh + tdi + employment_status + self_health + diet + qual_ea, data=analysis_set)
  
  hr_ww <- summary(model)$coefficients[1,2]; lower_ww <- summary(model)$conf.int[1,3]
  upper_ww <- summary(model)$conf.int[1,4]; z_ww <- summary(model)$coefficients[1,4]
  p_ww <- 2*pnorm(abs(z_ww),lower.tail=FALSE)
  
  result <- data.table(disease,n_events,fu_median,fu_q1,fu_q3,hr_ww,lower_ww,upper_ww,z_ww,p_ww)
  out <- rbind(out,result)
  if (n %% 50 == 0){print(paste0("Just finished model ",n," out of ",length(list),"!"))}
  n <- n+1
}

# Save out
write.csv(out,file='cox_ww_ref_regular.csv',row.names=F)