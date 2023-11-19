# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)

# Load censor data
censor_data <- fread(file='censor_202106.csv')

# Load weekend warrior data
ww <- fread(file='complete_time_052623.csv')

# Load withdrawals
withdrawals <- fread(file='w7089_20230425.csv')

# Merges
setkey(ww,sample_id); setkey(censor_data,sample_id)
censor_data[ww,':='(activity_group = i.activity_group, activity_group_median = i.activity_group_median,
                     accel_age = i.age_accel, sex = i.sex, race = i.race_category_adjust, tob = i.tob,
                     etoh = i.etoh_grams, tdi = i.tdi, employment_status = i.employment_status, 
                     self_health = i.self_health, diet = i.diet, qual_ea = i.qual_ea, accel_date = i.end_date)]

#### Fix censor data in censor file
# Load center categories
center <- fread(file='center0.csv')
center_lookup <- fread(file='enrollment_correspondences.csv')

# Add center value to dataset
setkey(censor_data,sample_id); setkey(center,sample_id)
censor_data[center,':='(center_code = i.value)]

setkey(censor_data,center_code); setkey(center_lookup,Code)
censor_data[center_lookup,':='(center_location = i.Region)]

# Now correct censor dates based on location
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(center_location=='England',phenotype_censor_date,
                                                         ifelse(center_location=='Scotland',pmin(phenotype_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                                pmin(phenotype_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]

# And set censor date to date of death for those who died
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(!is.na(death_date),pmin(death_date,phenotype_censor_date),phenotype_censor_date),origin='1970-01-01'))]

# Remove missing exposure data
censor_data <- censor_data[!is.na(activity_group)] #502485 - 412912 = 89573

# Remove withdrawals
censor_data <- censor_data[!(sample_id %in% withdrawals$V1)] #89573 - 0 = 89573

# Scope columns for PRS analysis
ww <- censor_data[,c('sample_id','activity_group','activity_group_median','accel_age','accel_date',
                      'sex','race','tob','etoh','tdi','employment_status',
                      'self_health','diet','qual_ea','phenotype_censor_date')]
# Write out 
write.csv(ww,file='cox_data_ww.csv',row.names = F)

### ADD BLANKING PERIOD SENSITIVITY ANALYSIS
# Create blanked date (2 years after accelerometer)
ww[,blanked_date := accel_date + 365.25*2]

# Write out
write.csv(ww,file='cox_data_ww_blank_2y.csv',row.names = F)

