# Script to create incidence waterfalls using pre-processed data

# Depends
library(data.table)
library(viridis)
library(stringr)
library(fdrtool)
library(EValue)

# Load processed results
reg <- fread('regular_150_ref_inactive.csv')
ww <- fread('ww_150_ref_inactive.csv')

# Sigs
reg <- reg[sig==1]
ww <- ww[sig==1]

# Find unique sigs
reg_but_not_ww <- reg[!(phecode %in% ww$phecode)]
ww_but_not_reg <- ww[!(phecode %in% reg$phecode)]

# Write out
write.csv(reg_but_not_ww,file='reg_but_not_ww.csv',row.names=F)
write.csv(ww_but_not_reg,file='ww_but_not_reg.csv',row.names=F)