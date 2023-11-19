# # Script to create incidence waterfalls using pre-processed data

# Depends
library(data.table)
library(viridis)
library(stringr)
library(fdrtool)
library(EValue)

# Loads
dict <- fread(file='phecode_definitions1.2.csv')
results <- fread(file='v1_cox_ww_weekend_ref_inactive_150.csv')

# Load color correspondences from Fig 2
col_corr <- fread(file='col_corr.csv')

# Convert output format to phecode
results[,phecode := as.numeric(str_remove(str_replace_all(disease,'\\_','\\.'),'phecode\\.'))]

# Join on phecode to get names
setkey(dict,phecode); setkey(results,phecode)
results[dict,':='(phenotype = i.phenotype,
                  category = i.category)]

# Apply post-hoc removal of congenital anomalies and pregnancy related conditions
results <- results[!(category %in% c('congenital anomalies','pregnancy complications'))]

# Apply post-hoc event filter (N=1754 conditions prior to exclusion, N=697 after exclusion)
results <- results[n_events>=120]

# Remove NAs
results <- results[!is.na(hr_ww)]

# Fix p-values
results[,p_ww := 2*pnorm(abs(z_ww),lower.tail = FALSE)]

# FDR
fdr_out <- fdrtool(x=results$p_ww,statistic='pvalue')

results[which(fdr_out$qval<0.01),sig := 1]

# E-values
results[,rare_dz := ifelse(n_events/89573 <= 0.15,1,0)]
evals <- data.table(e_point=NULL, e_null=NULL)
for (i in 1:nrow(results)){
  ev <- evalues.HR(est=results$hr_ww[i],lo=results$lower_ww[i],hi=results$upper_ww[i],rare=results$rare_dz[i])
  result <- data.table(e_point=ev[2,1],e_null=ev[2,!is.na(ev[2,])][2])
  evals <- rbind(evals,result)
  i <- i+1
}

results[,':='(e_point = evals$e_point,e_null = evals$e_null)]

# Set up plotting variables
setkey(results,p_ww)
correction <- nrow(results)

# Unite miscellaneous category
results[is.na(category) | category=='NULL']$category <- 'other'

# Create dummy variable to control category order
results[,category_order := ifelse(category=='other','a_other',category)]

# color
setkey(results,category); setkey(col_corr,all_categories)
results[col_corr,col := i.cat_col]

# shape
results[,shape := ifelse(hr_ww < 1,6,2)]

# Sort by category then randomly within category for better spread
setkeyv(results,c('category_order'))
results <- results[,.SD[sample(.N)],by='category']

# Create vector of midpoints per category for x-axis labels
x_counts <- results[,.N,by='category_order']

midpoint <- function(x){
  if(x %% 2==0){result <- x %/% 2
  } else {result <- (x %/% 2) + 1}
}

x_locs <- data.table(category=x_counts$category_order,x_count=sapply(x_counts$N,midpoint))

out <- rep(NA,length(x_locs$x_count))
for (i in 1:length(x_locs$x_count))
  if (i == 1){out[i] <- x_locs$x_count[i]
  } else {
    out[i] <- x_locs$x_count[i] + cumsum(x_counts$N)[i-1]
  }

x_locs$x_coord <- out

setkey(x_locs,category)

# Graphical y
results[,p_graphical := ifelse(p_ww < 1*10^-20,1*10^-20,p_ww)]

# Plot p-values
pdf(file='ww_weekend_150_ref_inactive.pdf',height=5,width=12,pointsize=5)
par(mar=c(15,3,1,6),oma=c(1,1,1,1))

plot(x=1:nrow(results),y=-log10(results$p_graphical),col=ifelse(!is.na(results$sig),results$col,paste0(results$col,'4D')),
     bty='n',xaxt='n',yaxt='n',xlim=c(0,nrow(results)),
     ylim=c(0,20),xlab='',ylab='',cex=2,pch=results$shape)

axis(1,at=x_locs$x_coord,cex.axis=2,labels=rep('',length(x_locs$x_coord)))
axis(2,cex.axis=2,at=seq(0,20,2),las=2,pos=-12)

mtext("-log(p)",2,line=1.5,cex=2)

segments(0,-log10(max(results[sig==1]$p_ww)),nrow(results),-log10(max(results[sig==1]$p_ww)),col='#bd0026',lty=5)

text(x = x_locs$x_coord,
     y = par("usr")[3] - 0.8,
     labels = str_remove_all(x_locs$category,'[A-z]\\_'),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 2)

dev.off()

# Plot effect sizes
setkey(results,hr_ww)

# color
sig_results <- results[sig==1]
sig_results[is.na(phenotype)]$phenotype <- 'Eating disorders'
sig_results$color <- rep('darkgray',nrow(sig_results))
color <- viridis(nrow(sig_results[sig==1]))
for (i in 1:length(color)){
  sig_results$color[i] <- color[i]
}

pdf(file='ww_weekend_150_ref_inactive_effect.pdf',height=7,width=10,pointsize=5)
par(mar=c(32,3,1,8),oma=c(1,1,1,1))

plot(x=1:nrow(sig_results),y=sig_results$hr_ww,col=sig_results$col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,nrow(sig_results)),
     ylim=c(0.9,1.1),xlab='',ylab='',cex=2.5)

axis(1,at=1:nrow(sig_results),cex.axis=2,labels=rep('',nrow(sig_results)))
axis(2,cex.axis=2.2,at=seq(0.9,1.1,0.05),las=2,pos=0.5)

mtext("Hazard ratio",2,line=0.5,cex=2.8)

segments(0.6,1,nrow(sig_results),1,col='#bd0026',lty=2)

segments(1:nrow(sig_results),sig_results$lower_ww,1:nrow(sig_results),sig_results$upper_ww,col=sig_results$col,
         lwd=2)

text(x = 1:nrow(sig_results),
     y = par("usr")[3] - 0.008,
     labels = sig_results$phenotype,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.8)

dev.off()

write.csv(results,'ww_weekend_150_ref_inactive.csv',row.names=F)

# Noodles
sig_good <- sig_results[hr_ww < 1]
sig_bad <- sig_results[hr_ww > 1]
write.csv(sig_bad,file='ww_weekend_150_ref_inactive_bad.csv',row.names=F)

