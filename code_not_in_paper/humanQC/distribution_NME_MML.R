source('mainFunctions_sub.R')

GR_merge=readRDS(GR_merge_file)
CpG_hg19=getCpgSitesH19()#26752702
subsetByOverlaps(CpG_hg19,GR_merge[GR_merge$dMML_pval<=0.1])#25105
subsetByOverlaps(CpG_hg19,NME_all)#21234156
MML_all=readRDS(MML_agnostic_all_file)
#Remove ASM filtered NME region
MML_all_dt=MML_all_dt[!is.na(MML)]
MML_all_dt=convert_GR(MML_all,direction="DT")
MML_all_dt_rmASM=MML_all_dt[paste0(region,Sample)%in%paste0(NME_all_dt$region,NME_all_dt$Sample)]
saveRDS(MML_all_dt_rmASM,"../downstream/output/human_analysis/CPEL_outputs/MML_all_dt_rmASM.rds")
NME_all=readRDS(NME_agnostic_all_file)
NME_all_dt=convert_GR(NME_all,direction="DT")
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.text.x=element_text(size=7),
                                 axis.text.y=element_text(size=7))
# Plotting NME vs MML ------------------------------------
NME_all_dt=NME_all_dt[!is.na(NME)]
MML_all_dt_rmASM=readRDS("../downstream/output/human_analysis/CPEL_outputs/MML_all_dt_rmASM.rds")
merged_dt=data.table(MML=MML_all_dt_rmASM$MML,NME=NME_all_dt$NME,Sample=NME_all_dt$Sample,hyper_var_fn=NME_all_dt$hyper_var_fn)#number of MML: 135135992, number of NME: 135135992
saveRDS(merged_dt,"../downstream/output/human_analysis/CPEL_outputs/NME_MML_merged_dt.rds")
#Aggregate NME, using quantiles, 0.05 and 0.95
merged_dt_agg=merged_dt[, list(NME=round(median(NME),digits=digits_round),
                                   Bottom25=round(quantile(NME,probs=0.25),digits=digits_round),
                                   top25=round(quantile(NME,probs=0.75),digits=digits_round)), 
                            by = list(MML = round(MML,digits=digits_round))]

merged_dt_agg$Bottom25= predict(loess(Bottom25~MML,merged_dt_agg),newdata=merged_dt_agg$MML)
merged_dt_agg$top25= predict(loess(top25~MML,merged_dt_agg),newdata=merged_dt_agg$MML)
pdf('../downstream/output/human_analysis/QC/NME_vs_MML_with_quantile_agnostic.pdf',width=5,height=5)
#Plotting
print(ggplot(merged_dt_agg,aes(x=MML, y=NME))+
        xlim(c(0,1))+ylim(c(0,1))+ggtitle("MML and NME relationship")+geom_smooth(method="loess",se=FALSE)+
        ylab("NME")+geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)+theme_glob)
dev.off()
subsample=3000000
pdf('../downstream/output/human_analysis/QC/NME_vs_MML_bin2d.pdf',width=5,height=5)
#Plotting
set.seed(1)
print(ggplot(merged_dt[sample(1:nrow(merged_dt),subsample)],aes(x=MML, y=NME))+
        xlim(c(0,1))+ylim(c(0,1))+ggtitle("MML and NME relationship")+
        #geom_smooth(method="loess",se=FALSE)+
        ylab("NME")+geom_bin2d(bins=100)+theme_glob+
        scale_fill_gradientn(name = "count", trans = "log10",colours=rainbow(6)))
dev.off()


#Fstat:
sqare_m_NME=lm(merged_dt$NME~ poly(merged_dt$MML, 2))
summary(sqare_m_NME)
NME_all_dt$NME_predict=predict(sqare_m_NME, NME_all_dt)
NME_all_dt$NME_corrected=NME_all_dt$NME-NME_all_dt$NME_predict
saveRDS(NME_all_dt,"../downstream/output/human_analysis/CPEL_outputs/NME_all_dt_corrected.rds")
# Call:
# lm(formula = merged_dt$NME ~ poly(merged_dt$MML, 2))

# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.87087 -0.02867  0.02082  0.06093  0.14306

# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)
# (Intercept)              4.553e-01  8.005e-06   56877   <2e-16 ***
# poly(merged_dt$MML, 2)1 -6.424e+02  9.305e-02   -6904   <2e-16 ***
# poly(merged_dt$MML, 2)2 -2.744e+03  9.305e-02  -29491   <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.09305 on 135135989 degrees of freedom
# Multiple R-squared:  0.8716,    Adjusted R-squared:  0.8716
# F-statistic: 4.587e+08 on 2 and 135135989 DF,  p-value: < 2.2e-16
linear_m_NME=lm(merged_dt$NME~poly(merged_dt$MML, 1))
summary(linear_m_NME)
# Call:
# lm(formula = merged_dt$NME ~ poly(merged_dt$MML, 1))

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.61068 -0.17772 -0.00985  0.20002  0.49693 

# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             4.553e-01  2.183e-05   20858   <2e-16 ***
# poly(merged_dt$MML, 1) -6.424e+02  2.537e-01   -2532   <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.2537 on 135135990 degrees of freedom
# Multiple R-squared:  0.04528,   Adjusted R-squared:  0.04528 
# F-statistic: 6.41e+06 on 1 and 135135990 DF,  p-value: < 2.2e-16

sqare_m=lm(merged_dt$NME~ poly(merged_dt$MML, 2))
summary(sqare_m)
# Call:
# lm(formula = GR_merge$dNME ~ poly(GR_merge$dMML, 2))

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.34018 -0.05208 -0.00898  0.05313  0.80038 

# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              1.530e-01  5.435e-05  2815.5   <2e-16 ***
# poly(GR_merge$dMML, 2)1  1.231e+02  9.922e-02  1241.1   <2e-16 ***
# poly(GR_merge$dMML, 2)2 -6.473e+01  9.922e-02  -652.4   <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.09922 on 3332741 degrees of freedom
# Multiple R-squared:  0.371,     Adjusted R-squared:  0.371 
# F-statistic: 9.83e+05 on 2 and 3332741 DF,  p-value: < 2.2e-16

linear_m=lm(GR_merge$dNME~poly(GR_merge$dMML, 1))
summary(linear_m)

# Call:
# lm(formula = GR_merge$dNME ~ poly(GR_merge$dMML, 1))

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.87531 -0.06744 -0.01448  0.05932  0.69739 

# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            1.530e-01  5.771e-05    2651   <2e-16 ***
# poly(GR_merge$dMML, 1) 1.231e+02  1.054e-01    1169   <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 0.1054 on 3332742 degrees of freedom
# Multiple R-squared:  0.2907,    Adjusted R-squared:  0.2907 
# F-statistic: 1.366e+06 on 1 and 3332742 DF,  p-value: < 2.2e-16
cor.test(GR_merge$dNME,GR_merge$dMML)

#         Pearson's product-moment correlation

# data:  GR_merge$dNME and GR_merge$dMML
# t = 1168.7, df = 3332742, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5383987 0.5399217
# sample estimates:
#       cor 
# 0.5391607 

cor.test(merged_dt$NME,merged_dt$MML)

#         Pearson's product-moment correlation

# data:  merged_dt$NME and merged_dt$MML
# t = -2531.7, df = 135135990, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2129591 -0.2126371
# sample estimates:
#        cor 
# -0.2127981 

#Correlation with CpG density

#correlation between NME, MML and dNME, dMML separating by sample
merged_dt_r2_linear=merged_dt[,list(R2=summary(lm(NME~poly(MML,1)))$r.squared),by=list(Sample)]
merged_dt_r2_poly=merged_dt[,list(R2=summary(lm(NME~poly(MML,2)))$r.squared),by=list(Sample)]
GR_merge_dt=convert_GR(GR_merge,direction="DT")
GR_merge_r2_linear=GR_merge_dt[dMML_pval<=pval_cutoff&dNME_pval<=pval_cutoff,
                                list(R2=ifelse(length(dMML)>10&length(dNME)>10,summary(lm(dNME~poly(dMML,1)))$r.squared,as.numeric(NA)),
                                corP=ifelse(length(dMML)>10&length(dNME)>10,cor(dNME,dMML),as.numeric(NA))
                                #corP=cor(dNME,dMML)                                
                                ),by=list(Sample)]
