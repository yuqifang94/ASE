source('mainFunctions_sub.R')
#Do it per sample
merged_dt=readRDS("../downstream/output/human_analysis/CPEL_outputs/NME_MML_merged_dt.rds")
merged_dt_sc=merged_dt[!is.na(hyper_var_fn)]
out=data.table()
for(sp in  unique(merged_dt_sc$Sample)){
    fn_all=unique(merged_dt_sc[Sample==sp]$hyper_var_fn)
    for(fn_in in unlist(strsplit(fn_all,';'))){
        scRNA_in=readRDS(fn_in)
        out=rbind(out,data.table(fn=fn_in,Sample=sp,
                                    R2_scRNA=summary(lm(log(scRNA_in$var)~log(scRNA_in$mean)))$r.squared,
                                    R2_methylation=summary(lm(merged_dt_sc[Sample==sp]$NME~merged_dt_sc[Sample==sp]$MML))$r.squared
                                    ))
    }
    
}
write.csv(out,'../downstream/output/human_analysis/QC/out_R2_scRNA_methylation.csv')