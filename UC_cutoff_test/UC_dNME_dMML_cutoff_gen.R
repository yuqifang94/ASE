# nolint start
#This file is to generating the overlap between dNME-selected region (or dMML-selected region) and UC-selected region under different UC, dNME and dMML cutoff
source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
dNME_in=mclapply(UC_in,function(x) x[,grepl('dNME',colnames(x))],mc.cores=7)
dMML_in=mclapply(UC_in,function(x) x[,grepl('dMML',colnames(x))],mc.cores=7)
UC_in_only=mclapply(UC_in,function(x) x[,grepl('UC',colnames(x))],mc.cores=7)
cutoff_selection_only<-function(dat_in,ts,ts_only=T,step_in=0.0001#,
    ){
   
    if(ts_only){return(mclapply(seq(0,1,step_in),function(cutoff,ts) {
        #cutoff=cutoff+step*diff
        cat("Processing:",cutoff,'\n')
        aid <- sapply(names(dat_in),function(i) {
                    names(which(rowSums(dat_in[[i]] > cutoff) > 0))
                    })          
        out=list()
        return(setdiff(aid[[ts]],unlist(aid[names(aid)!=ts])))

    },ts=ts,mc.cores=20
    ))
    }else{return(
        mclapply(seq(0,1,step_in),function(cutoff,ts) {
        #cutoff=cutoff+step*diff
        cat("Processing:",cutoff,'\n')

        return(names(which(rowSums(dat_in[[ts]] > cutoff) > 0)))
        },ts=ts,mc.cores=20
        ))
    }

}
assign_name<-function(x,step_in=0.0001) {names(x)=as.character(seq(0,1,step_in)); return(x)}


#generate for one and only one 
dNME_cutoff=list()
dMML_cutoff=list()
UC_cutoff=list()
for(ts in names(dNME_in)){
    cat("Processing:",ts,"\n")
    dNME_cutoff[[ts]]=cutoff_selection_only(dNME_in,ts)
    
    dMML_cutoff[[ts]]=cutoff_selection_only(dMML_in,ts)
    UC_cutoff[[ts]]=cutoff_selection_only(UC_in_only,ts)

}

saveRDS(list(UC_cutoff=lapply(UC_cutoff,assign_name),
            dMML_cutoff=lapply(dMML_cutoff,assign_name),
            dNME_cutoff=lapply(dNME_cutoff,assign_name)
            ),'../downstream/output/mouse_analysis/UC_dNME_olap/regions_ts_only.rds')
UC_in=readRDS(UC_merge_file)
#For each possible cutoff, find the overlap
region_rand_pool=lapply(UC_in,rownames)
rm(UC_in)
cutoff_ts_only_out=readRDS('../downstream/output/mouse_analysis/UC_dNME_olap/regions_ts_only.rds')
cutoff_ts_num_output_all=list()
set.seed(12345)
 pdf(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/UC_dNME_dMML_cutoff_num.pdf'))
for(ts in names(region_rand_pool)){
    
    tt1=proc.time()[[3]]
    
    cutoff_ts_only_out_ts=lapply(cutoff_ts_only_out,function(x) x[[ts]])
    ts_num=lapply(cutoff_ts_only_out_ts,function(x) unlist(lapply(x,length)))
    ts_num_dt=data.table(cutoff=names(ts_num$UC))
    ts_num_dt$dNME=ts_num$dNME[ts_num_dt$cutoff]
    ts_num_dt$dMML=ts_num$dMML[ts_num_dt$cutoff]
    ts_num_dt$UC=ts_num$UC[ts_num_dt$cutoff]
    ts_num_dt_mt=melt.data.table(ts_num_dt,id.vars='cutoff',variable.name='stat_type',value.name='region_count')
   ts_num_dt_mt$cutoff=as.numeric(ts_num_dt_mt$cutoff)
    print(ggplot(ts_num_dt_mt,aes(x=cutoff,y=region_count,color=stat_type))+geom_line()+
            ylab("Number of regions")+ggtitle(ts))
    cat("Processing:",ts,'\n')
    cutoff_ts_num_output=as.data.table(expand.grid(seq(0,max(ts_num_dt[UC>=round(mean(ts_num_dt$UC)/100,digits=2)]$cutoff),0.005),
                                        seq(0,max(ts_num_dt[dNME>=round(mean(ts_num_dt$dNME)/100,digits=2)]$cutoff),0.01),
                                        seq(0,max(ts_num_dt[dMML>=round(mean(ts_num_dt$dMML)/100,digits=2)]$cutoff),0.01)))
    colnames(cutoff_ts_num_output)=c("UC","dNME","dMML")

    cutoff_ts_num_output_all[[ts]]=
    cbind(cutoff_ts_num_output,
        fastDoCall("rbind",
                    mclapply(1:nrow(cutoff_ts_num_output),
                                function(x,cutoff_ts_only_out_ts,cutoff_ts_num_output) {
                                                dNME_region=cutoff_ts_only_out_ts$dNME_cutoff[[as.character(cutoff_ts_num_output$dNME[[x]])]]
                                                dMML_region=cutoff_ts_only_out_ts$dMML_cutoff[[as.character(cutoff_ts_num_output$dMML[[x]])]]
                                                UC_region=cutoff_ts_only_out_ts$UC_cutoff[[as.character(cutoff_ts_num_output$UC[[x]])]]
                                                #if(length(UC_region)>0){
                                                UC_region_rand=sample(region_rand_pool[[ts]],length(UC_region))
                                                return(data.table(        
                                                    dNME_specific=length(setdiff(dNME_region,UC_region)),
                                                    dMML_specific=length(setdiff(dMML_region,UC_region)),
                                                    UC_specific_dNME=length(setdiff(UC_region,dNME_region)),
                                                    UC_specific_dMML=length(setdiff(UC_region,dMML_region)),
                                                    shared_dNME=length(intersect(dNME_region,UC_region)),
                                                    shared_dMML=length(intersect(dMML_region,UC_region)),
                                                    dNME_specific_rand=length(setdiff(dNME_region,UC_region_rand)),
                                                    dMML_specific_rand=length(setdiff(dMML_region,UC_region_rand)),
                                                    UC_specific_dNME_rand=length(setdiff(UC_region_rand,dNME_region)),
                                                    UC_specific_dMML_rand=length(setdiff(UC_region_rand,dMML_region)),
                                                    shared_dNME_rand=length(intersect(dNME_region,UC_region_rand)),
                                                    shared_dMML_rand=length(intersect(dMML_region,UC_region_rand)),
                                                    dNME_cutoff=cutoff_ts_num_output$dNME[[x]],
                                                    dMML_cutoff=cutoff_ts_num_output$dMML[[x]],
                                                    UC_cutoff=cutoff_ts_num_output$UC[[x]]
                                                )
                                                )
                                                },
                                                cutoff_ts_only_out_ts=cutoff_ts_only_out_ts,cutoff_ts_num_output=cutoff_ts_num_output,mc.cores=20)))
  
    cat("Finish processing:",ts,'in',proc.time()[[3]]-tt1,'\n')

}
dev.off()
saveRDS(cutoff_ts_num_output_all,'../downstream/output/mouse_analysis/UC_dNME_olap/regions_ts_only_tb.rds')
cutoff_ts_num_output_all=list()

#For not one and onely one
dNME_cutoff=list()
dMML_cutoff=list()
UC_cutoff=list()
step=0.001
for(ts in names(dNME_in)){
    cat("Processing:",ts,"\n")
    dNME_cutoff[[ts]]=cutoff_selection_only(dNME_in,ts,ts_only=F,step_in=step)
    
    dMML_cutoff[[ts]]=cutoff_selection_only(dMML_in,ts,ts_only=F,step_in=step)
    UC_cutoff[[ts]]=cutoff_selection_only(UC_in_only,ts,ts_only=F,step_in=step)


}

saveRDS(list(UC_cutoff=lapply(UC_cutoff,assign_name,step_in=step),
            dMML_cutoff=lapply(dMML_cutoff,assign_name,step_in=step),
            dNME_cutoff=lapply(dNME_cutoff,assign_name,step_in=step)
            ),'../downstream/output/mouse_analysis/UC_dNME_olap/regions_non_ts_only_01.rds')


cutoff_ts_all_out=readRDS('../downstream/output/mouse_analysis/UC_dNME_olap/regions_non_ts_only_001.rds')

cutoff_ts_all_out=list(UC_cutoff=lapply(UC_cutoff,assign_name,step_in=step),
            dMML_cutoff=lapply(dMML_cutoff,assign_name,step_in=step),
            dNME_cutoff=lapply(dNME_cutoff,assign_name,step_in=step))
cutoff_ts_num_output_all=list()
set.seed(12345)
pdf(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/UC_dNME_dMML_cutoff_num_non_ts_only_0001.pdf'))
for(ts in names(region_rand_pool)){
    
    tt1=proc.time()[[3]]
    
    cutoff_ts_all_out_ts=lapply(cutoff_ts_all_out,function(x) x[[ts]])
    ts_num=lapply(cutoff_ts_all_out_ts,function(x) unlist(lapply(x,length)))
    ts_num_dt=data.table(cutoff=names(ts_num$UC))
    ts_num_dt$dNME=ts_num$dNME[ts_num_dt$cutoff]
    ts_num_dt$dMML=ts_num$dMML[ts_num_dt$cutoff]
    ts_num_dt$UC=ts_num$UC[ts_num_dt$cutoff]
    ts_num_dt_mt=melt.data.table(ts_num_dt,id.vars='cutoff',variable.name='stat_type',value.name='region_count')
   ts_num_dt_mt$cutoff=as.numeric(ts_num_dt_mt$cutoff)
    print(ggplot(ts_num_dt_mt,aes(x=cutoff,y=region_count,color=stat_type))+geom_line()+
            ylab("Number of regions")+ggtitle(ts))
    cat("Processing:",ts,'\n')
    #Expand grid based on number of regions
    # cutoff_ts_num_output=as.data.table(expand.grid(seq(0,max(ts_num_dt[UC>=round(mean(ts_num_dt$UC)/100,digits=2)]$cutoff),step),
    #                                     seq(0,max(ts_num_dt[dNME>=round(mean(ts_num_dt$dNME)/100,digits=2)]$cutoff),step*5),
    #                                     seq(0,max(ts_num_dt[dMML>=round(mean(ts_num_dt$dMML)/100,digits=2)]$cutoff),step*5)))
    num=100
    num_dMML_dNME=50
    cutoff_ts_num_output=as.data.table(expand.grid( ts_num_dt[unlist(lapply( quantile(ts_num_dt$UC,seq(0,1,1/num)),function(x) which.min(abs(ts_num_dt$UC-x))))]$cutoff,
                                                    ts_num_dt[unlist(lapply( quantile(ts_num_dt$dNME,seq(0,1,1/num_dMML_dNME)),function(x) which.min(abs(ts_num_dt$dNME-x))))]$cutoff,
                                                    ts_num_dt[unlist(lapply( quantile(ts_num_dt$dMML,seq(0,1,1/num_dMML_dNME)),function(x) which.min(abs(ts_num_dt$dMML-x))))]$cutoff
                                                        ))
    colnames(cutoff_ts_num_output)=c("UC","dNME","dMML")

    cutoff_ts_num_output_all[[ts]]=
    cbind(cutoff_ts_num_output,
        fastDoCall("rbind",
                    mclapply(1:nrow(cutoff_ts_num_output),
                                function(x,cutoff_ts_all_out_ts,cutoff_ts_num_output) {
                                                dNME_region=cutoff_ts_all_out_ts$dNME_cutoff[[as.character(cutoff_ts_num_output$dNME[[x]])]]
                                                dMML_region=cutoff_ts_all_out_ts$dMML_cutoff[[as.character(cutoff_ts_num_output$dMML[[x]])]]
                                                UC_region=cutoff_ts_all_out_ts$UC_cutoff[[as.character(cutoff_ts_num_output$UC[[x]])]]
                                                #if(length(UC_region)>0){
                                                UC_region_rand=sample(region_rand_pool[[ts]],length(UC_region))
                                                dMML_region_rand=sample(region_rand_pool[[ts]],length(dMML_region))
                                                dNME_region_rand=sample(region_rand_pool[[ts]],length(dNME_region))
                                                return(data.table(        
                                                    dNME_specific=length(setdiff(dNME_region,UC_region)),
                                                    dMML_specific=length(setdiff(dMML_region,UC_region)),
                                                    UC_specific_dNME=length(setdiff(UC_region,dNME_region)),
                                                    UC_specific_dMML=length(setdiff(UC_region,dMML_region)),
                                                    shared_dNME=length(intersect(dNME_region,UC_region)),
                                                    shared_dMML=length(intersect(dMML_region,UC_region)),
                                                    UC_rand_dNME_specific=length(setdiff(dNME_region,UC_region_rand)),
                                                    UC_rand_dMML_specific=length(setdiff(dMML_region,UC_region_rand)),
                                                    UC_rand_UC_specific_dNME=length(setdiff(UC_region_rand,dNME_region)),
                                                    UC_rand_UC_specific_dMML=length(setdiff(UC_region_rand,dMML_region)),
                                                    UC_rand_shared_dNME=length(intersect(dNME_region,UC_region_rand)),
                                                    UC_rand_shared_dMML=length(intersect(dMML_region,UC_region_rand)),

                                                    dNME_rand_dNME_specific=length(setdiff(dNME_region_rand,UC_region)),
                                                    dNME_rand_UC_specific_dNME=length(setdiff(UC_region,dNME_region_rand)),
                                                    dNME_rand_shared_dNME=length(intersect(dNME_region_rand,UC_region)),

                                                    dMML_rand_dMML_specific=length(setdiff(dMML_region_rand,UC_region)),
                                                    dMML_rand_UC_specific_dMML=length(setdiff(UC_region,dMML_region_rand)),
                                                    dMML_rand_shared_dMML=length(intersect(dMML_region_rand,UC_region)),

                                                    dNME_cutoff=cutoff_ts_num_output$dNME[[x]],
                                                    dMML_cutoff=cutoff_ts_num_output$dMML[[x]],
                                                    UC_cutoff=cutoff_ts_num_output$UC[[x]]
                                                )
                                                )
                                                },
                                                cutoff_ts_all_out_ts=cutoff_ts_all_out_ts,cutoff_ts_num_output=cutoff_ts_num_output,mc.cores=20)))
  
    cat("Finish processing:",ts,'in',proc.time()[[3]]-tt1,'\n')

}
dev.off()

saveRDS(cutoff_ts_num_output_all,'../downstream/output/mouse_analysis/UC_dNME_olap/regions_non_ts_only_tb_001.rds')
#Generate figures for each cutoff
cutoff_ts_num_output_all=readRDS('../downstream/output/mouse_analysis/UC_dNME_olap/regions_non_ts_only_tb_001.rds')
# nolint end