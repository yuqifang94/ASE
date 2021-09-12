source('mainFunctions_sub.R')


plot_out_cor=list()
for(cutoff_char in  c("005","02","03","04")){
    output_dir=paste0('../downstream/output/mouse_analysis/cutoff_testing/',cutoff_char,'/')
    plot_out_all=paste0(output_dir,'tissue_out_N17_kmeans_10run_filtered_all_region_',cutoff_char,'.rds')
   plot_cutoff=plot_correlation(readRDS(plot_out_all),pdf_fn=paste0(output_dir,cutoff_char,'_correlation_main.pdf'),plot_pdf=F)+
                ggtitle(paste0("cutoff:",cutoff_char))+theme(plot.title = element_text(size=36,hjust=0.5))
    plot_out_cor[[cutoff_char]]=plot_cutoff

                
}
ggsave(
    filename='../downstream/output/mouse_analysis/cutoff_testing/correlation_fig.pdf',
    plot=ggarrange(plotlist=plot_out_cor, ncol=2, nrow=2, common.legend = TRUE, legend="bottom"),
   width = 21,
   height =28,
)

#heatmap for values
d=readRDS(UC_merge_max_loc_file)
for(cutoff_char in  c("005","01","02","03","04")){3
    output_dir=paste0('../downstream/output/mouse_analysis/cutoff_testing/',cutoff_char,'/')
    cluster_assigned_dir=paste0(output_dir,'cluster_assigned/')
    figure_name=paste0(figure_path,'all_sc_N17_ft_kmeans_10run_filtered_all',cutoff_char,'.png')
    clu=readRDS(paste0(cluster_assigned_dir,'cluster_assginment_filtered_',cutoff_char,'.rds'))
    plot_heatmap_cluster(d,clu,figure_name,figure_width=1800,figure_height=2000,res=200)
}

