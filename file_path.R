# file locations ----------------------------------------------------------

#global pval cutoff file names
pval_cutoff=0.1 #Onuchic use 0.1
gff_in_file='../downstream/input/human_analysis/gff_in.rds'
variant_HetCpG_file='../downstream/input/human_analysis/variant_HetCpG_new.rds'
GR_file='../downstream/output/human_analysis/CPEL_outputs/GRs_final1.rds'
GR_allele_file='../downstream/output/human_analysis/CPEL_outputs/GRs_allele_final1.rds'
hetCpG_gff_file='../downstream/input/hetCpG_gff_final1.rds'
GR_merge_file="../downstream/output/human_analysis/CPEL_outputs/GR_merge_final12_ls.rds"
CpG_hg19_file='../downstream/input/human_analysis/CpG_hg19.rds'
variant_HetCpG_meta_file='../downstream/output/human_analysis/CPEL_outputs/variant_HetCpG_meta_final1_ls.rds'
genomic_features_file="../downstream/input/human_analysis/genomic_features2020.rds"
NME_agnostic_file="../downstream/output/human_analysis/CPEL_outputs/NME_allele_agnostic_merge_20k_homogeneous_excluding_dMML2.rds"
MML_agnostic_file="../downstream/output/human_analysis/CPEL_outputs/MML_allele_agnostic_merge_20k_homogeneous2.rds"
motif_gene_file='../downstream/output/human_analysis/motif_analysis/motif_all_JASPAR_default.rds' #For all SNP
NME_agnostic_DNase_file="../downstream/output/human_analysis/CPEL_outputs/allele_agnostic_hg19_DNase_NME_homogeneous_excluding_dMML.rds"
MML_agnostic_DNase_file="../downstream/output/human_analysis/CPEL_outputs/allele_agnostic_hg19_DNase_MML_homogeneous.rds"
NME_agnostic_ASM_file="../downstream/output/human_analysis/CPEL_outputs/NME_agnostic_ASM.rds"
MML_agnostic_ASM_file="../downstream/output/human_analysis/CPEL_outputs/MML_agnostic_ASM.rds"
NME_agnostic_comp_file="../downstream/output/human_analysis/CPEL_outputs/NME_agnostic_comp.rds"
MML_agnostic_comp_file="../downstream/output/human_analysis/CPEL_outputs/MML_agnostic_comp.rds"
NME_agnostic_all_file="../downstream/output/human_analysis/CPEL_outputs/NME_agnostic_all.rds"
MML_agnostic_all_file="../downstream/output/human_analysis/CPEL_outputs/MML_agnostic_all.rds"
variant_HetCpG_meta_dt_file='../downstream/output/human_analysis/CpG_density/variant_HetCpG_meta_dt.rds'
variant_HetCpG_meta_dt_uq_file='../downstream/output/human_analysis/CpG_density/variant_HetCpG_meta_dt_relative_dNME2uq.rds'
figure_path='../downstream/output/graphs_tables/'
#This is from Ken
DNase_hg19_file='../downstream/input/human_analysis/DNase_hg19_250bp.rds'
control_hg19_file='../downstream/input/human_analysis/DNase_hg19_250bp_control.rds'
JASPAR_motif_hg19_file='../downstream/output/human_analysis/motif_analysis/motif_JASPAR_hg19.rds'
#This is from Ken
DNase_mm10_file='../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp.rds'
control_mm10_file='../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp_control.rds'
JASPAR_motif_mm10_file='../downstream/input/mouse_analysis/motif_analysis/motif_JASPAR_mm10.rds'
#Mouse analysis
MML_matrix_file='../downstream/output/mouse_analysis/CPEL_outputs/MML_matrix_mouse_all_dedup_N2_all_regions.rds'
NME_matrix_file='../downstream/output/mouse_analysis/CPEL_outputs/NME_matrix_mouse_all_dedup_N2_all_regions.rds'
CG_density_mouse='../downstream/output/mouse_analysis/CpG_density_NME/region_CG_exp.rds'
mouse_DNase_control_gff_file='../downstream/output/mouse_analysis/CPEL_inputs/mm10_allele_agnostic_analysis_DNase_control.gff'
mouse_compliment_gff_file='../downstream/output/mouse_analysis/CPEL_inputs/mm10_allele_agnostic_analysis_compliment.gff'
UC_in_matrix_ls_file='../downstream/output/mouse_analysis/CPEL_outputs/UC_matrix_ls_N2_all_regions.rds'
UC_in_MDS_all_file='../downstream/output/mouse_analysis/CPEL_outputs/UC_MDS_N2_all_regions.rds'
mouse_enhancer_bin='../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds'
dir_cluster_in_01='../downstream/input/mouse_analysis/clustering/tissue_specific/uc_01/'
UC_merge_file='../downstream/output/mouse_analysis/CPEL_outputs/UC_merge_all_regions.rds'
UC_merge_max_loc_file='../downstream/output/mouse_analysis/CPEL_outputs/UC_merge_all_regions_max_loc.rds'
UC_merge_max_loc_01_file='../downstream/output/mouse_analysis/CPEL_outputs/UC_merge_all_regions_max_loc_01.rds'
dir_out_cluster01='../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_1/cluster_assigned/'
cluster_01_region_out_fn='../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_1/cluster_all_region_assignment_filtered_0_1.rds'
dmml_cor_file='../downstream/input/mouse_analysis/correlation_analysis/all_regions/fulldmmlcor.rds'
dnme_cor_file='../downstream/input/mouse_analysis/correlation_analysis/all_regions/fulldnmecor.rds'
dir_out_rds_correlation='../downstream/output/mouse_analysis/correlation/'
bin_enhancer_rds='../downstream/input/mouse_analysis/enhancer_selection/bin_enhancer.rds'
bin_enhancer_bed='../downstream/input/mouse_analysis/enhancer_selection/bin_enhancer.bed'
GO_01_dir='../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/'
