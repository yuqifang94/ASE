library(data.table)
library(rtracklayer)
library(stringr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Gmisc)
window_size=100
cluster_gap=200
sub_pattern<-function(pattern,ref){
    gr_pat=vmatchPattern(pattern,ref)
    fastDoCall('c',lapply(names(gr_pat)[unlist(lapply(gr_pat,length))>0],function(x) GRanges(seqnames=x, ranges=gr_pat[[x]])))
}
sub_hetCpG<-function(pattern,ref,vcf){
    mm10_hetCG=sub_pattern(pattern,ref)
    mm10_hetCG$hetCG_type=pattern
    fix=ifelse(pattern=='CN','end','start')
    olap_CG=findOverlaps(resize(mm10_hetCG,1,fix=fix),vcf)
    mm10_hetCG$refCG=NA
     if (pattern=="CN") {mm10_hetCG$refCG[queryHits(olap_CG)]=paste0("C",vcf$REF[subjectHits(olap_CG)])}
    else if(pattern=="NG") {mm10_hetCG$refCG[queryHits(olap_CG)]=paste0(vcf$REF[subjectHits(olap_CG)],"G")}
    mm10_hetCG$altCG=NA
    if (pattern=="CN") {mm10_hetCG$altCG[queryHits(olap_CG)]=paste0("C",vcf$ALT[subjectHits(olap_CG)])}
    else if(pattern=="NG") {mm10_hetCG$altCG[queryHits(olap_CG)]=paste0(vcf$ALT[subjectHits(olap_CG)],"G")}
    mm10_hetCG[which(mm10_hetCG$refCG=="CG"|mm10_hetCG$altCG=="CG")]
}
split_by_split<-function(row_in){
    if(row_in$split!=1){
        return(data.table(seqnames=row_in$seqnames,
        start=(row_in$start)+row_in$split_wd*(0:(row_in$split-1)),
        end=row_in$start+row_in$split_wd*(1:row_in$split)-1,
        width=row_in$split_wd-1,strand="*"))
    }else(return(row_in))
}
CpG_N_assign<-function(gr_in,CG,het_CG){
    olap_CG=findOverlaps(gr_in,CG)
    CG_DT=data.table(query=queryHits(olap_CG),start_pos=start(CG)[subjectHits(olap_CG)])
    CG_DT=CG_DT[which(start_pos %in% setdiff(start_pos,start(het_CG[het_CG$refCG=="CG"]))),list(N=length(start_pos),CpG=paste0("[",paste(start_pos,collapse=", "),"]")),by=query]
    CG_DT=CG_DT[order(CG_DT$query,decreasing=F)]
    gr_in=gr_in[CG_DT$query]
    gr_in$N=NA
    gr_in$N=CG_DT$N
    gr_in$CpGs=NA
    gr_in$CpGs=CG_DT$CpG
    return(gr_in)
}

vcf_in=fread('all_SNPs_DBA_2J_GRCm38.txt',header=F,stringsAsFactors=F)
vcf_in=vcf_in[,-"V4"]
colnames(vcf_in)=c("ID","chr","start","SNP")
vcf_in$REF=sub('/.*','',vcf_in$SNP)
vcf_in$ALT=sub('.*/','',vcf_in$SNP)
vcf_in=vcf_in[,-"SNP"]
vcf_in=makeGRangesFromDataFrame(vcf_in,end.field="start",keep.extra.columns=T)
#getting hetCpG
vcf_in$tri= as.character(Views(Mmusculus,vcf_in+1))
vcf_in$Het=grepl('CG', vcf_in$tri)
vcf_in=vcf_in+window_size
#Psedu Phasing
vcf_in_red=reduce(vcf_in)
#expand windows
#win_exp is expand the windows, default setting: 100
#read in fasta file
ref_genome=readDNAStringSet('../fasta/BL6DBA.masked.fa')
ref_genome=ref_genome[names(ref_genome)!="chrM"]
mm10_CG=sub_pattern("CG",ref_genome)
#Find the heterozygous CpG in ref or alt
mm10_hetCpG=c(sub_hetCpG("CN",ref_genome,vcf_in),sub_hetCpG("NG",ref_genome,vcf_in))
#Count overlap in CG
#subset the regions with N>=20
N_cutoff=20

vcf_in_red$N=countOverlaps(vcf_in_red,mm10_CG)
vcf_in_red_N_lrg=as.data.table(vcf_in_red[vcf_in_red$N>N_cutoff])
vcf_in_red=vcf_in_red[vcf_in_red$N<=N_cutoff]
vcf_in_red_N_lrg$split=ceiling(vcf_in_red_N_lrg$N/N_cutoff)
vcf_in_red_N_lrg$split_wd=round(vcf_in_red_N_lrg$width/(vcf_in_red_N_lrg$split))
#splitting function for each row

while(any(vcf_in_red_N_lrg$N>20)){

    vcf_in_red_N_lrg=fastDoCall("rbind",lapply(1:nrow(vcf_in_red_N_lrg),function(x) split_by_split(vcf_in_red_N_lrg[x])))
    vcf_in_red_N_lrg$N=countOverlaps(makeGRangesFromDataFrame(vcf_in_red_N_lrg),mm10_CG)
    vcf_in_red_N_lrg$split=ceiling(vcf_in_red_N_lrg$N/N_cutoff)
    vcf_in_red_N_lrg$split_wd=round(vcf_in_red_N_lrg$width/(vcf_in_red_N_lrg$split))
    vcf_in_red_N_lrg=vcf_in_red_N_lrg[vcf_in_red_N_lrg$N>0]
    vcf_in_red=c(vcf_in_red,makeGRangesFromDataFrame(vcf_in_red_N_lrg[vcf_in_red_N_lrg$N<=N_cutoff],keep.extra.columns=T))
    vcf_in_red_N_lrg=vcf_in_red_N_lrg[vcf_in_red_N_lrg$N>N_cutoff]
}
vcf_in_red$split=NULL
vcf_in_red$split_wd=NULL
#getting CpG and N 
vcf_in_red=CpG_N_assign(vcf_in_red,mm10_CG,mm10_hetCpG)
#Count overlapping CG for g1 and g2
olap_hetCG=findOverlaps(vcf_in_red,mm10_hetCpG)
hetCG_DT=data.table(query=queryHits(olap_hetCG),start_pos=start(mm10_hetCpG)[subjectHits(olap_hetCG)],
    refCG=mm10_hetCpG$refCG[subjectHits(olap_hetCG)],altCG=mm10_hetCpG$altCG[subjectHits(olap_hetCG)])
hetCGg1_DT=hetCG_DT[refCG=="CG",list(hetCpGg1=paste0("[",paste(start_pos,collapse=", "),"]")),by=query]
hetCGg1_DT=hetCGg1_DT[order(hetCGg1_DT$query,decreasing=F)]
hetCGg2_DT=hetCG_DT[altCG=="CG",list(hetCpGg2=paste0("[",paste(start_pos,collapse=", "),"]")),by=query]
hetCGg2_DT=hetCGg2_DT[order(hetCGg2_DT$query,decreasing=F)]
seqlengths(vcf_in_red)= seqlengths(ref_genome)[match(names(seqlengths(vcf_in_red)),names(seqlengths(ref_genome)))]
genome(vcf_in_red)="mm10"
vcf_in_red$hetCpGg1=NA
vcf_in_red$hetCpGg1[hetCGg1_DT$query]=hetCGg1_DT$hetCpGg1
vcf_in_red$hetCpGg2=NA
vcf_in_red$hetCpGg2[hetCGg2_DT$query]=hetCGg2_DT$hetCpGg2
export.gff3(vcf_in_red,'BL6DBA_het.cpelasm.gff')
#homozygous region
vcf_in_hom=gaps(vcf_in_red)
vcf_in_hom=vcf_in_hom[strand(vcf_in_hom)=="*"]
#getting CpG and N 
vcf_in_hom=CpG_N_assign(vcf_in_hom,mm10_CG,mm10_hetCpG)
vcf_in_hom$score=vcf_in_hom$N
vcf_in_hom$N=NULL
export.gff3(vcf_in_hom,'BL6DBA_hom.cpelasm.gff')
sed -i 's/%2c/,/g' BL6DBA_het.cpelasm.gff
sed -i 's/rtracklayer/\./g' BL6DBA_het.cpelasm.gff
sed -i 's/sequence_feature/\./g' BL6DBA_het.cpelasm.gff
sed -i 's/];$/]/g' BL6DBA_het.cpelasm.gff
sed -i 's/%2c/,/g' BL6DBA_hom.cpelasm.gff
sed -i 's/rtracklayer/\./g' BL6DBA_hom.cpelasm.gff
sed -i 's/sequence_feature/\./g' BL6DBA_hom.cpelasm.gff
sed -i 's/];$/]/g' BL6DBA_hom.cpelasm.gff


#use gaps to get the ranges between heterozygous regions
# vcf_in_red$CG=str_count(as.character(Views(Mmusculus,vcf_in_red)),'CG')
# vcf_in_red=vcf_in_red[vcf_in_red$CG>0]
# olap=findOverlaps(vcf_in,vcf_in_red)
# vcf_in$hap=subjectHits(olap)
# vcf_in$hap_size=vcf_in_red$hap_size[subjectHits(olap)]
# vcf_in=as.data.table(vcf_in)
# vcf_in=vcf_in[,-c("end",'width','strand','tri')]
# colnames(vcf_in)=c("CHROM","POS","ID",'REF',"ALT","HET","HAP","Hap-size")
# vcf_in$QUAL=100
# vcf_in$FILTER="PASS"
# vcf_in$INFO="HET=FALSE"
# vcf_in$INFO[vcf_in$HET]="HET=TRUE"
# vcf_in=vcf_in[,-"HET"]
# vcf_in$FORMAT="GT:FI:PS"
# vcf_in$FORMAT[vcf_in$Hap]="GT:FI"
# vcf_in$BL6DBA=paste0("0|1:1:",vcf_in$HAP)
# vcf_in$BL6DBA[vcf_in$Hap]="0/1:1"
# vcf_in=vcf_in[,-c("HAP","Hap-size")]
# write.table(vcf_in,file="BL6DBA.phased.noheader.vcf",sep='\t',row.names=F,col.names=T,quote=FALSE)
# #make gff file




