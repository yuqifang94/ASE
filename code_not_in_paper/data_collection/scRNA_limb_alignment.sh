#Should be in mate1
#Upon dissolution of the Single Cell 3’ Gel Bead in a GEM, primers containing (i) an Illumina R1 sequence (read
# 1 sequencing primer), (ii) a 16 bp 10x Barcode, (iii) a 10 bp Unique Molecular Identifier (UMI) and (iv) a polydT primer sequence are released and mixed with cell lysate and Master Mix. Incubation of the GEMs then
# produces barcoded, full-length cDNA from poly-adenylated mRNA. After incubation, the GEMs are broken
# and the pooled fractions are recovered.
#https://github.com/alexdobin/STAR/issues/768
#STAR --runThreadN 20 --genomeDir /users/yfang/yfang_dcs04/referenceGenome/mm10_10x_STAR/ --soloBarcodeMate 1 --clip3pNbases 39 0  --readFilesIn E10_5_ENCSR874BOF/FT-SA17497_S1_L004_R1_001.fastq.gz E10_5_ENCSR874BOF/FT-SA17497_S1_L004_R2_001.fastq.gz --soloType CB_UMI_Simple --soloCBstart 1   --soloCBlen 16   --soloUMIstart 17   --soloUMIlen 10 --readFilesCommand zcat --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloCBwhitelist /dcs04/feinberg/data/personal/yfang/ASE_clean_run_403/downstream/data/mouse_scRNA/scFetalLimb/10xV2WhiteList.txt 1>logFiles/E10_5_ENCSR874BOF.out 2>logFiles/E10_5_ENCSR874BOF.err

CPUS=20
REF=/users/yfang/yfang_dcs04/referenceGenome/mm10_10x_STAR/
BC=/dcs04/feinberg/data/personal/yfang/ASE_clean_run_403/downstream/data/mouse_scRNA/scFetalLimb/10xV2WhiteList.txt
for sample in E1*
do
     
     R2=$sample/*_R2_*fastq.gz
     R1=$sample/*_R1_*fastq.gz
     mkdir $sample/output
     STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R1 $R2 --runDirPerm All_RWX --readFilesCommand zcat \
          --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloUMIlen 10 --soloStrand Forward \
          --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
          --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
          --outSAMtype BAM SortedByCoordinate\
          --soloFeatures Gene GeneFull Velocyto --soloOutFileNames $sample/output/ genes.tsv barcodes.tsv matrix.mtx
done
sample=E10_5_ENCSR874BOF
R2=$sample/*_R2_*fastq.gz
R1=$sample/*_R1_*fastq.gz
[! -d "/path/to/dir" ] && mkdir $sample/output
STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R1 $R2 --runDirPerm All_RWX --readFilesCommand zcat\
     --soloBarcodeMate 1   --clip5pNbases 39 0 --soloStrand Forward \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \ 
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --outFilterScoreMin 30 \
     --outSAMtype BAM SortedByCoordinate\
     --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ genes.tsv barcodes.tsv matrix.mtxq
#V1 single end
#V2 PE Forward