if (!requireNamespace("MotifDb", quietly = TRUE))
{BiocManager::install("MotifDb")}
library(MotifDb)
#need libatlas for linux: ml atlas
if (!requireNamespace("motifbreakR", quietly = TRUE))
{BiocManager::install("motifbreakR")}
library(motifbreakR)
if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
{BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")}
library("BSgenome.Mmusculus.UCSC.mm10")
if (!requireNamespace("BiocParallel", quietly = TRUE))
{BiocManager::install("BiocParallel")}
library(BiocParallel)
library('org.Mm.eg.db')
library(VGAM)
if (!requireNamespace("TFBSTools", quietly = TRUE))
{ devtools::install_github("ge11232002/TFBSTools")}
library(TFBSTools)
file_in="SNP_DBA.rds"
method_motif="default"
cat('Reading in',file_in,'for: ',n,'with gap of',1/cut_number,'\n')
tt1=proc.time()[[3]]
variant_in=readRDS(paste('../downstream/output/',file_in,sep=''))
names(variant_in)=paste0("SNP",1:length(variant_in))
JASPAR_2020=c(readRDS('../downstream/output/JASPAR_2020_mouse_motif_PFM.rds'),
              readRDS('../downstream/input/JASPAR_2020_human_PFM.rds'))

meta_JASPAR=data.frame(providerName=names(JASPAR_2020),providerId=unlist(lapply(JASPAR_2020@listData,function(x) x@ID)),
                       datasource='JASPAR2020',geneSymbol=unlist(lapply(JASPAR_2020@listData,function(x) x@name)),stringsAsFactors = F)
#meta_JASPAR$geneId=mapIds(org.Hs.eg.db,meta_JASPAR$geneSymbol,'ENTREZID','SYMBOL')
meta_JASPAR$geneIdType=NA
meta_JASPAR$geneIdType[!is.na(meta_JASPAR$geneId)]='JASPAR'
 meta_JASPAR$organism='MusmucusHsapiens'
meta_JASPAR[,c('protein','proteinIdType','sequenceCount','bingdingSequence','bindingDomain','tfFamily','experiment','pubmedID')] =NA
mcols(JASPAR_2020)=meta_JASPAR
JASPAR_2020@listData=lapply(JASPAR_2020@listData,function(x) x@profileMatrix)
#function to separate them into  non-overlapping ?: split
tt2=proc.time()[[3]]
cat('Finishing reading in data in: ',tt2-tt1,'\n')
cat('Starting motifbreak\n')
pwmList=JASPAR_2020
snpList = variant_in
filterp = TRUE
threshold = 1e-4
method =method_motif
bkg = c(A=0.25, C=0.25, G=0.25, T=0.25)
verbose=T
BPPARAM =MulticoreParam(workers=20,progressbar = TRUE)
if (.Platform$OS.type == "windows" && inherits(BPPARAM,
                                                 "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n",
                   "Windows, please supply an alternative BPPARAM"))
  }
  cores <- bpnworkers(BPPARAM)
  num.snps <- length(snpList)
  if (num.snps < cores) {
    cores <- num.snps
  }
  if (is(BPPARAM, "SnowParam")) {
    bpstart(BPPARAM)
    cl <- bpbackend(BPPARAM)
    clusterEvalQ(cl, library("MotifDb"))
  }
  genome.package <- attributes(snpList)$genome.package
  if (requireNamespace(eval(genome.package), quietly = TRUE,
                       character.only = TRUE)) {
    genome.bsgenome <- eval(parse(text = paste(genome.package,
                                               genome.package, sep = "::")))
  }else {
    stop(paste0(eval(genome.package), " is the genome selected for this snp list and \n",
                "is not present on your environment. Please load it and try again."))
  }
pwms <- motifbreakR:::preparePWM(pwmList = pwmList, filterp = filterp,
                     scoreThresh = threshold, bkg = bkg, method = method)
saveRDS(pwms,'../downstream/output/mm10_JASPR_pwm.rds')

