library(SRAdb)
sqlfile <- file.path(system.file('D:/', package='SRAdb'), 'SRAmetadb.sqlite')
sra_con <- dbConnect(SQLite(),sqlfile)
scRNA_table=as.data.frame(read.csv('../../Meeting ppt/H1 data/scRNA.csv',header = T,stringsAsFactors = F))
scRNA_table$title2=unlist(lapply(strsplit(scRNA_table$Title,'\\.'),function(x) x[[1]][1]))
unique(scRNA_table$title2)
SRX=scRNA_table[scRNA_table$title2=='single cell H1 hESC cell [H1_Exp1','SRA.Accession']
run_SRR<- sraConvert( SRX, sra_con = sra_con )$run
