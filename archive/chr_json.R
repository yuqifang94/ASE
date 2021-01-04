chr_json<-function(id_in,json_in){
  chr=json_in$Chromosome[id_in]
  pos=json_in$Position[id_in]
  REF=json_in$`Reference Allele`[id_in]
  ALT=json_in$`Alternative Allele`[id_in]
  enc=json_in$`Is In Enhancer Region`[id_in]
  GWAS=json_in$`Is Near GWAS Variant`[id_in]
  HetCpG=json_in$`Is On Heterogenous CpG`[id_in]
  promoter=json_in$`Is In Promoter Region`[id_in]
  freq1000=json_in$`1000 Genomes Allele Frequency`[id_in]
  Anc_allele=json_in$`Ancestral Allele`[id_in]
  TS_analysis=json_in$`Tissue Specific Analysis`[[id_in]]
  rm(json_in)
  gc()
  if (!is.null(TS_analysis)){
    TS_analysis=as.data.frame(TS_analysis)
    TS_analysis=TS_analysis[which(TS_analysis$Experiment$`rdfs:label`=='Bisulfite-Seq'),]
    tissue=TS_analysis$Tissue$`foaf:name`
    subject=TS_analysis$Patient
    pvalue=TS_analysis$`Raw P-value`
    qvalue=TS_analysis$`FDR P-value`
    REF_unmet=TS_analysis$`Ref Allele Unmethylated CpG Count`
    REF_met=TS_analysis$`Ref Allele Methylated CpG Count`
    ALT_unmet=TS_analysis$`Alt Allele Unmethylated CpG Count`
    ALT_met=TS_analysis$`Alt Allele Methylated CpG Count`
    nsubj=length(ALT_met)
    df_out=data.frame(chr=rep(chr,nsubj),start=rep(pos,nsubj),end=rep(pos,nsubj),REF=rep(REF,nsubj),ALT=rep(ALT,nsubj),
                      enc=rep(enc,nsubj),GWAS=rep(GWAS,nsubj),HetCpG=rep(HetCpG,nsubj),promoter=rep(promoter,nsubj),freq1000=rep(freq1000,nsubj),
                      Anc_allele=rep(Anc_allele,nsubj),tissue=tissue,subject=subject,pvalule=pvalue,qvalue=qvalue,REF_unmet=REF_unmet,REF_met=REF_met,
                      ALT_unmet=ALT_unmet,ALT_met=ALT_met,stringsAsFactors = F)
  }else{
    
    df_out=data.frame(chr=chr,start=pos,end=pos,REF=REF,ALT=ALT,
                      enc=enc,GWAS=GWAS,HetCpG=HetCpG,promoter=promoter,freq1000=freq1000,
                      Anc_allele=Anc_allele,tissue="No tissue specific",subject="No tissue specific",pvalule="No tissue specific",
                      qvalue="No tissue specific",REF_unmet="No tissue specific",REF_met="No tissue specific",
                      ALT_unmet="No tissue specific",ALT_met="No tissue specific",stringsAsFactors = F)
  }
  return(df_out)
}

