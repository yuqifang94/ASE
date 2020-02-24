##############################Look at CGT, CGA, CGG, CGC -> CTA, CTT, CTC, CTG ################################################## 
#Try percent ASM
variant_values <- readRDS("E:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/motif/variant_values.rds")
#Ref=CG ALT=CT: CI:3.073674 3.205317
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('CGT','CGA','CGG','CGC')& 
                          variant_values$nucleo_3_alt %in% c('CTT','CTA','CTG','CTC')) |
                         (variant_values$nucleo_3_alt %in% c('CGT','CGA','CGG','CGC')& 
                            variant_values$nucleo_3_ref %in% c('CTT','CTA','CTG','CTC'))]=TRUE

variant_enrich(variant_values)


##############################Look at CGT, CGA, CGG, CGC -> CAA, CAT, CAC, CAG ################################################## 
#Try percent ASM
#Ref=CG ALT=CA: CI:3.263671 3.344500
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('CGT','CGA','CGG','CGC')& 
                          variant_values$nucleo_3_alt %in% c('CAT','CAA','CAG','CAC')) |
                         (variant_values$nucleo_3_alt %in% c('CGT','CGA','CGG','CGC')& 
                            variant_values$nucleo_3_ref %in% c('CAT','CAA','CAG','CAC'))]=TRUE

variant_enrich(variant_values)

##############################Look at CGT, CGA, CGG, CGC -> CCA, CCT, CCC, CCG ################################################## 
#Try percent ASM
#Ref=CG ALT=CC: CI:2.815779 2.940513
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('CGT','CGA','CGG','CGC')& 
                          variant_values$nucleo_3_alt %in% c('CCT','CCA','CCG','CCC')) |
                         (variant_values$nucleo_3_alt %in% c('CGT','CGA','CGG','CGC')& 
                            variant_values$nucleo_3_ref %in% c('CCT','CCA','CCG','CCC'))]=TRUE

variant_enrich(variant_values)

##############################Look at GCT, GCA, GCG, GCC -> GAA, GAT, GAC, GAG ################################################## 
#Try percent ASM
#Ref=GC ALT=GA: CI:1.056811 1.134475
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('GCT','GCA','GCG','GCC')& 
                          variant_values$nucleo_3_alt %in% c('GAT','GAA','GAG','GAC')) |
                         (variant_values$nucleo_3_alt %in% c('GCT','GCA','GCG','GCC')& 
                            variant_values$nucleo_3_ref %in% c('GAT','GAA','GAG','GAC'))]=TRUE

variant_enrich(variant_values)

##############################Look at GCT, GCA, GCG, GCC -> GTA, GTT, GTC, GTG ################################################## 
#Try percent ASM
#Ref=GC ALT=GT: CI:0.3628457 0.4065706
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('GCT','GCA','GCG','GCC')& 
                          variant_values$nucleo_3_alt %in% c('GTT','GTA','GTG','GTC')) |
                         (variant_values$nucleo_3_alt %in% c('GCT','GCA','GCG','GCC')& 
                            variant_values$nucleo_3_ref %in% c('GTT','GTA','GTG','GTC'))]=TRUE

variant_enrich(variant_values)

##############################Look at GCT, GCA, GCG, GCC -> GGA, GGT, GGC, GGG ################################################## 
#Try percent ASM
#Ref=GC ALT=GG: CI:1.133891 1.221736
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('GCT','GCA','GCG','GCC')& 
                          variant_values$nucleo_3_alt %in% c('GGT','GGA','GGG','GGC')) |
                         (variant_values$nucleo_3_alt %in% c('GCT','GCA','GCG','GCC')& 
                            variant_values$nucleo_3_ref %in% c('GGT','GGA','GGG','GGC'))]=TRUE

variant_enrich(variant_values)

##############################Look at GCT, GCA, GCG, GCC -> GAA, GAT, GAC, GAG ################################################## 
#Try percent ASM
#Ref=GC ALT=GA: CI:1.056811 1.134475
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('GCT','GCA','GCG','GCC')& 
                          variant_values$nucleo_3_alt %in% c('GAT','GAA','GAG','GAC')) |
                         (variant_values$nucleo_3_alt %in% c('GCT','GCA','GCG','GCC')& 
                            variant_values$nucleo_3_ref %in% c('GAT','GAA','GAG','GAC'))]=TRUE

variant_enrich(variant_values)

##############################Look at  TAA, TAT, TAC, TAG >TCT, TCA, TCG, TCC ################################################## 
#Ref=TC ALT=TA: CI:0.9778536 1.0583562, if remove TCG: 0.3430561 0.3966811
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('TCT','TCA','TCG','TCC')& 
                          variant_values$nucleo_3_alt %in% c('TAT','TAA','TAG','TAC')) |
                         (variant_values$nucleo_3_alt %in% c('TCT','TCA','TCG','TCC')& 
                            variant_values$nucleo_3_ref %in% c('TAT','TAA','TAG','TAC'))]=TRUE

variant_enrich(variant_values)

##############################Look at  TAA, TAT, TAC >TCT, TCA, TCC ################################################## 
#Ref=TC ALT=TA if remove TCG: 0.3430561 0.3966811
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('TCT','TCA','TCC')& 
                          variant_values$nucleo_3_alt %in% c('TAT','TAA','TAC')) |
                         (variant_values$nucleo_3_alt %in% c('TCT','TCA','TCC')& 
                            variant_values$nucleo_3_ref %in% c('TAT','TAA','TAC'))]=TRUE

variant_enrich(variant_values)


##############################Look at TAT, TAA, TAG, TAC -> TGA, TGT, TGC, TGG ################################################## 
#Ref=TG ALT=TA: CI:0.369119 0.399706
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('TAT','TAA','TAG','TAC')& 
                          variant_values$nucleo_3_alt %in% c('TGT','TGA','TGG','TGC')) |
                         (variant_values$nucleo_3_alt %in% c('TAT','TAA','TAG','TAC')& 
                            variant_values$nucleo_3_ref %in% c('TGT','TGA','TGG','TGC'))]=TRUE

variant_enrich(variant_values)

##############################Look at TAT, TAA, TAG, TAC -> TTA, TTT, TTC, TTG ################################################## 
#Ref=TT ALT=TA: CI: 0.3284033 0.3817967
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('TAT','TAA','TAG','TAC')& 
                          variant_values$nucleo_3_alt %in% c('TTT','TTA','TTG','TTC')) |
                         (variant_values$nucleo_3_alt %in% c('TAT','TAA','TAG','TAC')& 
                            variant_values$nucleo_3_ref %in% c('TTT','TTA','TTG','TTC'))]=TRUE

variant_enrich(variant_values)




##############################Look at ATT, ATA, ATG, ATC -> AGA, AGT, AGC, AGG for control ################################################## 
#Try percent ASM
#Ref=AT ALT=AG: CI:0.3292249 0.3725984
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('ATT','ATA','ATG','ATC')& 
                          variant_values$nucleo_3_alt %in% c('AGT','AGA','AGG','AGC')) |
                         (variant_values$nucleo_3_alt %in% c('ATT','ATA','ATG','ATC')& 
                            variant_values$nucleo_3_ref %in% c('AGT','AGA','AGG','AGC'))]=TRUE
variant_enrich(variant_values)

##############################Look at ATT, ATA, ATG, ATC -> ACA, ACT, ACC, ACG for control ################################################## 
#Try percent ASM
#Ref=AT ALT=AC: CI:0.3753827 0.4007720
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('ATT','ATA','ATG','ATC')& 
                          variant_values$nucleo_3_alt %in% c('ACT','ACA','ACG','ACC')) |
                         (variant_values$nucleo_3_alt %in% c('ATT','ATA','ATG','ATC')& 
                            variant_values$nucleo_3_ref %in% c('ACT','ACA','ACG','ACC'))]=TRUE
variant_enrich(variant_values)

##############################Look at ATT, ATA, ATG, ATC -> AAA, AAT, AAC, AAG for control ################################################## 
#Try percent ASM
#Ref=AT ALT=AA: CI: 0.3545607 0.4006557
variant_values$variant=FALSE
variant_values$variant[(variant_values$nucleo_3_ref %in% c('ATT','ATA','ATG','ATC')& 
                          variant_values$nucleo_3_alt %in% c('AAT','AAA','AAG','AAC')) |
                         (variant_values$nucleo_3_alt %in% c('ATT','ATA','ATG','ATC')& 
                            variant_values$nucleo_3_ref %in% c('AAT','AAA','AAG','AAC'))]=TRUE
variant_enrich(variant_values)