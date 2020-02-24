# Function to set genome and chromosome lengths to a GR object
setGenomeLengths <- function(GR,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Get genome info
 
  genome.seqinfo <- Seqinfo(genome="hg19")
  genome.seqinfo <- genome.seqinfo[chrsOfInterest]
  GR <- GR[seqnames(GR) %in% chrsOfInterest]
  genome(GR) <- genome(genome.seqinfo)
  seqlevels(GR) <- seqlevels(genome.seqinfo)
  seqlengths(GR) <- seqlengths(genome.seqinfo)
  
  return(GR)
}

# Function that returns subjects and tissues
getSubjectsTissues <- function(inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # H9
  H9_subject <- rep("H9",1)
  H9_subject_labels <- rep("CL1",1)
  H9_tissues <- c(
    "42_embryonic_stem_cell_single"
  )             
  
  # HUES64
  HUES64_subject <- rep("HUES64",4)
  HUES64_subject_labels <- rep("CL2",4)
  HUES64_tissues <- c(
    "stem_27_undifferentiated_paired",
    "ectoderm_paired",
    "endoerm_27_paired",
    "mesoderm_23_paired"
  )             
  
  # skin03
  skin03_subject <- rep("skin03",2)
  skin03_subject_labels <- rep("D11",2)
  skin03_tissues <- c(
    "foreskin_keratinocyte_paired",
    "foreskin_melanocyte_paired"
  )             
  
  # STL001
  stl001_subject <- rep("STL001",11)
  stl001_subject_labels <- rep("D5",11)
  stl001_tissues <- c(
    "Adipose_single",
    "Small_Intestine_single",
    "Bladder_single",
    "Gastric_single",
    "Left_Ventricle_single",
    "Lung_single",
    "Psoas_Muscle_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Spleen_single",
    "Thymus_single"
  )             
  
  # STL002
  stl002_subject <- rep("STL002",11)
  stl002_subject_labels <- rep("D6",11)
  stl002_tissues <- c(
    "Adipose_single",
    "Adrenal_Gland_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Lung_single",
    "Ovary_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Small_Intestine_single",
    "Spleen_single"
  )  
  
  # STL003
  stl003_subject <- rep("STL003",13)
  stl003_subject_labels <- rep("D7",13)
  stl003_tissues <- c(
    "Adipose_Tissue_single",
    "Adrenal_Gland_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Left_Ventricle_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Right_Atrium_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Small_Intestine_single",
    "Spleen_single"
  )             
  
  # STL011
  stl011_subject <- rep("STL011",1)
  stl011_subject_labels <- rep("D8",1)
  stl011_tissues <- c(
    "Liver_single"
  )             
  
  # Create single vectors
  subjects <- c(H9_subject,HUES64_subject,skin03_subject,
                stl001_subject,stl002_subject,stl003_subject,stl011_subject)
  tissues <- c(H9_tissues,HUES64_tissues,skin03_tissues,
               stl001_tissues,stl002_tissues,stl003_tissues,stl011_tissues)
  
  # Return list
  return(list("Subject"=subjects,"Tissues"=tissues))
}

# Function to load all JuliASM output
resultsCpelAsm <- function(inDir,chrsOfInterest=paste("chr",1:22,sep="")){

	# H9
	H9_subject <- rep("H9",1)
	H9_subject_labels <- rep("CL1",1)
	H9_tissues <- c(
	  "42_embryonic_stem_cell_single"
	)             
	H9_tissue_labels <- c(
	  "Embryonic Stem Cell"
	)
	H9_gtex_labels <- c(
	  ""
	)
	H9_diff_labels <- c(
	  "Undifferentiated"
	)

	# HUES64
	HUES64_subject <- rep("HUES64",4)
	HUES64_subject_labels <- rep("CL2",4)
	HUES64_tissues <- c(
	  "stem_27_undifferentiated_paired",
	  "ectoderm_paired",
	  "endoerm_27_paired",
	  "mesoderm_23_paired"
	)             
	HUES64_tissue_labels <- c(
	  "Embyonic Stem Cell",
	  "Ectoderm",
	  "Endoderm",
	  "Mesoderm"
	)
	HUES64_gtex_labels <- c("","","","")
	HUES64_diff_labels <- c(
	  "Undifferentiated",
	  "Semidifferentiated",
	  "Semidifferentiated",
	  "Semidifferentiated"
	)

	# skin03
	skin03_subject <- rep("skin03",2)
	skin03_subject_labels <- rep("D11",2)
	skin03_tissues <- c(
	  "foreskin_keratinocyte_paired",
	  "foreskin_melanocyte_paired"
	)             
	skin03_tissue_labels <- c(
	  "Foreskin Keratinocyte",
	  "Foreskin Melanocyte"
	)
	skin03_gtex_labels <- c("","")
	skin03_diff_labels <- rep("Differentiated",length(skin03_subject))
	
	# # HuFGM02
	# HuFGM02_subject <- rep("HuFGM02",2)
	# HuFGM02_subject_labels <- rep("D4",2)
	# HuFGM02_tissues <- c(
	#   "brain_cerebellum_tissue_paired",
	#   "brain_germinal_matrix_tissue_paired"
	# )             
	# HuFGM02_tissue_labels <- c(
	#   "Brain Cerebellum",
	#   "Brain Germinal Matrix"
	# )
	# HuFGM02_gtex_labels <- c(
	#   "Brain",
	#   "Brain"
	# )
	# HuFGM02_diff_labels <- rep("Differentiated",length(HuFGM02_subject))

	# STL001
	stl001_subject <- rep("STL001",11)
	stl001_subject_labels <- rep("D5",11)
	stl001_tissues <- c(
	  "Adipose_single",
	  "Small_Intestine_single",
	  "Bladder_single",
	  "Gastric_single",
	  "Left_Ventricle_single",
	  "Lung_single",
	  "Psoas_Muscle_single",
	  "Right_Ventricle_single",
	  "Sigmoid_Colon_single",
	  "Spleen_single",
	  "Thymus_single"
	)             
	stl001_tissue_labels <- c(
	  "Adipose",
	  "Small Intestine",
	  "Bladder",
	  "Gastric",
	  "Left Ventricle",
	  "Lung",
	  "Psoas Muscle",
	  "Right Ventricle",
	  "Sigmoid Colon",
	  "Spleen",
	  "Thymus"
	)
	stl001_gtex_labels <- c(
	  "Adipose_Subcutaneous",
	  "",
	  "",
	  "",
	  "Heart_Left_Ventricle",
	  "Lung",
	  "",
	  "",
	  "Colon_Transverse",
	  "",
	  ""
	)
	stl001_diff_labels <- rep("Differentiated",length(stl001_subject))

	# STL002
	stl002_subject <- rep("STL002",11)
	stl002_subject_labels <- rep("D6",11)
	stl002_tissues <- c(
	  "Adipose_single",
	  "Adrenal_Gland_single",
	  "Aorta_single",
	  "Esophagus_single",
	  "Gastric_single",
	  "Lung_single",
	  "Ovary_single",
	  "Pancreas_single",
	  "Psoas_Muscle_single",
	  "Small_Intestine_single",
	  "Spleen_single"
	)  
	stl002_tissue_labels <- c(
	  "Adipose",
	  "Adrenal Gland",
	  "Aorta",
	  "Esophagus",
	  "Gastric",
	  "Lung",
	  "Ovary",
	  "Pancreas",
	  "Psoas Muscle",
	  "Small Intestine",
	  "Spleen"
	)  
	stl002_gtex_labels <- c(
	  "Adipose_Subcutaneous",
	  "Adrenal_Gland",
	  "Artery_Aorta",
	  "Esophagus_Mucosa",
	  "",
	  "Lung",
	  "Uterus",
	  "Pancreas",
	  "",
	  "",
	  ""
	)
	stl002_diff_labels <- rep("Differentiated",length(stl002_subject))

	# STL003
	stl003_subject <- rep("STL003",13)
	stl003_subject_labels <- rep("D7",13)
	stl003_tissues <- c(
	  "Adipose_Tissue_single",
	  "Adrenal_Gland_single",
	  "Aorta_single",
	  "Esophagus_single",
	  "Gastric_single",
	  "Left_Ventricle_single",
	  "Pancreas_single",
	  "Psoas_Muscle_single",
	  "Right_Atrium_single",
	  "Right_Ventricle_single",
	  "Sigmoid_Colon_single",
	  "Small_Intestine_single",
	  "Spleen_single"
	)             
	stl003_tissue_labels <- c(                 
	  "Adipose",
	  "Adrenal Gland",
	  "Aorta",
	  "Esophagus",
	  "Gastric",
	  "Left Ventricle",
	  "Pancreas",
	  "Psoas Muscle",
	  "Right Atrium",
	  "Right Ventricle",
	  "Sigmoid Colon",
	  "Small Intestine",
	  "Spleen"
	)
	stl003_gtex_labels <- c(                 
	  "Adipose_Subcutaneous",
	  "Adrenal_Gland",
	  "Artery_Aorta",
	  "Esophagus",
	  "",
	  "Heart_Left_Ventricle",
	  "Pancreas",
	  "",
	  "",
	  "",
	  "Colon_Transverse",
	  "",
	  ""
	)
	stl003_diff_labels <- rep("Differentiated",length(stl003_subject))
	
	# STL011
	stl011_subject <- rep("STL011",1)
	stl011_subject_labels <- rep("D8",1)
	stl011_tissues <- c(
	  "Liver_single"
	)             
	stl011_tissue_labels <- c(
	  "Liver"
	)
	stl011_gtex_labels <- c(
	  "Liver"
	)
	stl011_diff_labels <- rep("Differentiated",length(stl011_subject))
	
	# Create single vectors
	subjects <- c(H9_subject,HUES64_subject,skin03_subject,HuFGM02_subject,
	              stl001_subject,stl002_subject,stl003_subject,stl011_subject)
	tissues <- c(H9_tissues,HUES64_tissues,skin03_tissues,HuFGM02_tissues,
	             stl001_tissues,stl002_tissues,stl003_tissues,stl011_tissues)
	subject_labels <- c(H9_subject_labels,HUES64_subject_labels,skin03_subject_labels,
	                    HuFGM02_subject_labels,stl001_subject_labels,stl002_subject_labels,
	                    stl003_subject_labels,stl011_subject_labels) 
	tissue_labels <- c(H9_tissue_labels,HUES64_tissue_labels,skin03_tissue_labels,
	                   HuFGM02_tissue_labels,stl001_tissue_labels,stl002_tissue_labels,
	                   stl003_tissue_labels,stl011_tissue_labels)
	gtex_labels <- c(H9_gtex_labels,HUES64_gtex_labels,skin03_gtex_labels,
			             HuFGM02_gtex_labels,stl001_gtex_labels,stl002_gtex_labels,
			             stl003_gtex_labels,stl011_gtex_labels)
	diff_labels <- c(H9_diff_labels,HUES64_diff_labels,skin03_diff_labels,
	                 HuFGM02_diff_labels,stl001_diff_labels,stl002_diff_labels,
	                 stl003_diff_labels,stl011_diff_labels)

	# Read in all data
	GRs <- GRanges()
	for (i in 1:length(subjects)) {
	  # Print sample being loaded
	  print(paste("Loading sample:",subjects[i],tissues[i]))
	  # dmml
	  dmml <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_dmml_pvals.bedGraph",sep=""))
	  dmml$Subject <- subjects[i]
	  dmml$SubjectLabel <- subject_labels[i]
	  dmml$Tissue <- tissue_labels[i]
	  dmml$GTEx <- gtex_labels[i]
	  dmml$State <- diff_labels[i]
	  dmml$Statistic <- "dMML"
	  dmml$Value <- dmml$score
	  dmml$pvalue <- as.numeric(dmml$NA.)
	  dmml$score <- NULL
	  dmml$NA. <- NULL
	  dmml <- dmml[!duplicated(dmml[,c()])]
	  dmml <- resize(dmml, width(dmml) + 1, fix="end")
	  dmml <- resize(dmml, width(dmml) + 1, fix="start")
	  GRs <- append(GRs,dmml)
	  # dnme
	  dnme <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_dnme_pvals.bedGraph",sep=""))
	  dnme$Subject <- subjects[i]
	  dnme$SubjectLabel <- subject_labels[i]
	  dnme$Tissue <- tissue_labels[i]
	  dnme$GTEx <- gtex_labels[i]
	  dnme$State <- diff_labels[i]
	  dnme$Statistic <- "dNME"
	  dnme$Value <- dnme$score
	  dnme$pvalue <- as.numeric(dnme$NA.)
	  dnme$score <- NULL
	  dnme$NA. <- NULL
	  dnme <- dnme[!duplicated(dnme[,c()])]
	  dnme <- resize(dnme,width(dnme) + 1, fix="end")
	  dnme <- resize(dnme,width(dnme) + 1, fix="start")
	  GRs <- append(GRs,dnme)
	  # uc
	  uc <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_uc_pvals.bedGraph",sep=""))
	  uc$Subject <- subjects[i]
	  uc$SubjectLabel <- subject_labels[i]
	  uc$Tissue <- tissue_labels[i]
	  uc$GTEx <- gtex_labels[i]
	  uc$State <- diff_labels[i]
	  uc$Statistic <- "UC"
	  uc$Value <- uc$score
	  uc$pvalue <- as.numeric(uc$NA.)
	  uc$score <- NULL
	  uc$NA. <- NULL
	  uc <- uc[!duplicated(uc[,c()])]
	  uc <- resize(uc,width(uc) + 1, fix="end")
	  uc <- resize(uc,width(uc) + 1, fix="start")
	  GRs <- append(GRs,uc)
	}

	# Add significance to GRs
	GRs <- GRs[!is.na(GRs$pvalue)]
	GRs$ASM <- NA
	GRs[GRs$pvalue<=0.05]$ASM <- "Yes"
	GRs[GRs$pvalue>0.05]$ASM <- "No"

	# Add sample field
	GRs$Sample <- paste(GRs$Tissue,"-",GRs$Subject)
	
	# Add genome info 
	GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
	
	# Return
	return(GRs)
	
}

# Function to load Nat Comms CpelAsm output
resultsCpelAsmNatComms <- function(inDir,chrsOfInterest=paste("chr",1:22,sep="")){

  # STL003
  stl003_subject <- rep("STL003",10)
  stl003_subject_labels <- rep("D7",10)
  stl003_tissues <- c(
    "Adipose_Tissue_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Small_Intestine_single",
    "Spleen_single"
  )             
  stl003_tissue_labels <- c(                 
    "Adipose",
    "Aorta",
    "Esophagus",
    "Gastric",
    "Pancreas",
    "Psoas Muscle",
    "Right Ventricle",
    "Sigmoid Colon",
    "Small Intestine",
    "Spleen"
  )
  stl003_gtex_labels <- c(                 
    "Adipose_Subcutaneous",
    "Artery_Aorta",
    "Esophagus",
    "",
    "Pancreas",
    "",
    "",
    "Colon_Transverse",
    "",
    ""
  )
  stl003_diff_labels <- rep("Differentiated",length(stl003_subject))

  # Create single vectors
  subjects <- stl003_subject
  tissues <- stl003_tissues
  subject_labels <- stl003_subject_labels
  tissue_labels <- stl003_tissue_labels
  gtex_labels <- stl003_gtex_labels
  diff_labels <- stl003_diff_labels
  
  # Read in all data
  GRs <- GRanges()
  for (i in 1:length(subjects)) {
    # Print sample being loaded
    print(paste("Loading sample:",subjects[i],tissues[i]))
    # dmml
    dmml <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_dmml_pvals.bedGraph",sep=""))
    dmml$Subject <- subjects[i]
    dmml$SubjectLabel <- subject_labels[i]
    dmml$Tissue <- tissue_labels[i]
    dmml$GTEx <- gtex_labels[i]
    dmml$State <- diff_labels[i]
    dmml$Statistic <- "dMML"
    dmml$Value <- dmml$score
    dmml$pvalue <- as.numeric(dmml$NA.)
    dmml$score <- NULL
    dmml$NA. <- NULL
    dmml <- dmml[!duplicated(dmml[,c()])]
    dmml <- resize(dmml, width(dmml) + 1, fix="end")
    dmml <- resize(dmml, width(dmml) + 1, fix="start")
    GRs <- append(GRs,dmml)
    # dnme
    dnme <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_dnme_pvals.bedGraph",sep=""))
    dnme$Subject <- subjects[i]
    dnme$SubjectLabel <- subject_labels[i]
    dnme$Tissue <- tissue_labels[i]
    dnme$GTEx <- gtex_labels[i]
    dnme$State <- diff_labels[i]
    dnme$Statistic <- "dNME"
    dnme$Value <- dnme$score
    dnme$pvalue <- as.numeric(dnme$NA.)
    dnme$score <- NULL
    dnme$NA. <- NULL
    dnme <- dnme[!duplicated(dnme[,c()])]
    dnme <- resize(dnme,width(dnme) + 1, fix="end")
    dnme <- resize(dnme,width(dnme) + 1, fix="start")
    GRs <- append(GRs,dnme)
    # uc
    uc <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_uc_pvals.bedGraph",sep=""))
    uc$Subject <- subjects[i]
    uc$SubjectLabel <- subject_labels[i]
    uc$Tissue <- tissue_labels[i]
    uc$GTEx <- gtex_labels[i]
    uc$State <- diff_labels[i]
    uc$Statistic <- "UC"
    uc$Value <- uc$score
    uc$pvalue <- as.numeric(uc$NA.)
    uc$score <- NULL
    uc$NA. <- NULL
    uc <- uc[!duplicated(uc[,c()])]
    uc <- resize(uc,width(uc) + 1, fix="end")
    uc <- resize(uc,width(uc) + 1, fix="start")
    GRs <- append(GRs,uc)
  }
  
  # Add significance to GRs
  GRs <- GRs[!is.na(GRs$pvalue)]
  GRs$ASM <- NA
  GRs[GRs$pvalue<=0.05]$ASM <- "Yes"
  GRs[GRs$pvalue>0.05]$ASM <- "No"
  
  # Add sample field
  GRs$Sample <- paste(GRs$Tissue,"-",GRs$Subject)
  
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  
  # Return
  return(GRs)
  
}

# Function to load JuliASM output
resultsCpelAllele <- function(inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # H9
  H9_subject <- rep("H9",1)
  H9_subject_labels <- rep("CL1",1)
  H9_tissues <- c(
    "42_embryonic_stem_cell_single"
  )             
  H9_tissue_labels <- c(
    "Embryonic Stem Cell"
  )
  H9_gtex_labels <- c(
    ""
  )
  H9_diff_labels <- c(
    "Undifferentiated"
  )
  
  # HUES64
  HUES64_subject <- rep("HUES64",4)
  HUES64_subject_labels <- rep("CL2",4)
  HUES64_tissues <- c(
    "stem_27_undifferentiated_paired",
    "ectoderm_paired",
    "endoerm_27_paired",
    "mesoderm_23_paired"
  )
  HUES64_tissue_labels <- c(
    "Embyonic Stem Cell",
    "Ectoderm",
    "Endoderm",
    "Mesoderm"
  )
  HUES64_gtex_labels <- c("","","","")
  HUES64_diff_labels <- c(
    "Undifferentiated",
    "Semidifferentiated",
    "Semidifferentiated",
    "Semidifferentiated"
  )

  # skin03
  skin03_subject <- rep("skin03",2)
  skin03_subject_labels <- rep("D11",2)
  skin03_tissues <- c(
    "foreskin_keratinocyte_paired",
    "foreskin_melanocyte_paired"
  )
  skin03_tissue_labels <- c(
    "Foreskin Keratinocyte",
    "Foreskin Melanocyte"
  )
  skin03_gtex_labels <- c("","")
  skin03_diff_labels <- rep("Differentiated",length(skin03_subject))

  # STL001
  stl001_subject <- rep("STL001",11)
  stl001_subject_labels <- rep("D5",11)
  stl001_tissues <- c(
    "Adipose_single",
    "Small_Intestine_single",
    "Bladder_single",
    "Gastric_single",
    "Left_Ventricle_single",
    "Lung_single",
    "Psoas_Muscle_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Spleen_single",
    "Thymus_single"
  )
  stl001_tissue_labels <- c(
    "Adipose",
    "Small Intestine",
    "Bladder",
    "Gastric",
    "Left Ventricle",
    "Lung",
    "Psoas Muscle",
    "Right Ventricle",
    "Sigmoid Colon",
    "Spleen",
    "Thymus"
  )
  stl001_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "",
    "",
    "",
    "Heart_Left_Ventricle",
    "Lung",
    "",
    "",
    "Colon_Transverse",
    "",
    ""
  )
  stl001_diff_labels <- rep("Differentiated",length(stl001_subject))

  # STL002
  stl002_subject <- rep("STL002",11)
  stl002_subject_labels <- rep("D6",11)
  stl002_tissues <- c(
    "Adipose_single",
    "Adrenal_Gland_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Lung_single",
    "Ovary_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Small_Intestine_single",
    "Spleen_single"
  )
  stl002_tissue_labels <- c(
    "Adipose",
    "Adrenal Gland",
    "Aorta",
    "Esophagus",
    "Gastric",
    "Lung",
    "Ovary",
    "Pancreas",
    "Psoas Muscle",
    "Small Intestine",
    "Spleen"
  )
  stl002_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Esophagus_Mucosa",
    "",
    "Lung",
    "Uterus",
    "Pancreas",
    "",
    "",
    ""
  )
  stl002_diff_labels <- rep("Differentiated",length(stl002_subject))

  # STL003
  stl003_subject <- rep("STL003",13)
  stl003_subject_labels <- rep("D7",13)
  stl003_tissues <- c(
    "Adipose_Tissue_single",
    "Adrenal_Gland_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Left_Ventricle_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Right_Atrium_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Small_Intestine_single",
    "Spleen_single"
  )
  stl003_tissue_labels <- c(
    "Adipose",
    "Adrenal Gland",
    "Aorta",
    "Esophagus",
    "Gastric",
    "Left Ventricle",
    "Pancreas",
    "Psoas Muscle",
    "Right Atrium",
    "Right Ventricle",
    "Sigmoid Colon",
    "Small Intestine",
    "Spleen"
  )
  stl003_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Esophagus",
    "",
    "Heart_Left_Ventricle",
    "Pancreas",
    "",
    "",
    "",
    "Colon_Transverse",
    "",
    ""
  )
  stl003_diff_labels <- rep("Differentiated",length(stl003_subject))

  # STL011
  stl011_subject <- rep("STL011",1)
  stl011_subject_labels <- rep("D8",1)
  stl011_tissues <- c(
    "Liver_single"
  )
  stl011_tissue_labels <- c(
    "Liver"
  )
  stl011_gtex_labels <- c(
    "Liver"
  )
  stl011_diff_labels <- rep("Differentiated",length(stl011_subject))
  
  # Create single vectors
  subjects <- c(H9_subject,HUES64_subject,skin03_subject,stl001_subject,
                stl002_subject,stl003_subject,stl011_subject)
  tissues <- c(H9_tissues,HUES64_tissues,skin03_tissues,stl001_tissues,
               stl002_tissues,stl003_tissues,stl011_tissues)
  subject_labels <- c(H9_subject_labels,HUES64_subject_labels,skin03_subject_labels,
                      stl001_subject_labels,stl002_subject_labels,stl003_subject_labels,
                      stl011_subject_labels)
  tissue_labels <- c(H9_tissue_labels,HUES64_tissue_labels,skin03_tissue_labels,
                     stl001_tissue_labels,stl002_tissue_labels,stl003_tissue_labels,
                     stl011_tissue_labels)
  gtex_labels <- c(H9_gtex_labels,HUES64_gtex_labels,skin03_gtex_labels,
                   stl001_gtex_labels,stl002_gtex_labels,stl003_gtex_labels,
                   stl011_gtex_labels)
  diff_labels <- c(H9_diff_labels,HUES64_diff_labels,skin03_diff_labels,stl001_diff_labels,
                   stl002_diff_labels,stl003_diff_labels,stl011_diff_labels)
  
  # Read in all data
  GRs <- GRanges()
  for (i in 1:length(subjects)) {
    # Print sample being loaded
    print(paste("Loading sample:",subjects[i],tissues[i]))
    # mml1
    mml1 <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_mml1.bedGraph",sep=""))
    mml1$Subject <- subjects[i]
    mml1$SubjectLabel <- subject_labels[i]
    mml1$Tissue <- tissue_labels[i]
    mml1$GTEx <- gtex_labels[i]
    mml1$State <- diff_labels[i]
    mml1$Statistic <- "MML"
    mml1$Genome <- "1"
    mml1$Value <- mml1$score
    mml1$score <- NULL
    mml1$NA. <- NULL
    mml1$NA.1 <- NULL
    mml1 <- mml1[!duplicated(mml1[,c()])]
    mml1 <- resize(mml1, width(mml1) + 1, fix="end")
    mml1 <- resize(mml1, width(mml1) + 1, fix="start")
    GRs <- append(GRs,mml1)
    # mml2
    mml2 <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_mml2.bedGraph",sep=""))
    mml2$Subject <- subjects[i]
    mml2$SubjectLabel <- subject_labels[i]
    mml2$Tissue <- tissue_labels[i]
    mml2$GTEx <- gtex_labels[i]
    mml2$State <- diff_labels[i]
    mml2$Statistic <- "MML"
    mml2$Genome <- "2"
    mml2$Value <- mml2$score
    mml2$score <- NULL
    mml2$NA. <- NULL
    mml2$NA.1 <- NULL
    mml2 <- mml2[!duplicated(mml2[,c()])]
    mml2 <- resize(mml2, width(mml2) + 1, fix="end")
    mml2 <- resize(mml2, width(mml2) + 1, fix="start")
    GRs <- append(GRs,mml2)
    # nme1
    nme1 <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_nme1.bedGraph",sep=""))
    nme1$Subject <- subjects[i]
    nme1$SubjectLabel <- subject_labels[i]
    nme1$Tissue <- tissue_labels[i]
    nme1$GTEx <- gtex_labels[i]
    nme1$State <- diff_labels[i]
    nme1$Statistic <- "NME"
    nme1$Genome <- "1"
    nme1$Value <- nme1$score
    nme1$score <- NULL
    nme1$NA. <- NULL
    nme1$NA.1 <- NULL
    nme1 <- nme1[!duplicated(nme1[,c()])]
    nme1 <- resize(nme1,width(nme1) + 1, fix="end")
    nme1 <- resize(nme1,width(nme1) + 1, fix="start")
    GRs <- append(GRs,nme1)
    # nme2
    nme2 <- import.bedGraph(paste(inDir,subjects[i],"_",tissues[i],"_phased_nme2.bedGraph",sep=""))
    nme2$Subject <- subjects[i]
    nme2$SubjectLabel <- subject_labels[i]
    nme2$Tissue <- tissue_labels[i]
    nme2$GTEx <- gtex_labels[i]
    nme2$State <- diff_labels[i]
    nme2$Statistic <- "NME"
    nme2$Genome <- "2"
    nme2$Value <- nme2$score
    nme2$score <- NULL
    nme2$NA. <- NULL
    nme2$NA.1 <- NULL
    nme2 <- nme2[!duplicated(nme2[,c()])]
    nme2 <- resize(nme2,width(nme2) + 1, fix="end")
    nme2 <- resize(nme2,width(nme2) + 1, fix="start")
    GRs <- append(GRs,nme2)
  }
  
  # Add sample field
  GRs$Sample <- paste(GRs$Tissue,"-",GRs$Subject)
  
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  
  # Return
  return(GRs)
  
}
getGeneralFeats <- function(CpGdir,enhancerDir,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Define list of feature GRs
  outGR <- GRangesList()
  GRtemp <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  outGR[["genome-wide"]] <- setGenomeLengths(GRtemp)
  
  # below, need to enter the path to CpG islands BED file (BED file has header removed)
  cpg_islands <- readRDS(paste(CpGdir,"cpg_islands_hg19.rds",sep=""))
  outGR[["CpG island"]] <- setGenomeLengths(cpg_islands)
  
  # extract the shore defined by 2000 bp upstream and downstream of cpg islands
  shore1 <- flank(cpg_islands, 2000)
  shore2 <- flank(cpg_islands,2000,FALSE)
  shore1_2 <- reduce(c(shore1,shore2))
  
  # extract the features (ranges) that are present in shores only and not in
  # cpg_islands (ie., shores not overlapping islands)
  cpgi_shores <- setdiff(shore1_2, cpg_islands)
  outGR[["CpG shore"]] <- setGenomeLengths(cpgi_shores)
  
  # extract the shore defined by 4000 bp upstream and downstream of cpg islands
  shelves1 <- flank(cpg_islands, 4000)
  shelves2 <- flank(cpg_islands,4000,FALSE)
  shelves1_2 <- reduce(c(shelves1,shelves2))
  
  # create a set of ranges consisting CpG Islands, Shores
  island_shores <- c(cpg_islands,cpgi_shores)
  
  # extract the features (ranges) that are present in shelves only
  # and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
  cpgi_shelves <- setdiff(shelves1_2, island_shores)
  outGR[["CpG shelf"]] <- setGenomeLengths(cpgi_shelves)
  
  # Open sea
  open_sea <- setdiff(outGR[["genome-wide"]],c(outGR[["CpG island"]],outGR[["CpG shore"]],outGR[["CpG shelf"]]))
  outGR[["CpG open sea"]] <- setGenomeLengths(open_sea)
  
  # CG rich/poor
  cgRich <- reduce(c(outGR[["CpG island"]],outGR[["CpG shore"]])) 
  outGR[["CG rich"]] <- setGenomeLengths(cgRich)
  cgPoor <- setdiff(outGR[["genome-wide"]],c(outGR[["CpG island"]],outGR[["CpG shore"]]))
  outGR[["CG poor"]] <- setGenomeLengths(cgPoor)
  
  # Enhancers 
  #enhancers <- import.bed(paste(enhancerDir,"enhancers.bed",sep=""))[,c()]
  #outGR[["enhancer"]] <- setGenomeLengths(enhancers)
  
  # Other generic features
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- genes(txdb)
  outGR[["gene body"]] <- setGenomeLengths(genes)
  exons <- exons(txdb)
  outGR[["exon"]] <- setGenomeLengths(exons[,c()])
  introns <- intronicParts(txdb)
  outGR[["intron"]] <- setGenomeLengths(introns[,c()])
  intergenic <- setdiff(outGR[["genome-wide"]],outGR[["gene body"]],ignore.strand=TRUE)
  outGR[["intergenic"]] <- setGenomeLengths(intergenic)
  proms <- promoters(genes,upstream=2000,downstream=2000)
  outGR[["promoter"]] <- setGenomeLengths(proms)
  
  # (Non)-transcriptional
  #trans <- reduce(c(outGR[["exon"]],outGR[["promoter"]],outGR[["enhancer"]]))
  #outGR[["transcriptional"]] <- setGenomeLengths(trans)
  #nontrans <- setdiff(outGR[["genome-wide"]],trans,ignore.strand=TRUE)
  #outGR[["non-transcriptional"]] <- setGenomeLengths(nontrans)
  
  # Gene name mapping
  geneBodyNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["gene body"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["gene body"]]$gene_name <- geneBodyNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["promoter"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["promoter"]]$gene_name <- promNameMap$SYMBOL
  
  # Return
  return(outGR)
  
}
# Function to read in a GFF file
readSubGff <- function(inDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Import GFF
  outGR <- import.gff(paste(inDir,sub,"_het.juliasm.gff",sep=""))
  
  # Retain a the required subset
  outGR <- outGR[,c("N","CpGs","hetCpGg1","hetCpGg2")]
  
  # Add subject metadata column
  outGR$Subject <- sub
  
  # Make N numeric and filter
  outGR$N <- as.numeric(outGR$N)
  outGR <- outGR[outGR$N<=20]
  
  # Add genome info 
  outGR <- setGenomeLengths(outGR,chrsOfInterest=chrsOfInterest)
  
  # Return GR with all haplotypes
  return(outGR)
  
}

# Function to read in all GFF files
readAllGff <- function(inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Loop over all GFF files
  outGR <- GRanges()
  subjects <- c("H9","HUES64","skin03","HuFGM02","STL001","STL002","STL003","STL011")
  for (sub in subjects) {
    # Import GFF file
    cat('importing',sub,'\n')
    tmpGR <- readSubGff(inDir,sub,chrsOfInterest)
    # Append to output GR
    outGR <- append(outGR,tmpGR)
  }
  
  # Add data and ASM fields
  outGR$ASM <- 0
  outGR$Data <- 0
  
  # Return GR with all haplotypes
  return(outGR)
  
}

# Function for coverting rownames to granges
rownames2Granges<-function(rn){
  stp=strsplit(rn,':')[[1]]
  chr=stp[1]
  start_end=strsplit(stp[2],'-')[[1]]
  start=as.numeric(start_end[1])
  end=as.numeric(start_end[2])
  return(GRanges(seqnames=chr,ranges=IRanges(start,end)))
}

# Function to read in enhancer data in Rdata format (very slow)
readEnhancer <- function(enhancerDir){
  
  #Function to convert colnames to granges
  load(paste(enhancerDir,"enhancers_intersect.RData",sep=""))
  enhancer_gr_all <- do.call('c',lapply(rownames(max_states),rownames2Granges))
  
  # Return
  return(enhancer_gr_all)
}

# Function to get general genomic features GR list
getGeneralFeats_CpG <- function(CpGdir,enhancerDir='',chrsOfInterest=paste("chr",1:22,sep="")){

  # Features included
  featureNickNames <- c("genome-wide","CpG island","CpG shore","CpG shelf","CpG open sea",
                        "gene body","exon","intron","intergenic")
  
  # Define list of feature GRs
  outGR <- GRangesList()
  GRtemp <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  
  outGR[["genome-wide"]] <- setGenomeLengths(GRtemp)
  CpG_all <- readRDS(paste(CpGdir,"CpG_hg19.rds",sep=""))
  CpG_all<-setGenomeLengths(CpG_all)
  # below, need to enter the path to CpG islands BED file (BED file has header removed)
  cpg_islands <- readRDS(paste(CpGdir,"cpg_islands.rds",sep=""))
  cpg_islands<-subsetByOverlaps(CpG_all,cpg_islands)
  outGR[["CpG island"]] <- setGenomeLengths(cpg_islands)
  
  # extract the shore defined by 2000 bp upstream and downstream of cpg islands
  shore1 <- flank(cpg_islands, 2000)
  shore2 <- flank(cpg_islands,2000,FALSE)
  shore1_2 <- reduce(c(shore1,shore2))
  
  # extract the features (ranges) that are present in shores only and not in
  # cpg_islands (ie., shores not overlapping islands)
  cpgi_shores <- setdiff(shore1_2, cpg_islands)
  olap=findOverlaps(CpG_all,cpgi_shores)
  cpgi_shores<-subsetByOverlaps(CpG_all,cpgi_shores)
  outGR[["CpG shore"]] <- setGenomeLengths(cpgi_shores)
  
  # extract the shore defined by 4000 bp upstream and downstream of cpg islands
  shelves1 <- flank(cpg_islands, 4000)
  shelves2 <- flank(cpg_islands,4000,FALSE)
  shelves1_2 <- reduce(c(shelves1,shelves2))
  
  # create a set of ranges consisting CpG Islands, Shores
  island_shores <- c(cpg_islands,cpgi_shores)
  
  # extract the features (ranges) that are present in shelves only
  # and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
  cpgi_shelves <- setdiff(shelves1_2, island_shores)
  cpgi_shelves<-subsetByOverlaps(CpG_all,cpgi_shelves)
  outGR[["CpG shelf"]] <- setGenomeLengths(cpgi_shelves)
  
  # Open sea
  open_sea <- setdiff(outGR[["genome-wide"]],c(outGR[["CpG island"]],outGR[["CpG shore"]],outGR[["CpG shelf"]]))
  open_sea<-subsetByOverlaps(CpG_all,open_sea)
  outGR[["CpG open sea"]] <- setGenomeLengths(open_sea)
  
  # Enhancers 
  #enhancers <- import.bed(paste(enhancerDir,"enhancers.bed",sep=""))[,c()]
  
  #outGR[["enhancer"]] <- setGenomeLengths(enhancers)
  
  # Other generic features
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- genes(txdb)
  outGR[["gene body"]] <- setGenomeLengths(genes)
  exons <- exons(txdb)
  outGR[["exon"]] <- setGenomeLengths(exons[,c()])
  introns <- intronicParts(txdb)
  outGR[["intron"]] <- setGenomeLengths(introns[,c()])
  intergenic <- setdiff(outGR[["genome-wide"]],outGR[["gene body"]],ignore.strand=TRUE)
  outGR[["intergenic"]] <- setGenomeLengths(intergenic)
  proms <- promoters(genes,upstream=2000,downstream=2000)
  outGR[["promoter"]] <- setGenomeLengths(proms)
  
  # Gene name mapping
  geneBodyNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["gene body"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["gene body"]]$gene_name <- geneBodyNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["promoter"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["promoter"]]$gene_name <- promNameMap$SYMBOL
  
  # Return
  return(outGR)
  
}

# Function to cross CpelASM output with genomic features
crossCpelWithFeatures <- function(cpelGrs,featGrs){
  
  # Get overlaps
  cpelWithFeatsGrs <- GRanges()
  for(feature in names(featGrs)){
    print(paste("Working on ",feature,sep=""))
    # Get feature GR
    GRs_tmp <- subsetByOverlaps(cpelGRs,featGrs[[feature]],type="within")
    GRs_tmp$Feature <- feature
    cpelWithFeatsGrs <- append(cpelWithFeatsGrs,GRs_tmp)
  }
  
  # Retain the relevant part of the data frame
  cpelWithFeatsDf <- as.data.frame(cpelWithFeatsGrs)
  cpelWithFeatsDf <- cpelWithFeatsDf[,6:ncol(cpelWithFeatsDf)]
  cpelWithFeatsDf$Feature <- factor(cpelWithFeatsDf$Feature,levels=c("genome-wide","CpG island","CpG shore",
                                    "CpG shelf","CpG open sea","gene body","exon","intron","intergenic",
                                    "promoter","enhancer"))
  
  # Return crossed data frame
  return(cpelWithFeatsDf)
  
}

# Function to get position of all CpG sites in hg19
getCpgSitesH19 <- function(chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cpgr <- do.call(c,lapply(1:length(chrsOfInterest), function(x) GRanges(names(Hsapiens)[x],IRanges(cgs[[x]],width=2))))
  
  # Set genome and seqlengths
  cpgr <- setGenomeLengths(cpgr)
  
  # Return CpG site GR
  return(cpgr)
  
}

# Function to compute CpG density in each haplotype of GR passed
cpgDensity <- function(GR,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Normalize by the number of nucleotides
  GR$CpgDens <- GR$N/(end(GR)-start(GR)+1)
    
  # Return GR
  return(GR)

}

# Function to get break-down on ASM events for each N
gffFormat <- function(inDir,cpelGRs,testStat,subjects=c("H9","HUES64","skin03","HuFGM02","STL001","STL002","STL003","STL011")){

  # Read in gffGR
  gffGR <- readAllGff(inDir)
  subjects <- c("H9","HUES64","skin03","HuFGM02","STL001","STL002","STL003","STL011")
  
  # Loop over all subjects
  outGR <- GRanges()
  for (subject in subjects) {
    # Get haplotypes pertaining to subject from GFF file
    tmpGR <- gffGR[gffGR$Subject==subject]
    # Find for each haplotype the number of times it had data
    olaps <- findOverlaps(cpelGRs[(cpelGRs$Subject==subject)&(cpelGRs$Statistic==testStat)],tmpGR,type="within")
    dataTable <- as.data.frame(table(subjectHits(olaps)))
    tmpGR[as.numeric(dataTable[,1])]$Data <- as.numeric(dataTable[,2])
    # Count many of the ones with data resulted in ASM
    olaps <- findOverlaps(cpelGRs[(cpelGRs$ASM=="Yes")&(cpelGRs$Subject==subject)&(cpelGRs$Statistic==testStat)],tmpGR,type="within")
    dataTable <- as.data.frame(table(subjectHits(olaps)))
    tmpGR[as.numeric(dataTable[,1])]$ASM <- as.numeric(dataTable[,2])
    # Append to output
    outGR <- append(outGR,tmpGR)
  }
  
  # Convert to binary
  outGR$ASM_bin <- NA
  outGR[outGR$ASM>0]$ASM_bin <- "At least one ASM"
  outGR[outGR$ASM==0]$ASM_bin <- "No ASM"
  outGR$ASM_bin <- factor(outGR$ASM_bin,levels=c("No ASM","At least one ASM"))

  # Return
  return(outGR)
  
}

# Function to test for enrichment for a generic feature in GR
testEnrichmentFeature <- function(dataGR,featureGR){
  
  # Find ranges overlapping with feature
  olaps <- findOverlaps(dataGR,featureGR,type="any",select="all")
  indFeature <- queryHits(olaps)
  featureData <- dataGR[indFeature]
  complementaryData <- dataGR[-indFeature]
    
  # Enrichment of dMML-HASM in feature
  featureDmml <- featureData[featureData$Statistic=="dMML"]
  complementaryDmml <- complementaryData[complementaryData$Statistic=="dMML"]
  contTableDmml <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableDmml) <- c("Feature","Complementary")
  contTableDmml[1,]$ASM <- sum(featureDmml$ASM=="Yes")
  contTableDmml[1,]$nonASM <- sum(featureDmml$ASM=="No")
  contTableDmml[2,]$ASM <- sum(complementaryDmml$ASM=="Yes")
  contTableDmml[2,]$nonASM <- sum(complementaryDmml$ASM=="No")
  dmmlFisher <- fisher.test(contTableDmml)
  
  # Enrichment of dNME-HASM in feature
  featureDnme <- featureData[featureData$Statistic=="dNME"]
  complementaryDnme <- complementaryData[complementaryData$Statistic=="dNME"]
  contTableDnme <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableDnme) <- c("Feature","Complementary")
  contTableDnme[1,]$ASM <- sum(featureDnme$ASM=="Yes")
  contTableDnme[1,]$nonASM <- sum(featureDnme$ASM=="No")
  contTableDnme[2,]$ASM <- sum(complementaryDnme$ASM=="Yes")
  contTableDnme[2,]$nonASM <- sum(complementaryDnme$ASM=="No")
  dnmeFisher <- fisher.test(contTableDnme)
  
  # Enrichment of UC-HASM in feature
  featureUc <- featureData[featureData$Statistic=="UC"]
  complementaryUc <- complementaryData[complementaryData$Statistic=="UC"]
  contTableUc <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableUc) <- c("Feature","Complementary")
  contTableUc[1,]$ASM <- sum(featureUc$ASM=="Yes")
  contTableUc[1,]$nonASM <- sum(featureUc$ASM=="No")
  contTableUc[2,]$ASM <- sum(complementaryUc$ASM=="Yes")
  contTableUc[2,]$nonASM <- sum(complementaryUc$ASM=="No")
  ucFisher <- fisher.test(contTableUc)
  
  # Return list of Fisher's test
  return(list(dmmlFisher,dnmeFisher,ucFisher))
  
}

imprintingGtexEnrichment <- function(gtexDir,cpgDir,enhancerDir,cpelGRs,testStat="dMML"){
  
  # Read in GTEx file
  gtex <- read.csv(file=paste(gtexDir,"GTEx Portal.csv",sep=""))
  gtex <- gtex[gtex$Classification=="imprinted",]
  
  # Get promoters
  genomicFeatures <- getGeneralFeats(cpgDir,enhancerDir)
  proms <- genomicFeatures[["promoter"]]
  rm(list=c("genomicFeatures"))
  
  # Get overlaps
  olaps <- findOverlaps(cpelGRs,proms,type="within",select="all")
  cpelGRs$geneSymbol <- NA
  cpelGRs[queryHits(olaps)]$geneSymbol <- proms[subjectHits(olaps)]$gene_name
  rm(list=c("olaps","proms"))
  
  # Get haps that overlap with a promoter with testStat data
  promHaps <- cpelGRs[(!is.na(cpelGRs$geneSymbol))&(cpelGRs$GTEx!="")&(cpelGRs$Statistic==testStat),]
  
  # Cross imprinting information from GTEx with CpelAsm
  asmInProms <- GRanges()
  contTable <- data.frame(ASM=c(0,0),nonASM=c(0,0))
  rownames(contTable) <- c("Imprinted","Non-imprinted")
  for(i in 1:length(promHaps)) {
    
    # Check if promoter is imprinted at all
    if(promHaps[i]$geneSymbol %in% gtex$Gene.Symbol){
      
      # Promoter is imprinted for at least a tissue. Check what tissues
      impTissues <- gtex[gtex$Gene.Symbol==promHaps[i]$geneSymbol,]$Tissue.Name
      
      # Check if pair (gene,tissue) is imprinted
      if(promHaps[i]$GTEx %in% impTissues){
        # Promoter is imprinted for the tissue
        if(promHaps[i]$ASM=="Yes"){
          contTable$ASM[1] <- contTable$ASM[1]+1
          asmInProms <- append(asmInProms,promHaps[i])
        } else{
          contTable$nonASM[1] <- contTable$nonASM[1]+1
        }
      } else {
        # Promoter is not imprinted for the tissue
        if(promHaps[i]$ASM=="Yes"){
          contTable$ASM[2] <- contTable$ASM[2]+1
        } else{
          contTable$nonASM[2] <- contTable$nonASM[2]+1
        }
      } 
    } else {
      # Promoter is not imprinted for any tissue
      if(promHaps[i]$ASM=="Yes"){
        contTable$ASM[2] <- contTable$ASM[2]+1
      } else {
        contTable$nonASM[2] <- contTable$nonASM[2]+1
      }
    }
    
  }
  
  # Return list with results and ASM events in promoters
  return(list(test=fisher.test(contTable),GR=asmInProms))
  
}


# Function to compute within subject co-occurrence statistic
testStat <- function(inTable){
  
  # Return proportion of ASM>1 given there is data in more than one tissue and at least one ASM event
  sum(inTable$Yes>1)/sum((inTable$Yes>0)&(inTable$Data>1))
  
}

# Function to to permute subject ASM table
permuteTable <- function(inTable){
  
  # Copy input table
  outTable <- inTable
  outTable$Yes <- 0 
  
  # Add Sigs until it equals the number of total number of sigs
  while(sum(outTable$Yes)<sum(inTable$Yes)){
    ind <- sample(1:nrow(outTable),1)
    outTable$Yes[ind] <- min(outTable$Yes[ind]+1,outTable$Data[ind])
  }
  
  # return permuted table
  return(outTable)
  
}

# Function that returns null statistic for co-occurrence permutation test
parallelExecutionPermtest <- function(cluster,subject_table,perms){
  
  parSapply(cluster,1:perms,function(i,...){testStat(permuteTable(subject_table))})
  
}

# Function for co-occurrence permutation test
cooccurrenceTest <- function(cpelGRs,testStatistic,perms=10,ncpus=4){
  
  # Initialize output table
  uniqueSubjectLabels <- unique(cpelGRs$Subject)
  cooccurrenceTable <- data.frame(Subject=uniqueSubjectLabels,permTestStat=NA,pValue=NA)
  
  for(i in 1:length(uniqueSubjectLabels)) {
    
    # Do test for subject 
    subject <- uniqueSubjectLabels[i]
    GRs_subject <- cpelGRs[(cpelGRs$Subject==subject)&(cpelGRs$Statistic==testStatistic)]
    
    # Skip if only one tissue
    if(subject %in% c("HUES64","STL001","STL002","STL003")){
      # Print subject label
      print(subject)
    } else {  
      next
    }
    
    # Prepare table
    GR_haps <- unique(GRs_subject[,c()])
    GR_haps$Haplotype <- paste("Hap-",seq(1:length(GR_haps)),sep="")
    olaps <- findOverlaps(GRs_subject,GR_haps,type="equal",ignore.strand=TRUE,select="all")
    GRs_subject$Haplotype <- NA
    GRs_subject[queryHits(olaps)]$Haplotype <- GR_haps[subjectHits(olaps)]$Haplotype
    subject_table <- xtabs(~ASM+Haplotype,data=GRs_subject)
    subject_table <- as.data.frame.matrix(t(subject_table))
    subject_table$Data <- subject_table$No+subject_table$Yes
    subject_table$No <- NULL
    
    # Compute statistic for permutation test
    obsStat <- testStat(subject_table)
    
    # put objects in place that might be needed for the code
    cluster <- makeCluster(ncpus)
    clusterExport(cluster,c("permuteTable","testStat"))
    
    # Set a different seed on each member of the cluster (just in case)
    clusterSetRNGStream(cluster)
    
    #... then parallel replicate...
    nullStats <- parallelExecutionPermtest(cluster,subject_table,perms)
    
    # stop cluster
    stopCluster(cluster)
    
    # Compute p-value permutation test
    pVal <- sum(nullStats>=obsStat)/length(nullStats)
    
    # Store result permutation test
    cooccurrenceTable[i,2] <- obsStat
    cooccurrenceTable[i,3] <- pVal
    
  }
  
  # Return table
  return(cooccurrenceTable)
  
}

# Function to plot a GR using Gviz
plotGR <- function(CpGdir,enhancerDir,GR,startHight,highSize=500,reverseStrand=FALSE,chr="chr11",lim=c(2010000,2022500)){
  
  # Get genome
  gen <- "hg19"

  # GRanges to intersect with and keep the relevant data
  windowGR <- GRanges(seqnames=chr,ranges=IRanges(start=lim[1],end=lim[2]),strand="*")
  
  # Subset GR
  GR <- subsetByOverlaps(GR,windowGR)

  # Create data tracks
  dmmlTrack <- DataTrack(GR[GR$Statistic=="dMML",c("Value")],name="dMML")
  dnmeTrack <- DataTrack(GR[GR$Statistic=="dNME",c("Value")],name="dNME")
  ucTrack <- DataTrack(GR[GR$Statistic=="UC",c("Value")],name="UC")
  
  # Gene track 1
  bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="grch37.ensembl.org",
                path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
  biomTrack <- BiomartGeneRegionTrack(genome=gen,chromosome=chr,start=lim[1],end=lim[2],name="ENSEMBL",biomart=bm)

  # Add CpG island annotation track
  genomicFeatures <- getGeneralFeats(CpGdir,enhancerDir)
  cpgIslands <- genomicFeatures[["CpG island"]]
  cpgIslands <- subsetByOverlaps(cpgIslands,windowGR)
  islandTrack <- AnnotationTrack(cpgIslands,name="CpG islands")
  
  # Chromosome information tracks
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome=gen,chromosome=chr)
  
  # Highlight
  ht <- HighlightTrack(trackList=list(dmmlTrack,dnmeTrack,ucTrack),start=startHight,width=highSize,chromosome=chr)
  
  # Return plot
  plotTracks(list(itrack,gtrack,biomTrack,islandTrack,ht),from=lim[1],to=lim[2],
             transcriptAnnotation="symbol",type=c("gradient"),stacking="squish",reverseStrand=reverseStrand,
             collapseTranscripts = "meta")
  
}

# Function to plot PCA
generatePcaPlot <- function(inDir,pcaDataDf,testStatistic){
  
  # Get matrix
  pcaDataDfUc <- pcaDataDf[pcaDataDf$Statistic==testStatistic,][,c(1,2,4,5)]
  sigHapNames <- unique(pcaDataDfUc$HapName[pcaDataDfUc$ASM=="Yes"])
  pcaDataDfUc <- pcaDataDfUc[pcaDataDfUc$HapName %in% sigHapNames,]
  pcaDataDfUc$Subject <- gsub('.*- ', '',pcaDataDfUc$Sample)
  pcaDataDfUc$Tissue <- gsub(' -.*', '',pcaDataDfUc$Sample)
  pcaDataDfUc <- subset(pcaDataDfUc,pcaDataDfUc$Subject %in% c("STL001","STL002","STL003"))
  pcaMatrix <- dcast(data=pcaDataDfUc,formula=HapName~Sample,fun.aggregate=mean,value.var="Value")
  
  # Keep only haplotypes with enough data
  pcaMatrix <- pcaMatrix[apply(pcaMatrix, 1, function(x) sum(!is.na(x)))>=20,]
  pcaMatrix <- pcaMatrix[,apply(pcaMatrix, 2, function(x) sum(!is.na(x)))>50]
  rownames(pcaMatrix) <- pcaMatrix$HapName
  pcaMatrix <- pcaMatrix[,-1]
  
  # PCA
  res.comp <- imputePCA(pcaMatrix,ncp=2,method='EM')
  res.pca <- prcomp(res.comp$completeObs,center=TRUE,scale.=TRUE)
  loadsPpca <- as.data.frame(res.pca$rotation[,1:2])
  loadsPpca$Subject <- gsub('.*- ', '', rownames(loadsPpca))
  loadsPpca$Tissue <- gsub(' -.*', '', rownames(loadsPpca))
  ggplot(loadsPpca,aes(x=PC1,y=PC2,color=Subject)) + geom_point() + 
    theme(panel.border=element_blank(),panel.background=element_blank(),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          axis.line=element_line(colour="black"))
  ggsave(paste(inDir,"PCA ",testStatistic,".eps",sep=""),plot=last_plot(),device="eps",width=4,height=2.5,dpi=300)
  
}

# Function to extract genotype (needs VariantAnnotation)
extractGenotype <- function(vcfDir,rngs,sub){
  
  # Subset on genome coordinates: 'file' must have a Tabix index
  param <- ScanVcfParam(which=rngs)
  vcf <- readVcf(file=paste(vcfDir,sub,".phased.vcf.gz",sep=""),"hg19",param)
  
  # Return ranges with variants
  return(rowRanges(vcf))
  
}
readSubVcf <- function(inDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Get entire genome GR object
  genomeGr <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  genomeGr <- setGenomeLengths(genomeGr)
  
  # Read VCF file and 
  vcfSNP <- readVcf(file=paste(inDir,sub,".phased.vcf.gz",sep=""),"hg19",genomeGr)
  
  # Return GR with all haplotypes
  return(rowRanges(vcfSNP))
  
}
# Function to get variant frequency table
getSubjectVariantTableDiff <- function(vcfDir,testStat,sub,asm){
  
  # Obtain temporary GR for significant statistics and subject
  tmpGR <- cpelGRs[(cpelGRs$Statistic==testStat)&(cpelGRs$Subject==sub)&(cpelGRs$ASM==asm)]
  # Extract variants from VCF file corresponding to subject
  varsDiff <- extractGenotype(vcfDir,tmpGR,sub)
  # Extract
  variants <- paste(as.character(varsDiff$REF),as.character(unlist(varsDiff$ALT)),sep="-")
  # Combine same variants (ALT/REF -> REF/ALT)
  variants[variants %in% c("A-C","C-A")] <- "A-C"
  variants[variants %in% c("A-G","G-A")] <- "A-G"
  variants[variants %in% c("A-T","T-A")] <- "A-T"
  variants[variants %in% c("C-G","G-C")] <- "C-G"
  variants[variants %in% c("C-T","T-C")] <- "C-T"
  variants[variants %in% c("G-T","T-G")] <- "G-T"
  # Create table of frequency counts
  variants <- as.data.frame(table(variants))
  colnames(variants) <- c("Variant","Frequency")
  # Add ASM
  variants$ASM <- asm
  # Add subject
  variants$Subject <- sub  
  # Return variant table
  return(variants)
  
}

# Function that returns frequency of each SNP for each sample and for ASM/non-ASM 
getAllVariantTablesDiff <- function(vcfDir,cpelGRs,testStat){
  
  # Loop over each subject
  singleVariantsTable <- data.frame(Variant=c(),Frequency=c(),ASM=c(),Subject=c())
  for(sub in unique(cpelGRs$Subject)){
    # Print subject
    print(sub)
    
    # Get variant frequency tables
    print("ASM Haplotypes")
    variantsYes <- getSubjectVariantTableDiff(vcfDir,testStat,sub,"Yes")
    print("Non-ASM haplotypes")
    variantsNo <- getSubjectVariantTableDiff(vcfDir,testStat,sub,"No")
    
    # Append to output table
    singleVariantsTable <- rbind(singleVariantsTable,variantsYes,variantsNo)
  }
  
  # Return complete table
  return(singleVariantsTable)
  
}

# Function to assign genome
assignGenomeFile <- function(row) {
  
  genome = NA
  if((row[1]=="REF")&(row[2] %in% c("0|1","0/1"))){
    genome = 1
  } else if((row[1]=="ALT")&(row[2] %in% c("0|1","0/1"))){
    genome = 2
  } else if((row[1]=="REF")&(row[2] %in% c("1|0","1/0"))) {
    genome = 2
  } else if((row[1]=="ALT")&(row[2] %in% c("1|0","1/0"))){
    genome = 1
  } else {
    # nothing
  }
  
  return(genome)
}

# Function that returns bool if variant results in a CpG site
hetCpgSite <- function(row) {
  
  # Check if we have C or G, otherwise return false
  if(row[3] %in% c("C","G")){
    # continue
  } else {
    return(FALSE)
  }
  
  # Initialize output
  hetCpg <- FALSE
  context <- toString(getSeq(Hsapiens,row[1],start=as.numeric(row[2])-1,end=as.numeric(row[2])+1,strand="+"))
  if((row[3]=="C")&(substr(context,3,3)=="G")){
    hetCpg <- TRUE
  } else if((row[3]=="G")&(substr(context,1,1)=="C")){
    hetCpg <- TRUE
  } else {
    # nothing
  }
  
  # Return binary vector
  return(hetCpg)
  
}
het_CpG_df<-function(var_gr){
  #Get plus one location
  plus_loc=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,1)))
  minus_loc=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,-1)))
  #get dinucleotide for ref, alt, plus and minus
  var_gr$REF_plus=paste_nucleotide(var_gr$REF,plus_loc,'plus')
  var_gr$REF_minus=paste_nucleotide(var_gr$REF,minus_loc,'minus')
  var_gr$ALT_plus=paste_nucleotide(var_gr$ALT,plus_loc,'plus')
  var_gr$ALT_minus=paste_nucleotide(var_gr$ALT,minus_loc,'minus')
  #get trinucleotide
  var_gr$REF_tri=paste_trinucleotide(minus_loc,var_gr$REF,plus_loc)
  var_gr$ALT_tri=paste_trinucleotide(minus_loc,var_gr$ALT,plus_loc)
  #check if heterogygouze: note rowSum =2 have trinucleotide form CGG with ref =G alt =C
  var_gr$npmCG=rowSums(as.data.frame(var_gr)[,c('REF_plus','REF_minus','ALT_plus','ALT_minus')]=='CG')
  var_gr$HetCpg=var_gr$npmCG>0
  return(var_gr)
}
paste_nucleotide<-function(ref,seq,direction){
  df=data.frame(ref=ref,seq=unlist(seq))
  if(direction=='plus'){return(paste(df$ref,df$seq,sep=''))}
  else if(direction=='minus'){return(paste(df$seq,df$ref,sep=''))}
  else{print('wrong direction')}
}
paste_trinucleotide<-function(minus,ref,plus){
  df=data.frame(minus=unlist(minus),ref=ref,plus=unlist(plus))
  return(paste(df$minus,df$ref,df$plus,sep=''))
}
# Function to extract genotype information for genome1 and genome2 (needs VariantAnnotation)
extractGenotypePerAllele <- function(vcfDir,rngs,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Subset on genome coordinates: 'file' must have a Tabix index
  param <- ScanVcfParam(which=rngs)
  vcf <- unique(readVcf(file=paste(vcfDir,sub,".phased.vcf.gz",sep=""),"hg19",param))
  
  # Extract genome-1/2 information
  gt <- as.vector(geno(vcf)$GT)
  
  # Extract Granges information
  vcf <- rowRanges(vcf)
  
  # Add unique SNP id
  vcf$GT <- gt
  vcf$snpId <- paste(sub,seq(1:length(vcf)),sep="-")
  
  # Keep only relevant variables
  vcf <- vcf[,c("REF","ALT","GT","snpId")]
  vcf$REF <- as.character(vcf$REF)
  vcf$ALT <- as.character(unlist(vcf$ALT))
  names(vcf)=NULL
  # Delete labels
  vcf=het_CpG_df(vcf)
  # Convert to data frame
  vcfDf <- as.data.frame(vcf)
  colnames(vcfDf) <- c("chr","pos","end","width","strand","REF","ALT","GT","snpId",'REF_plus',
                       'REF_minus','ALT_plus','ALT_minus','REF_tri','ALT_tri','npmCG',"HetCpg")
  #keep <- names(vcfDf) %in% c("chr","pos","REF","ALT","GT","snpId",'HetCpg')
  vcfDf <- vcfDf[,-which(colnames(vcfDf)%in%c('width','strand'))]
  
  # Melt: origin = ref/alt, variant= ATCG
  vcfDf <- melt(data=vcfDf,id.vars=c("chr","pos","GT","snpId",'REF_plus',
                                     'REF_minus','ALT_plus','ALT_minus','REF_tri','ALT_tri','npmCG','HetCpg'),
                measure.vars=c("REF","ALT"),
                variable.name="Origin",value.name="Variant")
  #reshape the vcfDf, ALT to ALT, REF to REF
  vcfDf$plus[vcfDf$Origin=='REF']=vcfDf$REF_plus[vcfDf$Origin=='REF']
  vcfDf$plus[vcfDf$Origin=='ALT']=vcfDf$ALT_plus[vcfDf$Origin=='ALT']
  vcfDf$ALT_plus=NULL
  vcfDf$REF_plus=NULL
  vcfDf$minus[vcfDf$Origin=='REF']=vcfDf$REF_minus[vcfDf$Origin=='REF']
  vcfDf$minus[vcfDf$Origin=='ALT']=vcfDf$ALT_minus[vcfDf$Origin=='ALT']
  vcfDf$REF_minus=NULL
  vcfDf$ALT_minus=NULL
  vcfDf$tri[vcfDf$Origin=='REF']=vcfDf$REF_tri[vcfDf$Origin=='REF']
  vcfDf$tri[vcfDf$Origin=='ALT']=vcfDf$ALT_tri[vcfDf$Origin=='ALT']
  vcfDf$REF_tri=NULL
  vcfDf$ALT_tri=NULL
  vcfDf$HetAllele=FALSE
  vcfDf$HetAllele[vcfDf$HetCpg]=TRUE
  #one is CG the other one is not
  vcfDf$HetAllele[vcfDf$plus=='CG'| vcfDf$minus=='CG']=FALSE
  #CGG and CCG special case
  #vcfDf$HetAllele[vcfDf$npmCG==2]=TRUE
  # Assign genome file: separate case 0|1 or 1|0 etc
  vcfDf$Genome <- apply(vcfDf[,c("Origin","GT")],1,assignGenomeFile)

  # Add bool about heterozygous CpG site
  #vcfDf$HetCpg <- apply(vcfDf[,c("chr","pos","Variant")],1,hetCpgSite)
    # vcfDf$HetCpg <- lapply(vcfDf[,c("chr","pos","Variant")],FUN=hetCpgSite)
  
  # Add bool about heterozygous CpG site
    # plan(multiprocess)
    # vcfDf$HetCpg <- future_lapply(vcfDf[,c("chr","pos","Variant")],FUN=hetCpgSite)
  
  # Convert it back to GRanges object
  vcf <- makeGRangesFromDataFrame(vcfDf,ignore.strand=TRUE,start.field="pos",
                                  end.field="pos",keep.extra.columns=TRUE)
  #olap=findOverlaps(vcf,rngs,type='within')
  # Return ranges with variants
  return(vcf)
  
}

# Function to extract genotype information for genome1 and genome2 (needs VariantAnnotation)
crossVcfEtGffEtCpelAsm <- function(vcfDir,inDir,cpelAllele,sub,tissue,chrsOfInterest=paste("chr",1:22,sep="")){
  cat(paste('Processing tissue:',tissue,'for subject',sub,'\n'))
  t1=proc.time()[3]
  # Get genotype for each allele overlapping a haplotype,#Here only look for heterozygous CpG on one allele not the other one, need a correction method
  cpelAllele=cpelAllele[(cpelAllele$Subject==sub)&(cpelAllele$Tissue==tissue)]#Sample, tissue specific 
  outGR <- extractGenotypePerAllele(vcfDir,cpelAllele,sub,chrsOfInterest=chrsOfInterest)
  
  gffSub <- readSubGff(inDir,sub,chrsOfInterest=chrsOfInterest)
  
  #  Assign Hap Id to GFF GR
  gffSub$HapId <- paste("Hap",seq(1:length(gffSub)),sep="-")
  
  # Find overlaps between VCF and GFF
  olaps <- findOverlaps(outGR,gffSub,type="within")
  outGR$N[queryHits(olaps)] <- gffSub$N[subjectHits(olaps)]
  outGR$HapId[queryHits(olaps)] <- gffSub$HapId[subjectHits(olaps)]
  
  # Cross resulting GR with MML of genome 1
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="1")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$MML1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$MML1[outGR$Genome=="2"] <- NA
  
  # Cross resulting GR with MML of genome 2
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="2")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$MML2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$MML2[outGR$Genome=="1"] <- NA
  
  # Consolidate MML1 and MML2 columns into single column
  outGR$MML <- rowSums(cbind(outGR$MML1,outGR$MML2), na.rm=T)
  outGR$MML1 <- outGR$MML2 <- NULL
  
  # Cross resulting GR with NME of genome 1
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="1")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$NME1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$NME1[outGR$Genome=="2"] <- NA
  
  # Cross resulting GR with NME of genome 2
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="2")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$NME2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$NME2[outGR$Genome=="1"] <- NA

  # Consolidate NME1 and NME2 columns into single column
  outGR$NME <- rowSums(cbind(outGR$NME1,outGR$NME2), na.rm=T)
  outGR$NME1 <- outGR$NME2 <- NULL
  cat(paste('Processing tissue:',tissue,'for subject',sub,'in',proc.time()[3]-t1,'\n'))
  outGR$Sample=paste(sub,tissue,sep='-')
  # Return
  return(outGR)
  
}

GR_allele_correct <-function(outGR,cpelAllele,sub,tissue){
  outGR$MML=NA
  outGR$NME=NA
  cpelAllele=cpelAllele[(cpelAllele$Subject==sub)&(cpelAllele$Tissue==tissue)]#Sample, tissue specific 
  # Cross resulting GR with MML of genome 1
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="1")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$MML1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$MML1[outGR$Genome=="2"] <- NA
  
  # Cross resulting GR with MML of genome 2
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="2")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$MML2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$MML2[outGR$Genome=="1"] <- NA
  
  # Consolidate MML1 and MML2 columns into single column
  outGR$MML <- rowSums(cbind(outGR$MML1,outGR$MML2), na.rm=T)
  outGR$MML1 <- outGR$MML2 <- NULL
  
  # Cross resulting GR with NME of genome 1
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="1")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$NME1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$NME1[outGR$Genome=="2"] <- NA
  
  # Cross resulting GR with NME of genome 2
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="2")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$NME2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$NME2[outGR$Genome=="1"] <- NA
  
  # Consolidate NME1 and NME2 columns into single column
  outGR$NME <- rowSums(cbind(outGR$NME1,outGR$NME2), na.rm=T)
  outGR$NME1 <- outGR$NME2 <- NULL
  outGR$Sample=paste(sub,tissue,sep='-')
  # Return
  return(outGR)
}
add_ASM<-function(allele_GR,ASM_GR){
  allele_GR_out=GRanges()
  for (sp in unique(allele_GR$Sample)){
    ASM_sub=ASM_GR[ASM_GR$Sample==sp]
    allele_sub=allele_GR[allele_GR$Sample==sp]
    allele_sub$ASM=NA
    olap=findOverlaps(allele_sub,ASM_sub)
    allele_sub$ASM[queryHits(olap)]=ASM_sub$ASM[subjectHits(olap)]
    allele_GR_out=c(allele_GR_out,allele_sub)
  }
  return(allele_GR_out) 
}
allele_diff<-function(allele_gr_in){
  out_gr=GRanges()
  for(sample in unique(allele_gr_in$Sample)){
    cat('Processing:',sample,'\n')
    genome1=sort(allele_gr_in[allele_gr_in$Genome==1 & allele_gr_in$Sample==sample])
    
    genome2=sort(allele_gr_in[allele_gr_in$Genome==2 & allele_gr_in$Sample==sample])
    gr=granges(genome1)
    gr$diff=genome1$Value-genome2$Value
    gr$Sample=sample
    gr$Statistic=unique(genome1$Statistic)
    out_gr=c(out_gr,gr)
  }
  return(out_gr)
}
diff_mean_sample<-function(diff_gr){
  out=data.frame(sample=as.character(unique(diff_gr$Sample)),diff=NA)
  for (sp in out$sample){
    out$diff[out$sample==sp]=mean(diff_gr$diff[diff_gr$Sample==sp])
    
  }
  
  return(out)
}
hetCpG_check<-function(varSub,GRASM,stat,ASM=TRUE){
  GRASM=GRASM[GRASM$Statistic==paste('d',stat,sep='')]
  if(ASM){GRASM=GRASM[GRASM$ASM=='Yes']}
  het_CpG=varSub[varSub$HetCpg]
  het_CpG=subsetByOverlaps(het_CpG,GRASM,type='within')
  het_CpG=het_CpG[het_CpG$npmCG!=2]
  het_CpG_homallele=het_CpG[!het_CpG$HetAllele]
  het_CpG_hetallele=het_CpG[het_CpG$HetAllele]
  het_CpG_homallele=sort(het_CpG_homallele)
  het_CpG_hetallele=sort(het_CpG_hetallele)
  cat('Identical check between hom and het:',identical(het_CpG_hetallele$snpId,het_CpG_homallele$snpId),'\n')
  out_gr=granges(het_CpG_homallele)
  out_gr$nonCpG=elementMetadata(het_CpG_hetallele)[,stat]
  out_gr$cpg=elementMetadata(het_CpG_homallele)[,stat]
  out_gr$diff= out_gr$nonCpG-out_gr$cpg
  out_gr$snpId=het_CpG_hetallele$snpId
  out_gr$Sample=unique(varSub$Sample)
  out_gr$stat=stat
  return(out_gr)
}
hetCpG_check_subject<-function(varSub,GRASM,stat,ASM=TRUE,sub){
  out=GRanges()
  for (i in 1:length(varSub)){
    var_sub_tissue=varSub[[i]]
    tissue=unique(gsub(paste(unique(GRASM$Subject),'-',sep=''),'',var_sub_tissue$Sample))
    cat('Checking:',tissue,'\n')
    GRASM_sub=GRASM[GRASM$Tissue==tissue]
    
    out=c(out,hetCpG_check(var_sub_tissue,GRASM_sub,stat,ASM=ASM))
  }
  out$sub=sub
  return(out)
}

GR_allele_correct_sub<-function(varSub,GR_allele,sub){
  out=GRangesList()
  for (i in 1:length(varSub)){
    var_sub_tissue=varSub[[i]]
    tissue=unique(gsub(paste(sub,'-',sep=''),'',var_sub_tissue$Sample))
    cat('Correcting:',tissue,'\n')
    out[[i]]=GR_allele_correct( var_sub_tissue,GR_allele,sub,tissue)
  }
  return(out)
}
calc_OR<-function(feature,GR_all,GR_hit,CpG){
  olap_gr=findOverlaps(GR_all,GR_hit,type='equal')
  GR_not_hit=GR_all[-queryHits(olap_gr)]
  CpG=subsetByOverlaps(CpG,GR_all,type='within')
  GR_hit=subsetByOverlaps(CpG,GR_hit,type='within')
  GR_not_hit=subsetByOverlaps(CpG,GR_not_hit,type='within')
  GR_hit_feature=length(subsetByOverlaps(GR_hit,feature))
  
  GR_hit_notfeature=length(GR_hit)-GR_hit_feature
  GR_not_hit_feature=length(subsetByOverlaps(GR_not_hit,feature))
  GR_not_hit_not_feature=length(GR_not_hit)-GR_not_hit_feature
  cont_table=matrix(c(GR_hit_feature,GR_hit_notfeature,GR_not_hit_feature,GR_not_hit_not_feature),nrow=2)
  colnames(cont_table)=c('Hit','Not Hit')
  rownames(cont_table)=c('Feature','Not feature')
  print(fisher.test(cont_table))
  return(cont_table)
}

gff_to_CpG_loop<-function(gff_gr_in,chrsOfInterest=paste("chr",1:22,sep="")){
  gff_gr_out=GRanges()
  for(chr in chrsOfInterest){
    gff_gr_chr=gff_gr_in[seqnames(gff_gr_in)==chr]
    for(sub in unique(gff_gr_chr$Subject)){gff_gr_out=c(gff_gr_out,gff_to_CpG(gff_gr_chr[gff_gr_chr$Subject==sub],chr,sub))}
  }
  return(setGenomeLengths(gff_gr_out))
}

gff_to_CpG<-function(gff_gr_in,chr,subject){
  loc=unique(do.call('c',lapply(gff_gr_in$CpGs,get_CpG_gff)))
  gr_df=data.frame(seqnames=chr,start=loc,end=loc,Subject=subject)
  return(makeGRangesFromDataFrame(gr_df,keep.extra.columns = T))
}
get_CpG_gff<-function(loc){
  loc=gsub('\\[','',loc)
  loc=gsub('\\]','',loc)
  loc=gsub(' ','',loc)
  loc=as.numeric(loc)
  return(loc)
}
