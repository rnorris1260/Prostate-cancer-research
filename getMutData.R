
if (nrow(BRCA2) != 0){ rbind(HMF_PCa_BRCA2, BRCA2)}

files <- list.files("~/HMF/somatic/DR-252-update1")
for (file in files){
   exome <- read.csv(file)
   BRCA2 <- subset(exome, exome$gene == "BRCA2")
   rbind(HMF_PCa_BRCA2, BRCA2)
 }


for (file in files){
  exome <- read.csv(file)
  if ("BRCA2" %in% exome$gene){
    print(file)
    exome <- read.csv(file)
    BRCA2 <- data.frame(subset(exome, exome$gene == "BRCA2"))
    print(nrow(BRCA2))
    HMF_PCa_BRCA2 <- rbind(HMF_PCa_BRCA2, BRCA2)
    print("mutations added to df")
  }
}

### Driver

read.delim(paste(file,"/",file,".isf.gene_data.csv", sep=""), sep=",")

setwd("~/HMF/somatic/DR-252-update1")

n = 1
for (file in files){
  driver <- read.delim(paste(file,"/purple/",file,".purple.driver.catalog.somatic.tsv",sep=""), sep="\t")
  genes <- driver$gene
  if ("BRCA2" %in% genes == T){
    print(n)
    print(file)
    #write.csv(driver, paste("BRCA2/",file,".csv",sep=""))
    # print(nrow(BRCA2))
    # HMF_PCa_BRCA2dr <- rbind(HMF_PCa_BRCA2dr, BRCA2)
    # print("mutations added to df")
    n = n+1
  }
}

setwd("~/HMF/exomeAnnotation/")

HMF_exome1 <- read.csv("ACTN01020017T.exome.annot.csv")
HMF_PCa_BRCA2 <- data.frame(subset(HMF_exome1, HMF_exome1$gene == "BRCA2"))

## get list of patients with and without BRCA2 alterations
BRCA2_files <- list.files("~/HMF/BRCA2/")
## remove .csv
BRCA2_files <- gsub('.{0,4}$', '', BRCA2_files)
non_BRCA2_files <- setdiff(files, BRCA2_files)

HMF_PCa_BRCA2dr <- data.frame(subset(HMF_exome1, HMF_exome1$gene == "BRCA2"))

for (file in files){
  df <- read.csv(paste("BRCA2/",file,sep=""))
  brca2 <- data.frame(subset(df, df$gene == "BRCA2"))
  HMF_PCa_BRCA2dr <- rbind(HMF_PCa_BRCA2dr, brca2)
}

driver2 <- read.delim("CPCT02010399T/purple/CPCT02010399T.purple.driver.catalog.somatic.tsv",sep="\t")
BRCA2dr2 <- subset(driver2, driver2$gene == "BRCA2")

driver1 <- read.delim("ACTN01020021T/purple/ACTN01020021T.purple.driver.catalog.somatic.tsv", sep="\t")

genes1 <- driver1$gene                      
"BRCA2" %in% driver$gene

setwd("~/HMF/somatic/DR-252-update1")

### empty df
HMF_PCa_empty <- HMF_PCa_BRCA2
new_col <- "Patient_ID"
HMF_PCa_empty[1,] <- NA
HMF_PCa_empty[, 'Patient_ID'] <- NA
HMF_PCa_empty <- HMF_PCa_empty[0,]
HMF_PCa_empty <- HMF_PCa_empty[,c(18,1:17)]

getMutFiles <- function(gene_ID){
  HMF_PCa_gene <- HMF_PCa_empty
  files <- list.files("~/HMF/somatic/DR-252-update1")
  for (file in files){
    driver <- read.delim(paste("~/HMF/somatic/DR-252-update1/",file,"/purple/",file,".purple.driver.catalog.somatic.tsv",sep=""), sep="\t")
    genes <- driver$gene
    if (gene_ID %in% genes == T){
      print(file)
      write.csv(driver, paste("~/HMF/",gene_ID,"/",file,".csv",sep=""))
    }
  }
  files <- list.files(paste("~/HMF/",gene_ID,sep=""))
  for (file in files){
    df <- read.csv(paste("~/HMF/",gene_ID,"/",file,sep=""))
    gene_data <- data.frame(subset(df, df$gene == gene_ID))
    ID <- substring(file, 1, nchar(file) -4)
    gene_data$Patient_ID <- ID
    gene_data <- gene_data[,c(19,2:18)]
    HMF_PCa_gene <- rbind(HMF_PCa_gene, gene_data)
  }
  return(HMF_PCa_gene) 
}

BRCA2 <- getMutFiles("BRCA2")
KMT2D <- getMutFiles("KMT2D")
APC <- getMutFiles("APC")
CTNNB1 <- getMutFiles("CTNNB1")


exome_annot_1 <- read.csv("~/HMF/exomeAnnotation/ACTN01020017T.exome.annot.csv")
exome_annot_empty <- exome_annot_1[0,]

getDrivers <- function(gene_ID){
  drivers_EA <- exome_annot_empty
  files <- list.files(paste("~/HMF/",gene_ID,sep=""))
  for (file in files){
    ID <- substring(file, 1, nchar(file) -4)
    print(ID)
    df <- read.csv(paste("~/HMF/",gene_ID,"/",file,sep=""))
    gene_data <- data.frame(subset(df, df$gene == gene_ID))
    tc <- unique(gene_data$transcript)
    print(tc)
    exome_annot <- read.csv(paste("~/HMF/exomeAnnotation/",ID,".exome.annot.csv",sep=""))
    exome_annot_tc <- subset(exome_annot, exome_annot$transcript == tc)
      if (nrow(exome_annot_tc) > 1){
    exome_annot_tc <- subset(exome_annot_tc, exome_annot_tc$annotation != "intron_variant")
      }
    print(nrow(exome_annot_tc))
    print(unique(exome_annot_tc$annotation))
    drivers_EA <- rbind(drivers_EA, exome_annot_tc)
  }
  return(drivers_EA)
}

df_test <- read.csv(paste("~/HMF/CTNNB1/",ID_test,".csv",sep=""))
gene_data_test <- data.frame(subset(df_test, df_test$gene == "CTNNB1"))
transcript <- unique(gene_data_test$transcript)

ID_test <-"DRUP01070222T"
test_exome_annot <- read.csv(paste("~/HMF/exomeAnnotation/",ID_test,".exome.annot.csv",sep=""))
test_exome_annot_sub <- subset(test_exome_annot, test_exome_annot$transcript == "ENST00000349496")


APC_drivers_EA <- getDrivers("APC")
CTNNB1_drivers_EA <- getDrivers("CTNNB1")

write.csv(APC_drivers_EA, "APC_driver_mutations.csv")
write.csv(CTNNB1_drivers_EA, "CTNNB1_driver_mutations.csv")


BRCA2_files <- list.files("~/HMF/BRCA2/")

# remove ".csv" from filenames
for (i in 1:length(BRCA2_files)){
  BRCA2_files[i] = substr(BRCA2_files[i], 1, nchar(BRCA2_files[i])-4)
}

All_files <- list.files("~/HMF/somatic/DR-252-update1")
`%nin%` = Negate(`%in%`)

non_BRCA2_files <- setdiff(All_files, BRCA2_files) 




### gene expression

APC_drivers_edit <- read.csv("APC_driver_mutations.csv")
CTNNB1_drivers_edit <- read.csv("CTNNB1_driver_mutations.csv")

# patient lists

APC_patients <- unique(APC_drivers_edit$sampleId)
CTNNB1_patients <- unique(CTNNB1_drivers_edit$sampleId)

APC_CTNNB1_patients <- append(APC_patients, CTNNB1_patients)
APC_CTNNB1_patients <- unique(APC_CTNNB1_patients)  # 2 patients have mutations in both AP and CTNNB1
non_APC_CTNNB1 <- setdiff(All_files, APC_CTNNB1_patients)

gene_ex1 <- read.csv("~/HMF/isofox/data_isofox/ACTN01020017T/ACTN01020017T.isf.gene_data.csv")
gene_ex_empty <- gene_ex1[0,]

gene_ex_empty[1,] <- NA
gene_ex_empty[, 'Patient_ID'] <- NA
gene_ex_empty <- gene_ex_empty[0,]

get_geneEx <- function(patients, gene_ID){
  geneEx_df <- gene_ex_empty
  for (pat in patients){
    if (file.exists(paste("~/HMF/isofox/data_isofox/",pat,"/",pat,".isf.gene_data.csv", sep="")) == T){ 
    df <- read.csv(paste("~/HMF/isofox/data_isofox/",pat,"/",pat,".isf.gene_data.csv", sep=""))
    gene_df <- subset(df, df$GeneName == gene_ID)
    gene_df$Patient_ID <- pat
    geneEx_df <- rbind(geneEx_df, gene_df)
    }}
  geneEx_df <- geneEx_df[, c(14,2,8:13)]
  return(geneEx_df)
}

mut_AXIN2 <- get_geneEx(APC_CTNNB1_patients, "AXIN2")  
mut_NKD1 <- get_geneEx(APC_CTNNB1_patients, "NKD1")
mut_RNF43 <- get_geneEx(APC_CTNNB1_patients, "RNF43")
mut_ZNRF3 <- get_geneEx(APC_CTNNB1_patients, "ZNRF3")

NOmut_AXIN2 <- get_geneEx(non_APC_CTNNB1, "AXIN2")  
NOmut_NKD1 <- get_geneEx(non_APC_CTNNB1, "NKD1")
NOmut_RNF43 <- get_geneEx(non_APC_CTNNB1, "RNF43")
NOmut_ZNRF3 <- get_geneEx(non_APC_CTNNB1, "ZNRF3")

write.csv(mut_AXIN2, file = "mut_AXIN2.csv")
write.csv(mut_NKD1, file = "mut_NKD1.csv")
write.csv(mut_RNF43, file = "mut_RNF43.csv")
write.csv(mut_ZNRF3, file = "mut_ZNRF3.csv")

write.csv(NOmut_AXIN2, file = "NOmut_AXIN2.csv")
write.csv(NOmut_NKD1, file = "NOmut_NKD1.csv")
write.csv(NOmut_RNF43, file = "NOmut_RNF43.csv")
write.csv(NOmut_ZNRF3, file = "NOmut_ZNRF3.csv")






