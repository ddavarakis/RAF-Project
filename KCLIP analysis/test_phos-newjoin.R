suppressPackageStartupMessages({
  library(PhosR)
  library(dplyr)
  library(tidyverse)
  library(tibble)
  library(sva)
  library("readxl")
  library("writexl")
})
#set directory
#setwd("C:/Users/mtheo/Desktop/BRAF KCLIP 260221/Proccessing tables S1 - S6")
setwd("D:\\Documents\\AMaria\\Jimmy\\2nd Trimester\\Project\\RAF substrates 2023\\BRAF KCLIP 260221\\Analysis")

# Loading
# xlsx files
#filename = "Table S1_BRAFWT.xlsx" # Table S1_BRAFWT, Table S2_BRAFV600E, Table S3_BRAFWT, Table S4_BRAFK483M, Table S5_BRAFV600E, Table S6_BRAFV600EK483M

# read all excels from folder
file.list <- list.files(pattern='*.xlsx')

for (filename in file.list) {
  Proteomics <- read_excel(filename)
  names(Proteomics)<-str_replace_all(names(Proteomics), c(" " = "." , "," = "" ))
  #print(filename)
  if(filename == "Table S1_BRAFWT.xlsx"){
    #conditions ???? BRAFWT No ATPArN3	BRAFWT with ATPArN3	BRAFK483M No ATPArN3	BRAFK483M with ATPArN3	Bead control
    conditions <- data.frame(c(1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20), 
                             c(rep("BRAFWT No ATPArN3",4),
                               rep("BRAFWT with ATPArN3",4),
                               rep("BRAFK483M No ATPArN3",4),
                               rep("BRAFK483M with ATPArN3",4),
                               rep("Bead control",4)),
                             rep(1:4, 5),
                             rep(1,20))
  }
  if(filename == "Table S2_BRAFV600E.xlsx"){
    #conditions ???? BRAFV600E No ATPArN3	BRAFV600E with ATPArN3	BRAFV600EK483M No ATPArN3	BRAFV600EK483M with ATPArN3	Bead control
    
    c2 <- c(1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20)#+20
    conditions <- data.frame(c2, 
                             c(rep("BRAFV600E No ATPArN3",4),
                               rep("BRAFV600E with ATPArN3",4),
                               rep("BRAFV600EK483M No ATPArN3",4),
                               rep("BRAFV600EK483M with ATPArN3",4),
                               rep("Bead control",4)),
                             rep(1:4, 5),
                             rep(2,20))
  }
  
  #TAbles S3 - S6
  if(filename == "Table S3_BRAFWT.xlsx"){ 
    #conditions ???? BRAF WT without ATPArN3	BRAF WT with ATPArN3	Bead control
    c3 <- c(1,4,7,10,2,5,8,11,3,6,9,12)#+40
    conditions <- data.frame(c3, 
                             c(rep("BRAFWT No ATPArN3",4),
                               rep("BRAFWT with ATPArN3",4),
                               rep("Bead control",4)),
                             rep(1:4, 3),
                             rep(3,12)) # (4 trials , 3 conditions) 
  }
  if(filename == "Table S4_BRAFK483M.xlsx"){ 
    #conditions ???? BRAF K483M without ATPArN3	BRAF K483M with ATPArN3	Bead control
    #1,2,3,
    #4,5,6,
    #7,8,9,
    #10,11,12 
    c4 <- c(1,4,7,10,2,5,8,11,3,6,9,12)#+52
    conditions <- data.frame(c4,
                             c(rep("BRAFK483M No ATPArN3",4),
                               rep("BRAFK483M with ATPArN3",4),
                               rep("Bead control",4)),
                             rep(1:4, 3),
                             rep(4,12))
  }
  if(filename == "Table S5_BRAFV600E.xlsx"){ 
    #conditions ???? BRAFV600E No ATParN3	BRAFV600E with ATParN3	Bead control
    c5 <- c(1,4,7,10,2,5,8,11,3,6,9,12)#+64
    conditions <- data.frame(c5,
                             c(rep("BRAFV600E No ATPArN3",4),
                               rep("BRAFV600E with ATPArN3",4),
                               rep("Bead control",4)),
                             rep(1:4, 3),
                             rep(5,12))
  }
  if(filename == "Table S6_BRAFV600EK483M.xlsx"){ 
    #conditions ???? BRAFV600EK483M without ATPArN3	BRAFV600EK483M with ATPArN3	Bead control
    c6 <- c(1,4,7,10,2,5,8,11,3,6,9,12)#+76
    conditions <- data.frame(c6,
                             c(rep("BRAFV600EK483M No ATPArN3",4),
                               rep("BRAFV600EK483M with ATPArN3",4),
                               rep("Bead control",4)),
                             rep(1:4, 3),
                             rep(6,12))
  }
  colnames(conditions) <- c("ID","condition","trial","dataset")
  
  conditions$trial <- as.factor(conditions$trial)
  conditions$dataset <- as.factor(conditions$dataset)
  

  ## remove potential contaminants ####
  #Proteomics <- Proteomics[Proteomics$Reverse != "+" & Proteomics$Potential.contaminant != "+", ]
  if((filename == "Table S2_BRAFV600E.xlsx") || (filename == "Table S1_BRAFWT.xlsx") ){
    Proteomics <- Proteomics[is.na(Proteomics$Reverse)  & is.na(Proteomics$Potential.contaminant), ]
  }
  #remove NA gene names
  which(is.na(Proteomics$Gene.names) )
  Proteomics <- Proteomics[!is.na(Proteomics$Gene.names), ]
  

  ## formatting intensity ####
  patterns <- c("BRAF*", "Bead*")
  Proteomics.Intensity <- Proteomics[,grepl(paste(patterns, collapse="|"), colnames(Proteomics))]

  ## Setting up experimental design for phosphoproteomics ####
  if((filename == "Table S2_BRAFV600E.xlsx") || (filename == "Table S1_BRAFWT.xlsx")){
    experimental_design <- data.frame("ID" = sub("_.*", "", colnames(Proteomics.Intensity)), row.names  = colnames(Proteomics.Intensity))
    experimental_design$ID <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
    experimental_design <- merge(experimental_design, conditions, by = "ID", all.x = T)
  }else{
    experimental_design <- data.frame("ID" = sub("_.*", "", colnames(Proteomics.Intensity)), row.names  = colnames(Proteomics.Intensity))
    experimental_design$ID <- c(1,2,3,4,5,6,7,8,9,10,11,12)
    experimental_design <- merge(experimental_design, conditions, by = "ID", all.x = T)
    
  }

  ## Make Na those with zero intensity
  Proteomics.Intensity[Proteomics.Intensity == 0] <- NA

  
  if(filename == "Table S1_BRAFWT.xlsx"){
    S1_Proteomics.Intensity <- data.frame(Proteomics.Intensity)
    rownames(S1_Proteomics.Intensity) <- Proteomics$Gene.names
    S1_experimental_design <- data.frame(experimental_design)
    S1_conditions <- data.frame(conditions)
  }
  if(filename == "Table S2_BRAFV600E.xlsx"){
    #S2_Proteomics.Intensity_log <- data.frame(Proteomics.Intensity_log) 
    S2_Proteomics.Intensity <- data.frame(Proteomics.Intensity)
    rownames(S2_Proteomics.Intensity) <- Proteomics$Gene.names
    S2_experimental_design <- data.frame(experimental_design)
    S2_conditions <- data.frame(conditions)
  }
  if(filename == "Table S3_BRAFWT.xlsx"){
    #S3_Proteomics.Intensity_log <- data.frame(Proteomics.Intensity_log) 
    S3_Proteomics.Intensity <- data.frame(Proteomics.Intensity)
    rownames(S3_Proteomics.Intensity) <- Proteomics$Gene.names
    S3_experimental_design <- data.frame(experimental_design)
    S3_conditions <- data.frame(conditions)
  }
  if(filename == "Table S4_BRAFK483M.xlsx"){
    #S4_Proteomics.Intensity_log <- data.frame(Proteomics.Intensity_log) 
    S4_Proteomics.Intensity <- data.frame(Proteomics.Intensity)
    rownames(S4_Proteomics.Intensity) <- Proteomics$Gene.names
    S4_experimental_design <- data.frame(experimental_design)
    S4_conditions <- data.frame(conditions)
  }
  if(filename == "Table S5_BRAFV600E.xlsx"){
    #S5_Proteomics.Intensity_log <- data.frame(Proteomics.Intensity_log) 
    S5_Proteomics.Intensity <- data.frame(Proteomics.Intensity)
    rownames(S5_Proteomics.Intensity) <- Proteomics$Gene.names
    S5_experimental_design <- data.frame(experimental_design)
    S5_conditions <- data.frame(conditions)
  }
  if(filename == "Table S6_BRAFV600EK483M.xlsx"){
    #S6_Proteomics.Intensity_log <- data.frame(Proteomics.Intensity_log) 
    S6_Proteomics.Intensity <- data.frame(Proteomics.Intensity)
    rownames(S6_Proteomics.Intensity) <- Proteomics$Gene.names
    S6_experimental_design <- data.frame(experimental_design)
    S6_conditions <- data.frame(conditions)
  }
  #global_experimendal_desine = rbind(experimental_design,dataframe2)
  
}# End For
# correct the conditions
S2_experimental_design$ID <- as.numeric(S2_experimental_design$ID) +20
global_experimental_design = rbind(S1_experimental_design,S2_experimental_design)
S3_experimental_design$ID <- as.numeric(S3_experimental_design$ID) +40
global_experimental_design = rbind(global_experimental_design,S3_experimental_design)
S4_experimental_design$ID <- as.numeric(S4_experimental_design$ID) +52
global_experimental_design = rbind(global_experimental_design,S4_experimental_design)
S5_experimental_design$ID <- as.numeric(S5_experimental_design$ID) +64
global_experimental_design = rbind(global_experimental_design,S5_experimental_design)
S6_experimental_design$ID <- as.numeric(S6_experimental_design$ID) +76
global_experimental_design = rbind(global_experimental_design,S6_experimental_design)
# correct the conditions
S2_conditions$ID <- as.numeric(S2_conditions$ID) +20
global_conditions = rbind(S1_conditions,S2_conditions)
S3_conditions$ID <- as.numeric(S3_conditions$ID) +40
global_conditions = rbind(global_conditions,S3_conditions)
S4_conditions$ID <- as.numeric(S4_conditions$ID) +52
global_conditions = rbind(global_conditions,S4_conditions)
S5_conditions$ID <- as.numeric(S5_conditions$ID) +64
global_conditions = rbind(global_conditions,S5_conditions)
S6_conditions$ID <- as.numeric(S6_conditions$ID) +76
global_conditions = rbind(global_conditions,S6_conditions)

# Make row names separate column
S1_Proteomics.Intensity <- tibble::rownames_to_column(S1_Proteomics.Intensity, "Gene_name")
S2_Proteomics.Intensity <- tibble::rownames_to_column(S2_Proteomics.Intensity, "Gene_name")
S3_Proteomics.Intensity <- tibble::rownames_to_column(S3_Proteomics.Intensity, "Gene_name")
S4_Proteomics.Intensity <- tibble::rownames_to_column(S4_Proteomics.Intensity, "Gene_name")
S5_Proteomics.Intensity <- tibble::rownames_to_column(S5_Proteomics.Intensity, "Gene_name")
S6_Proteomics.Intensity <- tibble::rownames_to_column(S6_Proteomics.Intensity, "Gene_name")

# Check whether there are proteins in S2 that does not exist in S1
# There are 337 proteins in S2 that does not exist in S1
diff <- S2_Proteomics.Intensity %>% filter(!Gene_name %in% S1_Proteomics.Intensity$Gene_name)

#### join S1 + S2 ####
joined_df <- merge(x = S1_Proteomics.Intensity, y = S2_Proteomics.Intensity, 
                   by = "Gene_name", all=TRUE)

colnames(joined_df) <- gsub('x','_S1',names(joined_df))
colnames(joined_df) <- gsub('.y','_S2',names(joined_df))

#### join S3 ####
colnames(S3_Proteomics.Intensity) <- gsub('Bead','S3_Bead.control',names(S3_Proteomics.Intensity))
joined_df <- merge(x = joined_df, y = S3_Proteomics.Intensity, 
                   by = "Gene_name", all=TRUE)

#### join S4 ####
colnames(S4_Proteomics.Intensity) <- gsub('Bead','S4_Bead.control',names(S4_Proteomics.Intensity))

joined_df <- merge(x = joined_df, y = S4_Proteomics.Intensity, 
                   by = "Gene_name", all=TRUE)

#### join S5 ####
colnames(S5_Proteomics.Intensity) <- gsub('Bead','S5_Bead.control',names(S5_Proteomics.Intensity))

joined_df <- merge(x = joined_df, y = S5_Proteomics.Intensity, 
                   by = "Gene_name", all=TRUE)

#### join S6 ####
colnames(S6_Proteomics.Intensity) <- gsub('Bead','S6_Bead.control',names(S6_Proteomics.Intensity))

joined_df <- merge(x = joined_df, y = S6_Proteomics.Intensity, 
                   by = "Gene_name", all=TRUE)

write_xlsx(joined_df,"artifacts\\joined_df.xlsx")
write_xlsx(global_experimental_design,"artifacts\\global_experimental_design.xlsx")
global.Proteomics.Intensity <- joined_df

# Drop Gene_name column
global.Proteomics.Intensity <- global.Proteomics.Intensity[,!names(global.Proteomics.Intensity) %in% c("Gene_name")]

# Rename column names
newcolnames <- as.data.frame(t(global_experimental_design['condition']))
global.Proteomics.Intensity <- setNames(global.Proteomics.Intensity,newcolnames)

# Compute Log2 intensities
global.Proteomics.Intensity_log <- log2(as.matrix(global.Proteomics.Intensity))

################################# PCA before any preprocecing ################################# 
cols <- c('#FF0000', '#FFFF00', '#00FF00', '#00FFFF', '#0000FF', '#FF00FF')
condition_cols <- c('#0A0200', '#00D1D1','#8AFFFF', 
                    '#0000D1','#8A8AFF',
                    '#D100D1','#FF8AFF',
                    '#00D100', '#8AFF8A')
# For some reason it does not work!!!!
# Starts here
p1 = plotQC(global.Proteomics.Intensity_log, 
            labels = "", 
            panel = "pca", 
            grps = global_experimental_design$dataset
)+
  labs(color = "tables")+
  scale_color_manual(values=cols)+
  ggplot2::ggtitle("PCA - Before any processing")
p1
# Check the quality of data
aaa<-colSums(is.na(global.Proteomics.Intensity_log))
names(which(colSums(is.na(global.Proteomics.Intensity_log))>2000))
# ends here

# alternatively use sparsepca
#install.packages("sparsepca")
library("sparsepca")

# Compute SPCA
spca_res <- spca(global.Proteomics.Intensity_log, k=5, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
#print(spca_res)
#summary(spca_res)
spercentVar <- spca_res$sdev^2/sum(spca_res$sdev^2)
sPCAvalues <- data.frame(spca_res$transform)
sPCAvalues$condition <- global_experimental_design$condition
sPCAvalues$dataset <- global_experimental_design$dataset
sPCAvalues$trial <- global_experimental_design$trial

p1<-ggplot(sPCAvalues, aes(x = X1, y = X2, colour = dataset)) +
  geom_point(size = 2) + 
  scale_color_manual(values=cols)+
  #theme(plot.title = element_text(size=12)) +
  theme_bw() +
  labs(color = "tables")+
  ggplot2::labs(
    x = paste0("PC1: ", round(spercentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(spercentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("PCA before processing (per table)")
p1

p2<-ggplot(sPCAvalues, aes(x = X1, y = X2, colour = condition)) +
  geom_point(size = 2) + 
  theme_bw() +
  scale_color_manual(values=condition_cols)+
  #theme(plot.title = element_text(size=10)) +
  #geom_text(label=rownames(PCAvalues), check_overlap = T) + # value labels
  ggplot2::labs(
    x = paste0("PC1: ", round(spercentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(spercentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("PCA before processing (per condition)")
p2

library(ggpubr)
plot <- ggpubr::ggarrange(p2, p1, nrow = 1)
annotate_figure(plot, top = text_grob("PCA before processing", face = "bold", size = 14))

################################# Boxplots before preprocecing ################################# 
# Plot the consolidated dataset
par(cex.main = 1)
par(cex.lab=0.5)
par(cex.axis=0.4)
boxplot(global.Proteomics.Intensity_log,main="All Tables",cex = 0.4,las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)

# Plot S1, S3 and S4 tables
par(mfrow=c(2,3))
par(cex.axis=0.5)
boxplot(global.Proteomics.Intensity_log[,c(1:20)], main="S1 Table", las=2, cex=0.8)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)
par(cex.axis=0.5)
boxplot(global.Proteomics.Intensity_log[,c(41:52)],main="S3 Table",las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)
boxplot(global.Proteomics.Intensity_log[,c(53:64)],main="S4 Table",las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)

# Plot S2, S5 and S6 tables
par(cex.axis=0.4)
boxplot(global.Proteomics.Intensity_log[,c(21:40)],main="S2 Table",las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)
par(cex.axis=0.4)
boxplot(global.Proteomics.Intensity_log[,c(65:76)],main="S5 Table",las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)
boxplot(global.Proteomics.Intensity_log[,c(77:88)],main="S6 Table",las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)
par(mfrow=c(1,1))

## Processing through PhosR ####
################################################################################
## Proteomics with quantification for at least 50% of the replicates in at least one of the conditions are retained
ppe_filtered <- selectGrps(global.Proteomics.Intensity_log, global_experimental_design$condition, 0.5, n=1) 
dim(ppe_filtered)

### Insted of median scaling: quantile normalisation
#ppe_filtered <- PhosR::medianScaling(ppe_filtered, scale = FALSE, grps = global_conditions$condition, reorder = FALSE, assay = NULL)

################################# Quantile normalization ################################# 
library(preprocessCore)
ppe_filtered <- normalize.quantiles(ppe_filtered,copy=TRUE)
colnames(ppe_filtered) <- newcolnames

################################# PCA after Quantile normalization ################################# 
cols <- c('#FF0000', '#FFFF00', '#00FF00', '#00FFFF', '#0000FF', '#FF00FF')
p2 = plotQC(ppe_filtered, 
            labels = "", 
            panel = "pca", 
            grps = global_experimental_design$dataset)+
  labs(color = "tables")+
  scale_color_manual(values=cols)+
  ggplot2::ggtitle("PCA - After Quantile normalization")
p2


## Imputation of Proteomics ####
################################################################################
set.seed(123)
ppe_imputed_tmp <- scImpute(ppe_filtered, 0.5, global_experimental_design$condition)[,colnames(ppe_filtered)]
ppe_imputed_tmp <- tImpute(ppe_imputed_tmp)

## Quantification plots filtered vs imputed ####
################################################################################
plotQC(ppe_filtered, 
       labels=colnames(ppe_filtered), 
       panel = "quantify", grps = global_experimental_design$dataset)+
  labs(color = "tables")+
  scale_color_manual(values=cols)+
  theme(axis.title.x=element_blank()) +  #remove x axis labels
  ggplot2::ggtitle("Dataset normalized and filtered (before imputation)")

plotQC(ppe_imputed_tmp, labels=colnames(ppe_imputed_tmp), 
       panel = "quantify", grps = global_experimental_design$dataset)+
  labs(color = "tables")+
  scale_color_manual(values=cols)+
  theme(axis.title.x=element_blank()) +  #remove x axis labels
  ggplot2::ggtitle("Dataset after Imputation")

## Dendogram plot filtered vs imputed ####
################################################################################

plotQC(ppe_imputed_tmp, labels=colnames(ppe_imputed_tmp), 
       panel = "dendrogram", grps = global_experimental_design$dataset)+
  labs(color = "tables")+
  scale_color_manual(values=cols)+
  #theme_gray(base_size = 5) +
  ggplot2::ggtitle("Dataset after Imputation")

## Batch effect correction using SVA ####
## Note: tried PhosR batch correction - didn't work well

modcombat <- model.matrix(~1, data=global_experimental_design)
trial <- global_experimental_design$trial
ppe_imputed_tmp <- ComBat(dat=ppe_imputed_tmp, batch=trial, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

## PCA inspection following imputation ####
################################################################################

#ppe_filtered ppe_imputed_tmp
pca <- prcomp(t(ppe_imputed_tmp), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAvalues <- data.frame(pca$x)
#PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- global_experimental_design$condition
#PCAvalues$stimulation <- experimental_design$stimulation
#PCAvalues$time <- experimental_design$time
PCAvalues$dataset <- global_experimental_design$dataset
PCAvalues$trial <- global_experimental_design$trial

condition_cols <- c('#0A0200', '#00D1D1','#8AFFFF', 
                    '#0000D1','#8A8AFF',
                    '#D100D1','#FF8AFF',
                    '#00D100', '#8AFF8A')
p2 <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 1.5) + 
  #scale_color_viridis(discrete = T) +
  theme_bw() +
  scale_color_manual(values=condition_cols)+
  #geom_text(label=rownames(PCAvalues), check_overlap = T) + # value labels
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("PCA after batch effect adjustment(per condition)")
p2

p3 <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = dataset)) +
  geom_point(size = 1.5) + 
  #scale_color_viridis(discrete = T) +
  scale_color_manual(values=cols)+
  theme_bw() +
  labs(color = "tables")+
  #geom_text(label=rownames(PCAvalues), check_overlap = T) + # value labels
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("PCA after batch effect adjustment (per table)")
p3



## Saving output ####
################################################################################
saveRDS(ppe_imputed_tmp, paste("artifacts_phos_imputed.Rds"))
saveRDS(experimental_design, paste("artifacts_experimental_design.Rds"))

