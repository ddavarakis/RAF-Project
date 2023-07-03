
# KCLIP Dataset Analysis

suppressPackageStartupMessages({
  library(PhosR)
  library(dplyr)
  library(tidyverse)
  library(tibble)
  library(sva)
  library("readxl")
  library("writexl")
  library("sparsepca")
  library(preprocessCore)
  
  library(limma)
  library(Glimma)
  library(edgeR)
  
  library(gplots)
  library(RColorBrewer)
})
#set directory
setwd("C:\\Users\\mtheo\\Desktop\\RAF substrates 2023\\BRAF KCLIP 260221")

# Loading xlsx files

file.list <- c("Table S1_BRAFWT_MaxQuantRawData.xlsx", "Table S2_BRAFV600E_MaxQuantRawData.xlsx")
#filename <- "Table S1_BRAFWT_MaxQuantRawData.xlsx"
for (filename in file.list) {
  dataset <- read_excel(filename, skip=4)
  names(dataset)<-str_replace_all(names(dataset), c(" " = "." , "," = "" ))
  #print(filename)
  if(filename == "Table S1_BRAFWT_MaxQuantRawData.xlsx"){
    conditions <- data.frame(c(1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20), 
                             c(rep("BRAFWT No ATPArN3",4),
                               rep("BRAFWT with ATPArN3",4),
                               rep("BRAFK483M No ATPArN3",4),
                               rep("BRAFK483M with ATPArN3",4),
                               rep("Bead control",4)),
                             rep(1:4, 5),
                             rep(1,20))
  }
  if(filename == "Table S2_BRAFV600E_MaxQuantRawData.xlsx"){
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

  colnames(conditions) <- c("ID","condition","trial","dataset")
  
  conditions$trial <- as.factor(conditions$trial)
  conditions$dataset <- as.factor(conditions$dataset)
  

  ## remove potential contaminants ####
  dataset <- dataset[is.na(dataset$Reverse)  & is.na(dataset$Potential.contaminant), ]

  #remove NA gene names
  which(is.na(dataset$Gene.names) )
  dataset <- dataset[!is.na(dataset$Gene.names), ]
  

  ## formatting intensity ####
  patterns <- c("BRAF*", "Bead*", "Intensity *")
  dataset.Intensity <- dataset[,grepl(paste(patterns, collapse="|"), colnames(dataset))]
  dataset.Intensity[,1] <- NULL
  

  
  ## Setting up experimental design for phosphoproteomics ####
  
  experimental_design <- data.frame("ID" = sub("_.*", "", colnames(dataset.Intensity)), row.names  = colnames(dataset.Intensity))
  experimental_design$ID <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
  experimental_design <- merge(experimental_design, conditions, by = "ID", all.x = T)


  ## Make Na those with zero intensity
  dataset.Intensity[dataset.Intensity == 0] <- NA

  
  if(filename == "Table S1_BRAFWT_MaxQuantRawData.xlsx"){
    S1_dataset.Intensity <- data.frame(dataset.Intensity)
    rownames(S1_dataset.Intensity) <- dataset$Gene.names
    S1_experimental_design <- data.frame(experimental_design)
    S1_conditions <- data.frame(conditions)
  }
  if(filename == "Table S2_BRAFV600E_MaxQuantRawData.xlsx"){
    #S2_Proteomics.Intensity_log <- data.frame(Proteomics.Intensity_log) 
    S2_dataset.Intensity <- data.frame(dataset.Intensity)
    rownames(S2_dataset.Intensity) <- dataset$Gene.names
    S2_experimental_design <- data.frame(experimental_design)
    S2_conditions <- data.frame(conditions)
  }

}# End For
# correct the conditions
S2_experimental_design$ID <- as.numeric(S2_experimental_design$ID) +20
global_experimental_design = rbind(S1_experimental_design,S2_experimental_design)

# correct the conditions
S2_conditions$ID <- as.numeric(S2_conditions$ID) +20
global_conditions = rbind(S1_conditions,S2_conditions)


# Make row names separate column
S1_dataset.Intensity <- tibble::rownames_to_column(S1_dataset.Intensity, "Gene_name")
S2_dataset.Intensity <- tibble::rownames_to_column(S2_dataset.Intensity, "Gene_name")


# Check whether there are proteins in S2 that does not exist in S1
# There are 337 proteins in S2 that does not exist in S1
diff <- S2_dataset.Intensity %>% filter(!Gene_name %in% S1_dataset.Intensity$Gene_name)

#### join S1 + S2 ####
joined_df <- merge(x = S1_dataset.Intensity, y = S2_dataset.Intensity, 
                   by = "Gene_name", all=TRUE)

colnames(joined_df) <- gsub('x','_S1',names(joined_df))
colnames(joined_df) <- gsub('.y','_S2',names(joined_df))

global.kclip.Intensity <- joined_df

# Drop Gene_name column
global.kclip.Intensity <- global.kclip.Intensity[,!names(global.kclip.Intensity) %in% c("Gene_name")]

# Rename column names
newcolnames <- as.data.frame(t(global_experimental_design['condition']))
global.kclip.Intensity <- setNames(global.kclip.Intensity,newcolnames)

# Compute Log2 intensities
global.kclip.Intensity_log <- log2(as.matrix(global.kclip.Intensity))
# Insert gene names as rownames
rownames(global.kclip.Intensity_log)<-joined_df[,1]

################################# PCA before any preprocecing ################################# 
cols <- c('#FF0000', '#FFFF00', '#00FF00', '#00FFFF', '#0000FF', '#FF00FF')
condition_cols <- c('#0A0200', '#00D1D1','#8AFFFF', 
                    '#0000D1','#8A8AFF',
                    '#D100D1','#FF8AFF',
                    '#00D100', '#8AFF8A')

# Compute SPCA
spca_res <- spca(global.kclip.Intensity_log, k=5, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
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
boxplot(global.kclip.Intensity_log,main="All Tables",cex = 0.4,las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)

# Plot S1, ad S2 tables
par(mfrow=c(2,1))
par(cex.axis=0.5)
boxplot(global.kclip.Intensity_log[,c(1:20)], main="S1 Table", las=2, cex=0.8)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)
par(cex.axis=0.4)
boxplot(global.kclip.Intensity_log[,c(21:40)],main="S2 Table",las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)

par(mfrow=c(1,1))

## Processing through PhosR ####
################################################################################
## Proteomics with quantification for at least 50% of the replicates in at least one of the conditions are retained
ppe_filtered <- selectGrps(global.kclip.Intensity_log, global_experimental_design$condition, 0.5, n=1) 
dim(ppe_filtered)

### Instead of median scaling: quantile normalisation
#ppe_filtered <- PhosR::medianScaling(ppe_filtered, scale = FALSE, grps = global_conditions$condition, reorder = FALSE, assay = NULL)

################################# Quantile normalization ################################# 

ppe_filtered_norm <- normalize.quantiles(ppe_filtered,copy=TRUE)
colnames(ppe_filtered_norm) <- newcolnames
rownames(ppe_filtered_norm) <- rownames(ppe_filtered)

################################# PCA after Quantile normalization ################################# 
cols <- c('#FF0000', '#FFFF00', '#00FF00', '#00FFFF', '#0000FF', '#FF00FF')
p2 = plotQC(ppe_filtered_norm, 
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
ppe_imputed_tmp <- scImpute(ppe_filtered_norm, 0.5, global_experimental_design$condition)[,colnames(ppe_filtered_norm)]
ppe_imputed_tmp <- tImpute(ppe_imputed_tmp)

## Quantification plots filtered vs imputed ####
################################################################################
plotQC(ppe_filtered_norm, 
       labels=colnames(ppe_filtered_norm), 
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

colnames(ppe_imputed_tmp)<-str_replace_all(colnames(ppe_imputed_tmp), c(" " = "." , "," = "" ))

## Saving output ####
saveRDS(ppe_imputed_tmp, paste("KCLIP_After_Imputation.Rds"))
saveRDS(experimental_design, paste("experimental_design.Rds"))

df <- as.data.frame(ppe_imputed_tmp)
df$Gene.names <- rownames(ppe_imputed_tmp)
write_xlsx(df,"ppe_imputed_tmp.xlsx")

backup <- ppe_imputed_tmp # keep it for backup

### Gene Expression Analysis with Limma ##########
# Refer to: 
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#organising-sample-information
# https://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/limma_edger.html
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# https://blog.devgenius.io/differential-gene-expression-analysis-using-limma-step-by-step-358da9d41c4e
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/experimental_planning_considerations.html

# KCLIP dataset: Control is "No ATPArN3" 
# drop Bead.control columns
ppe_imputed_tmp <- df[ , -which(names(df) %in% c("Bead.control"))]
# drop Gene.names column
ppe_imputed_tmp$Gene.names <- NULL

# create a DGEList-object 
dge <- DGEList(counts = as.matrix(ppe_imputed_tmp))
dim(dge)
class(dge)
colnames(dge)

# Organising sample information
# sample-level information related to the experimental design 
# needs to be associated with the columns of the counts matrix.
# organize samples as : "No ATPArN3" as Control and "With ATPArN3" as Treatment
group <- c(rep(c("BRAFWT.Cntl","BRAFWT.Trt","BRAFK483M.Cntl","BRAFK483M.Trt"), c(4)),
              rep(c("BRAFV600E.Cntl","BRAFV600E.Trt","BRAFV600EK483M.Cntl","BRAFV600EK483M.Trt"), c(4)))
group
group <- as.factor(group)
group
dge$samples$group <- group
# organize samples per trial
trial <- rep(c(rep(1, 4),rep(2, 4),rep(3, 4),rep(4, 4)),2)
dge$samples$trial <- trial
# organize samples per dataset (S1 or S2)
ds <- c(rep(1, 16),rep(2, 16))
dge$samples$dataset <- ds

# Organizing gene annotations
# In the DGEList object, associate count and sample information with gene annotations.
genes <- rownames(ppe_imputed_tmp)
dge$genes <- genes

# Data Preprocessing


# Removing genes that are lowly expressed
# 0% of genes in this dataset have zero counts across all  samples.
table(rowSums(dge$counts==0)==32)
# Run it in case of having genes with zero counts accoss all samples
keep.exprs <- filterByExpr(dge, group=group)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)


# Normalising gene expression distributions
# density plot or boxplot  showing the per sample expression distributions, 
# is useful in determining whether any samples are dissimilar to others.
#
dge <- calcNormFactors(dge, method = "TMM")
dge$samples$norm.factors

#Differential expression analysis
# Creating a design matrix and contrasts
# design
design <- model.matrix(~0+group+trial)
design
colnames(design) <- gsub("group", "", colnames(design))
design
# contrasts
contr.matrix <- makeContrasts(
  BRAFWT = BRAFWT.Cntl - BRAFWT.Trt,
  BRAFK483M = BRAFK483M.Cntl - BRAFK483M.Trt,
  BRAFV600E = BRAFV600E.Cntl - BRAFV600E.Trt,
  BRAFV600EK483M = BRAFV600EK483M.Cntl - BRAFV600EK483M.Trt,
  levels = colnames(design))
contr.matrix

# Removing heteroscedascity from count data
par(mfrow=c(1,2))
v <- voom(dge, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
vfit
par(mfrow=c(1,1))

# Fitting linear models for comparisons of interest

# Examining the number of DE genes
summary(decideTests(efit))

# stricter definition on significance
tfit <- treat(vfit, lfc = 0.005)
dt <- decideTests(tfit)
summary(dt)
tfit
# Common genes between BRAFWT and BRAFK483M
# twenty of DE genes
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

# Common genes accross all
# twenty of DE genes
de.common <- which(dt[,1]!=0 & dt[,2]!=0 & dt[,3]!=0 & dt[,4]!=0)
length(de.common)
head(tfit$genes[de.common], n=20)
vennDiagram(dt[,1:4], circle.col=c("turquoise", "salmon", "blue", "magenta"))

#Write out the analysis results
dge_results = topTreat(tfit, n = Inf)
## Saving output ####
saveRDS(dge_results, paste("KCLIP_DGE_Analysis_Results.Rds"))
df <- as.data.frame(ppe_imputed_tmp)
df$Gene.names <- rownames(ppe_imputed_tmp)
write_xlsx(dge_results,"KCLIP_dge_results.xlsx")
write.fit(tfit, dt, file="KCLIP_tfit_results.txt")

# Examining individual DE genes from top to bottom
BRAFWT_results <- topTreat(tfit, coef=1, n=Inf)
head(BRAFWT_results)

library(xlsx)
write.xlsx(BRAFWT_results,"KCLIP_dge_results.xlsx",sheetName="BRAFWT_DEP",append = TRUE)
BRAFTK483M_results  <- topTreat(tfit, coef=2, n=Inf)
head(BRAFTK483M_results )
write.xlsx(BRAFTK483M_results,"KCLIP_dge_results.xlsx",sheetName="BRAFTK483M_DEP",append = TRUE)
BRAFV600E_results  <- topTreat(tfit, coef=3, n=Inf)
head(BRAFV600E_results )
write.xlsx(BRAFV600E_results,"KCLIP_dge_results.xlsx",sheetName="BRAFV600E_DEP",append = TRUE)
BRAFV600EK483M_results  <- topTreat(tfit, coef=4, n=Inf)
head(BRAFV600EK483M_results )
write.xlsx(BRAFV600EK483M_results,"KCLIP_dge_results.xlsx",sheetName="BRAFV600EK483M_DEP",append = TRUE)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(6,13))
#glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         #side.main="ENTREZID", counts=lcpm, groups=group, launch=TRUE)


BRAFWT.topgenes <- BRAFWT_results$ID[1:100]
i <- which(v$genes %in% BRAFWT.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
#heatmap.2(lcpm[i,], scale="row",
          #labRow=v$genes[i], labCol=group, 
          #col=mycol, trace="none", density.info="none", 
          #margin=c(8,6), lhei=c(2,10), dendrogram="column")
