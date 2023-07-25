
# Abiotin Dataset Analysis

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
  library(xlsx)
})
#set directory
working_dir <- "C:\\Users\\mtheo\\Desktop\\RAF substrates 2023\\"
working_dir <- "D:\\Documents\\AMaria\\Jimmy\\2nd Trimester\\Project\\RAF substrates 2023\\"
dataset_dir <- paste0(working_dir, "RAF1_ATPbiotin_Gavin")
setwd(dataset_dir)
results_folder <- paste0(working_dir,"Results\\")
figures_folder <- paste0(working_dir,"Figures\\abiotin\\")

# Loading xlsx file

filename <- "RAF1 kinase assay_protein groups.xlsx"

dataset <- read_excel(filename)
dim(dataset)
names(dataset)<-str_replace_all(names(dataset), c(" " = "." ,"_" = ".", "-" = ".", "," = "" ))
colnames(dataset)

## remove potential contaminants ####
dataset <- dataset[is.na(dataset$Reverse)  & is.na(dataset$Potential.contaminant), ]

#remove NA gene names
which(is.na(dataset$Gene.names) )
dataset <- dataset[!is.na(dataset$Gene.names), ]
dim(dataset)

## formatting intensity ####

dataset.Intensity <- dataset[c(
                               #"Intensity.DMSO.No.RAF1.03",
                               #"Intensity.DMSO.EF1.RAF1.05",
                               #"Intensity.DMSO.RAF1.kinase.01",
                               "Gene.names",
                               "Intensity.FSBA.No.RAF1.04",
                               "Intensity.FSBA.EF1.RAF1.06",
                               "FSBA.RAF1.kinase.02")]

colnames(dataset.Intensity) <- c(
                                 #"Intensity.DMSO.No.RAF1.03",
                                 #"Intensity.DMSO.EF1.RAF1.05",
                                 #"Intensity.DMSO.RAF1.kinase.01",
                                 "Gene.names",
                                 "Intensity.FSBA.No.RAF1.04",
                                 "Intensity.FSBA.EF1.RAF1.06",
                                 "Intensity.FSBA.RAF1.kinase.02")
colnames(dataset.Intensity)
dim(dataset.Intensity)

dataset.Intensity <- tibble::column_to_rownames(dataset.Intensity, "Gene.names")

# Delete all rows having 0s in all columns
dataset.Intensity = dataset.Intensity[rowSums(dataset.Intensity[])>0,]
dim(dataset.Intensity)
# All values of Intensity.FSBA.No.RAF1.04 are 0 and the ComBat() function does not run !!!!
# Therefore replace 0 with 0.000000001
dataset.Intensity$Intensity.FSBA.No.RAF1.04[dataset.Intensity$Intensity.FSBA.No.RAF1.04 == 0] <- 0.000000001

# !!! dummy columns inserted only for running scImpute function !!!
dataset.Intensity$Intensity.FSBA.No.RAF1.04.dummy <- dataset.Intensity$Intensity.FSBA.No.RAF1.04
dataset.Intensity$Intensity.FSBA.EF1.RAF1.06.dummy <- dataset.Intensity$Intensity.FSBA.EF1.RAF1.06
dataset.Intensity$Intensity.FSBA.RAF1.kinase.02.dummy <- dataset.Intensity$Intensity.FSBA.RAF1.kinase.02

originalcolnames <- colnames(dataset.Intensity)

## Setting up experimental design_ab ####
experimental_design_ab <- data.frame(c(1,2,3,4,5,6),
                         c(
                           #"Intensity.DMSO.No.RAF1.03",
                           #"Intensity.DMSO.EF1.RAF1.05",
                           #"Intensity.DMSO.RAF1.kinase.01",
                           "Intensity.FSBA.No.RAF1.04",
                           "Intensity.FSBA.EF1.RAF1.06",
                           "Intensity.FSBA.RAF1.kinase.02",
                           "Intensity.FSBA.No.RAF1.04.dummy",
                           "Intensity.FSBA.EF1.RAF1.06.dummy",
                           "Intensity.FSBA.RAF1.kinase.02.dummy"
                           ),
                         c("run1", "run1", "run1","run2", "run2", "run2")
                         )

colnames(experimental_design_ab) <- c("ID","condition", "run")


## Make Na those with zero intensity
dataset.Intensity[dataset.Intensity == 0] <- NA

colnames(dataset.Intensity) <- gsub("Intensity.FSBA.RAF1.kinase.02.dummy", "Intensity.FSBA.RAF1.kinase.02", colnames(dataset.Intensity))
colnames(dataset.Intensity) <- gsub("Intensity.FSBA.EF1.RAF1.06.dummy", "Intensity.FSBA.EF1.RAF1.06", colnames(dataset.Intensity))
colnames(dataset.Intensity) <- gsub("Intensity.FSBA.No.RAF1.04.dummy", "Intensity.FSBA.No.RAF1.04", colnames(dataset.Intensity))


# Compute Log2 intensities
dataset.Intensity_log <- log2(as.matrix(dataset.Intensity))

n <-6
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cond_colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
run_colors <- c('#FF0000', '#0000FF')

# Compute SPCA
spca_res <- spca(dataset.Intensity_log, k=5, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
#print(spca_res)
#summary(spca_res)
spercentVar <- spca_res$sdev^2/sum(spca_res$sdev^2)
sPCAvalues <- data.frame(spca_res$transform)
sPCAvalues$condition <- experimental_design_ab$condition
sPCAvalues$run <- experimental_design_ab$run

p<-ggplot(sPCAvalues, aes(x = X1, y = X2, colour = condition)) +
  geom_point(size = 2) + 
  scale_color_manual(values=cond_colors)+
  #theme(plot.title = element_text(size=12)) +
  theme_bw() +
  labs(color = "Condition")+
  ggplot2::labs(
    x = paste0("PC1: ", round(spercentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(spercentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("Abiotin - before processing (per condition)")
p
# save plot to disk
filename = paste0(figures_folder,"abiotin-pca-before-processing-per-group.png")
png(file=filename,
    width=600, height=350)
p
dev.off()

p<-ggplot(sPCAvalues, aes(x = X1, y = X2, colour = run)) +
  geom_point(size = 2) + 
  scale_color_manual(values=run_colors)+
  #theme(plot.title = element_text(size=12)) +
  theme_bw() +
  labs(color = "Run")+
  ggplot2::labs(
    x = paste0("PC1: ", round(spercentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(spercentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("Abiotin before processing (per replicate)")
p
# save plot to disk
filename = paste0(figures_folder,"abiotin-pca-before-processing-per-replicate.png")
png(file=filename,
    width=600, height=350)
p
dev.off()

################################# Boxplots before preprocecing ################################# 
# Plot the consolidated dataset
par(cex.main = 1)
par(cex.lab=0.5)
par(cex.axis=0.4)
boxplot(dataset.Intensity_log,main="Abiotin dataset - before processing",cex = 0.4,las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)


################################# Quantile normalization ################################# 

ppe_filtered_norm <- normalize.quantiles(dataset.Intensity_log,copy=TRUE)
colnames(ppe_filtered_norm) <- colnames(dataset.Intensity)
rownames(ppe_filtered_norm) <- rownames(dataset.Intensity_log)

################################# PCA after Quantile normalization ################################# 
p = plotQC(ppe_filtered_norm, 
           labels = "", 
           panel = "pca", 
           grps = experimental_design_ab$run)+
  labs(color = "Runs")+
  scale_color_manual(values=run_colors)+
  ggplot2::ggtitle("Abiotin - PCA - After Quantile normalization (per replicate)")
p
# save plot to disk
filename = paste0(figures_folder,"abiotin-pca-after-normalization-per-replicate.png")
png(file=filename,
    width=600, height=350)
p
dev.off()
p = plotQC(ppe_filtered_norm, 
           labels = "", 
           panel = "pca", 
           grps = experimental_design_ab$condition)+
  labs(color = "tables")+
  scale_color_manual(values=cond_colors)+
  ggplot2::ggtitle("Abiotin - PCA - After Quantile normalization (per group)")
p
# save plot to disk
filename = paste0(figures_folder,"abiotin-pca-after-normalization-per-group.png")
png(file=filename,
    width=600, height=350)
p
dev.off()


## Imputation  ####
################################################################################

set.seed(123)
ppe_imputed_tmp_ab <- scImpute(ppe_filtered_norm, 0.5, colnames(ppe_filtered_norm))[,colnames(ppe_filtered_norm)]
ppe_imputed_tmp_ab <- tImpute(ppe_imputed_tmp_ab)

## Quantification plots filtered vs imputed ####
################################################################################
plotQC(ppe_filtered_norm, 
       labels=colnames(ppe_filtered_norm), 
       panel = "quantify", grps = experimental_design_ab$condition)+
  labs(color = "condition")+
  scale_color_manual(values=cond_colors)+
  theme(axis.title.x=element_blank()) +  #remove x axis labels
  ggplot2::ggtitle("Abiotin - normalized and filtered (before imputation)")

plotQC(ppe_imputed_tmp_ab, labels=colnames(ppe_imputed_tmp_ab), 
       panel = "quantify", grps = experimental_design_ab$condition)+
  labs(color = "tables")+
  scale_color_manual(values=cond_colors)+
  theme(axis.title.x=element_blank()) +  #remove x axis labels
  ggplot2::ggtitle("Abiotin -  after Imputation")

## Dendogram plot filtered vs imputed ####
################################################################################

plotQC(ppe_imputed_tmp_ab, labels=colnames(ppe_imputed_tmp_ab), 
       panel = "dendrogram", grps = experimental_design_ab$condition)+
  labs(color = "tables")+
  scale_color_manual(values=cond_colors)+
  #theme_gray(base_size = 5) +
  ggplot2::ggtitle("Abiotin - after Imputation")

## Batch effect correction ####

modcombat <- model.matrix(~1, data=experimental_design_ab)
run <- experimental_design_ab$run
ppe_imputed_tmp_ab <- ComBat(dat=ppe_imputed_tmp_ab, batch=run, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)



## PCA inspection following imputation ####
################################################################################

#ppe_filtered ppe_imputed_tmp_ab
pca <- prcomp(t(ppe_imputed_tmp_ab), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAvalues <- data.frame(pca$x)
#PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- experimental_design_ab$condition
PCAvalues$run <- experimental_design_ab$run

p <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 1.5) + 
  #scale_color_viridis(discrete = T) +
  theme_bw() +
  #labs(color = "Condition")+
  scale_color_manual(values=cond_colors)+
  #geom_text(label=rownames(PCAvalues), check_overlap = T) + # value labels
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("Abiotin - PCA after batch effect adjustment(per condition)")
p
# save plot to disk
filename = paste0(figures_folder,"abiotin-pca-after-normalization-per-condition.png")
png(file=filename,
    width=600, height=350)
p
dev.off()


p <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = run)) +
  geom_point(size = 1.5) + 
  #scale_color_viridis(discrete = T) +
  scale_color_manual(values=run_colors)+
  theme_bw() +
  labs(color = "Run")+
  #geom_text(label=rownames(PCAvalues), check_overlap = T) + # value labels
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("Abiotin PCA after batch effect adjustment (per run)")
p
# save plot to disk
filename = paste0(figures_folder,"abiotin-pca-after-normalization-per-run.png")
png(file=filename,
    width=600, height=350)
p
dev.off()

colnames(ppe_imputed_tmp_ab)<-str_replace_all(colnames(ppe_imputed_tmp_ab), c(" " = "." , "," = "" ))

## Saving output ####
saveRDS(ppe_imputed_tmp_ab, paste0(results_folder,"Abiotin_After_Preprocessing.Rds"))
saveRDS(experimental_design_ab, paste0(results_folder,"Abiotin_experimental_design.Rds"))

df <- as.data.frame(ppe_imputed_tmp_ab)
df$Gene.names <- rownames(ppe_imputed_tmp_ab)
write_xlsx(df,paste0(results_folder,"abiotin_ppe_imputed_tmp.xlsx"))

### Gene Expression Analysis with Limma ##########
# Refer to: 
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#organising-sample-information
# https://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/limma_edge_abr.html
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# https://blog.devgenius.io/differential-gene-expression-analysis-using-limma-step-by-step-358da9d41c4e
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/design_abmatrices.html
# https://hbctraining.github.io/dge_ab_workshop_salmon_online/lessons/experimental_planning_considerations.html

backup_abiotin <- ppe_imputed_tmp_ab # keep it for backup
colnames(ppe_imputed_tmp_ab) <- originalcolnames

dge_ab <- DGEList(counts = ppe_imputed_tmp_ab)
dim(dge_ab)
class(dge_ab)
colnames(dge_ab)

# Organising sample information
# experimental design_ab
# For simplicity, the "Intensity." substring from the column names of dge_abList-object 
# if Intensity = 11, If LFQ = 15
samplenames <- substring(colnames(dge_ab), 11, nchar(colnames(dge_ab)))
samplenames
colnames(dge_ab) <- samplenames


group <- rep(c("FSBA.No.RAF1","FSBA.EF1.RAF1","FSBA.RAF1.kinase"), 2)
group <- as.factor(group)
group
dge_ab$samples$group <- group
# organize samples per run (2nd run is dummy)
dge_ab$samples$run <- run

# Organizing gene annotations
# In the dge_abList object, associate count and sample information with gene annotations.
genes <- rownames(ppe_imputed_tmp_ab)
dge_ab$genes <- genes
dim(dge_ab)

# Data Preprocessing
# Data pre-processing
# Transformations from the raw-scale
# FIG start
ab_cpm <- cpm(dge_ab)
ab_lcpm <- cpm(dge_ab, log=TRUE)

ab_L <- mean(dge_ab$samples$lib.size) * 1e-6
ab_M <- median(dge_ab$samples$lib.size) * 1e-6
c(ab_L, ab_M)
summary(ab_lcpm)
# FIG end

# Removing genes that are lowly expressed
# 0% of genes in this dataset have zero counts across all  samples.
table(rowSums(dge_ab$counts==0)==6)

# FIG start
# Produce Figure
ab_lcpm.cutoff <- log2(10/ab_M + 2/ab_L)
nsamples <- ncol(dge_ab)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(ab_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,14), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=ab_lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(ab_lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
ab_lcpm <- cpm(dge_ab, log=TRUE)
plot(density(ab_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,14), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=ab_lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(ab_lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
par(mfrow=c(1,1))
# FIG end

# Normalising gene expression distributions
# density plot or boxplot  showing the per sample expression distributions, 
# is useful in determining whether any samples are dissimilar to others.
#
dge_ab <- calcNormFactors(dge_ab, method = "TMM")
dge_ab$samples$norm.factors

# FIG start
# Plot DGE after normalization
ab_lcpm <- cpm(dge_ab, log=TRUE)
boxplot(ab_lcpm, las=2, col=col, main="")
title(main="Abiotin - Normalized data (trimmed mean of M-values)",ylab="Log-cpm")
par(mfrow=c(1,1))
# FIG end


#Differential expression analysis
# Creating a design_ab matrix and contrasts
# design_ab
design_ab <- model.matrix(~0+group)
design_ab
colnames(design_ab) <- gsub("group", "", colnames(design_ab))
design_ab
# contrasts
contr.matrix.ab <- makeContrasts(
  FSBA.EF1.RAF1 = FSBA.No.RAF1 - FSBA.EF1.RAF1, 
  FSBA.RAF1.kinase = FSBA.No.RAF1 - FSBA.RAF1.kinase, 
  levels = colnames(design_ab))
contr.matrix.ab

# Removing heteroscedascity from count data
# FIG start
par(mfrow=c(1,2))
v <- voom(dge_ab, design_ab, plot=TRUE)
v

vfit_ab <- lmFit(v, design_ab)
vfit_ab <- contrasts.fit(vfit_ab, contrasts=contr.matrix.ab)
efit_ab <- eBayes(vfit_ab)
plotSA(efit_ab, main="Final model: Mean-variance trend")
vfit_ab
mtext("Abiotin dataset", side = 3, line = -25, outer = TRUE,  font = 2)
par(mfrow=c(1,1))
# FIG end

# Fitting linear models for comparisons of interest

# Examining the number of DE genes
summary(decideTests(efit_ab))
# -----

# stricter definition on significance
tfit_ab <- treat(vfit_ab, lfc=0.5)
dt_ab <- decideTests(tfit_ab)
summary(dt_ab)
tfit_ab
# Regulated genes (up and down)
de <- which(dt_ab[,2]!=0)


# Find the upregulated genes
de.up.ab <- which(dt_ab[,1]==1)
length(de.up.ab)
tfit_ab$genes[de.up.ab]
# Find the downregulated genes
de.down.ab <- which(dt_ab[,1]==-1)
length(de.down.ab)
tfit_ab$genes[de.down.ab]

par(mfrow=c(1,1))
# twenty of DE genes
# common between FSBA.EF1.RAF1 and FSBA.RAF1.kinase
de.common.ab <- which(dt_ab[,1]!=0 & dt_ab[,2]!=0)
length(de.common.ab)
head(tfit_ab$genes[de.common.ab], n=20)
vennDiagram(dt_ab[,1:2], circle.col=c("turquoise", "salmon"))
# Top 10 genes for FSBA.RAF1.kinase
m <- topTable(tfit_ab,coef=2,conf=0.95,sort="none")
dim(m)
m

#Write out the analysis results
dge_ab_results = topTreat(tfit_ab, n = Inf)
## Saving output ####
saveRDS(dge_ab_results,  paste0(results_folder,"Abiotin_DGE_Analysis_Results.Rds"))
write_xlsx(dge_ab_results,paste0(results_folder,"Abiotin_dge_results.xlsx"))
write.fit(tfit_ab, dt_ab, file=paste0(results_folder,"Abiotin_tfit_results.txt"))

#test_df = read.table('Abiotin_tfit_ab_results.txt',sep='\t', header=TRUE)

# Examining individual DE genes from top to bottom
filename <- paste0(results_folder,"Abiotin_dge_results.xlsx")
FSBA.EF1.RAF1_results <- topTreat(tfit_ab, coef=1, n=Inf)
head(FSBA.EF1.RAF1_results)
write.xlsx(FSBA.EF1.RAF1_results,filename,sheetName="FSBA.EF1.RAF1",append = TRUE)

FSBA.RAF1.kinase_results  <- topTreat(tfit_ab, coef=2, n=Inf)
head(FSBA.RAF1.kinase_results )
write.xlsx(FSBA.RAF1.kinase_results,filename,sheetName="FSBA.RAF1.kinase",append = TRUE)

# FIG start
# Useful graphical representations of differential expression results
plotMD(tfit_ab, column=1, status=dt_ab[,1], main=colnames(tfit_ab)[1], 
       xlim=c(9.7,12))
plotMD(tfit_ab, column=2, status=dt_ab[,2], main=colnames(tfit_ab)[2], 
       xlim=c(9.7,12))
# Take snapshot from glMDPlot !!!
glMDPlot(tfit_ab, coef=1, status=dt_ab, main=colnames(tfit_ab)[1],
         side.main="ENTREZID", counts=ab_lcpm, groups=group, launch=TRUE, html = "MD-Plot-1")
glMDPlot(tfit_ab, coef=2, status=dt_ab, main=colnames(tfit_ab)[2],
         side.main="ENTREZID", counts=ab_lcpm, groups=group, launch=TRUE, html = "MD-Plot-2")
# FIG end

# FIG start
fsba.ef1.raf1.topgenes <- FSBA.EF1.RAF1_results$ID[1:100]
i <- which(v$genes %in% fsba.ef1.raf1.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(ab_lcpm[i,], scale="row",
          labRow=v$genes[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
# FIG end
