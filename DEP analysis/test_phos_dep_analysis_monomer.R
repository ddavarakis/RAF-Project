
# Monomer - Dimer Dataset Analysis

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
dataset_dir <- paste0(working_dir, "BRAF monomer-dimer screens")
setwd(dataset_dir)
results_folder <- paste0(working_dir,"Results\\")
figures_folder <- paste0(working_dir,"Figures\\monomer\\")


all_interactors_filename = "All interactors.xlsx"

dset <- read_excel(all_interactors_filename, sheet ="proteinGroups")
names(dset)<-str_replace_all(names(dset), c(" " = "." ,"_" = "." , "," = "" ))
## remove potential contaminants ####
dset <- dset[is.na(dset$Reverse)  & is.na(dset$Contaminant), ]

#remove NA gene names
which(is.na(dset$Gene.names) )
dset <- dset[!is.na(dset$Gene.names), ]
# ----------- Monomer - Dimer dset ---------------
dset.Intensity <- dset[,grepl("Intensity *", colnames(dset))]
dset.Intensity[,1] <- NULL
colnames(dset.Intensity)

## Make Na those with zero intensity
dset.Intensity[dset.Intensity == 0] <- NA
#rownames(dset.Intensity) <- dset$Gene.names
originalcolnames <- colnames(dset.Intensity)

colnames(dset.Intensity) <- gsub("Intensity.24.", "RAF1.BRAFV600E.Vem.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.23.", "RAF1.BRAFV600E.Vem.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.22.", "RAF1.BRAFV600E.Sor.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.21.", "RAF1.BRAFV600E.Sor.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.20.", "RAF1.BRAFV600E.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.19.", "RAF1.BRAFV600E.DMSO.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.18.", "RAF1.BRAF.Vem.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.17.", "RAF1.BRAF.Vem.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.16.", "RAF1.BRAF.Sor.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.15.", "RAF1.BRAF.Sor.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.14.", "RAF1.BRAF.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.13.", "RAF1.BRAF.DMSO.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.12.", "RAF1.BRAFV600E.Vem.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.11.", "RAF1.BRAFV600E.Vem.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.10.", "RAF1.BRAFV600E.Sor.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.9.", "RAF1.BRAFV600E.Sor.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.8.", "RAF1.BRAFV600E.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.7.", "RAF1.BRAFV600E.DMSO.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.6.", "RAF1.BRAF.Vem.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.5.", "RAF1.BRAF.Vem.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.4.", "RAF1.BRAF.Sor.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.3.", "RAF1.BRAF.Sor.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.2.", "RAF1.BRAF.Dim.", colnames(dset.Intensity))
colnames(dset.Intensity) <- gsub("Intensity.1.", "RAF1.BRAF.DMSO.", colnames(dset.Intensity))

newcolnames <- colnames(dset.Intensity)
newcolnames

# Assume that: *_α: Control inside the experiment, *_b: 0m *_c: 5m, *_d: 15m, *_e: 30m, *_f:120m
timepoint <- as.factor(rep(c("Cntl","0m", "5m", "15m", "30m", "120m"), c(24)))
timepoint

# create group (=experiment) information
# Προσοχη: Βλεπω με ποια σειρα ειναι τα experiments στο αρχειο (δες samplenames)
# Και αναλογα ονοματισε τα experiments/groups
# ???? CHECK THEM
# 
group <- as.factor(rep(c("RAF1.BRAF.DMSO","RAF1.BRAFV600E.Sor.Dim", "RAF1.BRAFV600E.Vem",
                         "RAF1.BRAFV600E.Vem.Dim","RAF1.BRAF.DMSO", "RAF1.BRAF.Dim",
                         "RAF1.BRAF.Sor", "RAF1.BRAF.Sor.Dim", "RAF1.BRAF.Vem",
                         "RAF1.BRAF.Vem.Dim", "RAF1.BRAFV600E.DMSO", "RAF1.BRAF.Dim",
                         "RAF1.BRAFV600E.Dim", "RAF1.BRAFV600E.Sor", "RAF1.BRAFV600E.Sor.Dim",
                         "RAF1.BRAFV600E.Vem", "RAF1.BRAFV600E.Vem.Dim", "RAF1.BRAF.Sor", 
                         "RAF1.BRAF.Sor.Dim", "RAF1.BRAF.Vem", "RAF1.BRAF.Vem.Dim", 
                         "RAF1.BRAFV600E.DMSO", "RAF1.BRAFV600E.Dim", "RAF1.BRAFV600E.Sor"
), 
c(6,6,6,
  6,6,6,
  6,6,6, 
  6,6,6,
  6,6,6,
  6,6,6,
  6,6,6, 
  6,6,6
)
))

# Assume that Run1: Intensities from 1-12 and Run2: Intensities from 13-24
run <- as.factor(rep(c("Run1","Run1", "Run1",
                       "Run1","Run2", "Run2",
                       "Run2", "Run2", "Run2",
                       "Run2", "Run2", "Run1",
                       "Run2", "Run2", "Run2",
                       "Run2", "Run2", "Run1", 
                       "Run1", "Run1", "Run1", 
                       "Run1", "Run1", "Run1"
), 
c(6,6,6,
  6,6,6,
  6,6,6, 
  6,6,6,
  6,6,6,
  6,6,6,
  6,6,6, 
  6,6,6
)
))
condition <- newcolnames
ID <- c(1:144)

experimental_design_md <- data.frame(ID, condition, timepoint, group, run)
experimental_design_md


# Compute Log2 intensities
dset.Intensity_log <- log2(as.matrix(dset.Intensity))
# Insert gene names as rownames
rownames(dset.Intensity_log)<-dset$Gene.names

# get 12 distinct colors for group coloring
n <-12
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
group_colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
run_colors <- c('#FF0000', '#0000FF')

# Compute SPCA
spca_res <- spca(dset.Intensity_log, k=5, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
#print(spca_res)
#summary(spca_res)
spercentVar <- spca_res$sdev^2/sum(spca_res$sdev^2)
sPCAvalues <- data.frame(spca_res$transform)
sPCAvalues$condition <- experimental_design_md$condition
sPCAvalues$treatment <- experimental_design_md$treatment
sPCAvalues$group <- experimental_design_md$group
sPCAvalues$run <- experimental_design_md$run

p<-ggplot(sPCAvalues, aes(x = X1, y = X2, colour = group)) +
  geom_point(size = 2) + 
  scale_color_manual(values=group_colors)+
  #theme(plot.title = element_text(size=12)) +
  theme_bw() +
  labs(color = "Replicates")+
  ggplot2::labs(
    x = paste0("PC1: ", round(spercentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(spercentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("Monomer-dimer before processing (per group)")
p
# save plot to disk
filename = paste0(figures_folder,"monomer-pca-before-processing-per-group.png")
png(file=filename,
    width=600, height=350)
p
dev.off()

p<-ggplot(sPCAvalues, aes(x = X1, y = X2, colour = run)) +
  geom_point(size = 2) + 
  scale_color_manual(values=run_colors)+
  #theme(plot.title = element_text(size=12)) +
  theme_bw() +
  labs(color = "Replicates")+
  ggplot2::labs(
    x = paste0("PC1: ", round(spercentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(spercentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("Monomer-dimer before processing (per replicate)")
p
# save plot to disk
filename = paste0(figures_folder,"monomer-pca-before-processing-per-replicate.png")
png(file=filename,
    width=600, height=350)
p
dev.off()

################################# Boxplots before preprocecing ################################# 
# Plot the consolidated dataset
par(cex.main = 1)
par(cex.lab=0.5)
par(cex.axis=0.4)
boxplot(dset.Intensity_log,main="Monomer-dimer dataset - before processing",cex = 0.4,las=2)
mtext(side = 2, text = expression(bold("log(Intensity)")),line=1.5,cex=0.7)


## Processing through PhosR ####
################################################################################
## Proteomics with quantification for at least 50% of the replicates in at least one of the conditions are retained
ppe_filtered <- selectGrps(dset.Intensity_log, colnames(dset.Intensity_log), 0.5, n=1) 
dim(ppe_filtered)

### Instead of median scaling: quantile normalisation
#ppe_filtered <- PhosR::medianScaling(ppe_filtered, scale = FALSE, grps = global_conditions$condition, reorder = FALSE, assay = NULL)

################################# Quantile normalization ################################# 

ppe_filtered_norm <- normalize.quantiles(ppe_filtered,copy=TRUE)
colnames(ppe_filtered_norm) <- newcolnames
rownames(ppe_filtered_norm) <- rownames(ppe_filtered)

################################# PCA after Quantile normalization ################################# 
p = plotQC(ppe_filtered_norm, 
            labels = "", 
            panel = "pca", 
            grps = experimental_design_md$run)+
  labs(color = "tables")+
  scale_color_manual(values=run_colors)+
  ggplot2::ggtitle("PCA - After Quantile normalization (per replicate)")
p
# save plot to disk
filename = paste0(figures_folder,"monomer-pca-after-normalization-per-replicate.png")
png(file=filename,
    width=600, height=350)
p
dev.off()
p = plotQC(ppe_filtered_norm, 
           labels = "", 
           panel = "pca", 
           grps = experimental_design_md$group)+
  labs(color = "tables")+
  scale_color_manual(values=group_colors)+
  ggplot2::ggtitle("PCA - After Quantile normalization (per group)")
p
# save plot to disk
filename = paste0(figures_folder,"monomer-pca-after-normalization-per-group.png")
png(file=filename,
    width=600, height=350)
p
dev.off()


## Imputation of Proteomics ####
################################################################################
set.seed(123)
ppe_imputed_tmp_md <- scImpute(ppe_filtered_norm, 0.5, colnames(dset.Intensity_log))[,colnames(ppe_filtered_norm)]
ppe_imputed_tmp_md <- tImpute(ppe_imputed_tmp_md)

## Quantification plots filtered vs imputed ####
################################################################################
plotQC(ppe_filtered_norm, 
       labels=colnames(ppe_filtered_norm), 
       panel = "quantify", grps = experimental_design_md$group)+
  labs(color = "tables")+
  scale_color_manual(values=group_colors)+
  theme(axis.title.x=element_blank()) +  #remove x axis labels
  ggplot2::ggtitle("Monomer-dimer - normalized and filtered (before imputation) per group")

plotQC(ppe_imputed_tmp_md, labels=colnames(ppe_imputed_tmp_md), 
       panel = "quantify", grps = experimental_design_md$group)+
  labs(color = "tables")+
  scale_color_manual(values=group_colors)+
  theme(axis.title.x=element_blank()) +  #remove x axis labels
  ggplot2::ggtitle("Monomer-dimer - normalized and filtered (after imputation) per group")

## Dendogram plot filtered vs imputed ####
################################################################################

plotQC(ppe_imputed_tmp_md, labels=colnames(ppe_imputed_tmp_md), 
       panel = "dendrogram", grps = experimental_design_md$group)+
  labs(color = "tables")+
  scale_color_manual(values=group_colors)+
  #theme_gray(base_size = 5) +
  ggplot2::ggtitle("Monomer Dimer dendrogram - after imputation")

## Batch effect correction using SVA ####
## Note: tried PhosR batch correction - didn't work well

modcombat <- model.matrix(~1, data=experimental_design_md)
run <- experimental_design_md$run
ppe_imputed_tmp_md <- ComBat(dat=ppe_imputed_tmp_md, batch=run, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

## PCA inspection following imputation ####
################################################################################

#ppe_filtered ppe_imputed_tmp_md
pca <- prcomp(t(ppe_imputed_tmp_md), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAvalues <- data.frame(pca$x)
#PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- experimental_design_md$condition
PCAvalues$treatment <- experimental_design_md$treatment
PCAvalues$group <- experimental_design_md$group
PCAvalues$run <- experimental_design_md$run
p <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 1.5) + 
  #scale_color_viridis(discrete = T) +
  theme_bw() +
  scale_color_manual(values=group_colors)+
  #geom_text(label=rownames(PCAvalues), check_overlap = T) + # value labels
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("Monomer-dimer PCA after batch effect adjustment(per group)")
p
# save plot to disk
filename = paste0(figures_folder,"monomer-pca-after-batch-per-group.png")
png(file=filename,
    width=600, height=350)
p
dev.off()


p <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = run)) +
  geom_point(size = 1.5) + 
  #scale_color_viridis(discrete = T) +
  scale_color_manual(values=run_colors)+
  theme_bw() +
  labs(color = "tables")+
  #geom_text(label=rownames(PCAvalues), check_overlap = T) + # value labels
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))+
  ggplot2::ggtitle("Monomer-dimer PCA after batch effect adjustment(per replicate)")
p
# save plot to disk
filename = paste0(figures_folder,"monomer-pca-after-batch-per-group.png")
png(file=filename,
    width=600, height=350)
p
dev.off()

colnames(ppe_imputed_tmp_md)<-str_replace_all(colnames(ppe_imputed_tmp_md), c(" " = "." , "," = "" ))

## Saving output ####
saveRDS(ppe_imputed_tmp_md, paste0(results_folder,"Monomer_dimer_After_Preprocessing.Rds"))
saveRDS(experimental_design_md, paste0(results_folder,"Monomer_dimer_experimental_design_md.Rds"))

df <- as.data.frame(ppe_imputed_tmp_md)
df$Gene.names <- rownames(ppe_imputed_tmp_md)
write_xlsx(df,paste0(results_folder,"monomer_dimer_ppe_imputed_tmp.xlsx"))

backup_monomer_dimer <- ppe_imputed_tmp_md # keep it for backup

### Gene Expression Analysis with Limma ##########
# Refer to: 
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#organising-sample-information
# https://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/limma_edge_mdr.html
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# https://blog.devgenius.io/differential-gene-expression-analysis-using-limma-step-by-step-358da9d41c4e
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/design_mdmatrices.html
# https://hbctraining.github.io/dge_md_workshop_salmon_online/lessons/experimental_planning_considerations.html

#ppe_imputed_tmp_md <- backup

colnames(ppe_imputed_tmp_md) <- originalcolnames

dge_md <- DGEList(counts = ppe_imputed_tmp_md)
dim(dge_md)
class(dge_md)
colnames(dge_md)



# Organising sample information
# experimental design_md
# For simplicity, the "Intensity." substring from the column names of dge_mdList-object 
# if Intensity = 11, If LFQ = 15
samplenames <- substring(colnames(dge_md), 11, nchar(colnames(dge_md)))
samplenames
colnames(dge_md) <- samplenames

# Assume that: *_α: Control inside the experiment, *_b: 0m *_c: 5m, *_d: 15m, *_e: 30m, *_f:120m
timepoint <- as.factor(rep(c("Cntl","0m", "5m", "15m", "30m", "120m"), c(24)))
timepoint
dge_md$samples$timepoint <- timepoint

# create group (=experiment) information
# Προσοχη: Βλεπω με ποια σειρα ειναι τα experiments στο αρχειο (δες samplenames)
# Και αναλογα ονοματισε τα experiments/groups
# ???? CHECK THEM
# 
group <- as.factor(rep(c("RAF1.BRAF.DMSO","RAF1.BRAFV600E.Sor.Dim", "RAF1.BRAFV600E.Vem",
                         "RAF1.BRAFV600E.Vem.Dim","RAF1.BRAF.DMSO", "RAF1.BRAF.Dim",
                         "RAF1.BRAF.Sor", "RAF1.BRAF.Sor.Dim", "RAF1.BRAF.Vem",
                         "RAF1.BRAF.Vem.Dim", "RAF1.BRAFV600E.DMSO", "RAF1.BRAF.Dim",
                         "RAF1.BRAFV600E.Dim", "RAF1.BRAFV600E.Sor", "RAF1.BRAFV600E.Sor.Dim",
                         "RAF1.BRAFV600E.Vem", "RAF1.BRAFV600E.Vem.Dim", "RAF1.BRAF.Sor", 
                         "RAF1.BRAF.Sor.Dim", "RAF1.BRAF.Vem", "RAF1.BRAF.Vem.Dim", 
                         "RAF1.BRAFV600E.DMSO", "RAF1.BRAFV600E.Dim", "RAF1.BRAFV600E.Sor"
), 
c(6,6,6,
  6,6,6,
  6,6,6, 
  6,6,6,
  6,6,6,
  6,6,6,
  6,6,6, 
  6,6,6
)
))

dge_md$samples$group <- group

# Assume that Run1: Intensities from 1-12 and Run2: Intensities from 13-24
run <- as.factor(rep(c("Run1","Run1", "Run1",
                       "Run1","Run2", "Run2",
                       "Run2", "Run2", "Run2",
                       "Run2", "Run2", "Run1",
                       "Run2", "Run2", "Run2",
                       "Run2", "Run2", "Run1", 
                       "Run1", "Run1", "Run1", 
                       "Run1", "Run1", "Run1"
), 
c(6,6,6,
  6,6,6,
  6,6,6, 
  6,6,6,
  6,6,6,
  6,6,6,
  6,6,6, 
  6,6,6
)
))
run
dge_md$samples$run <- run
# check the experiment
dge_md$samples

# Organizing gene annotations
# In the dge_mdList object, associate count and sample information with gene annotations.
genes <- rownames(ppe_imputed_tmp_md)
dge_md$genes <- genes
dim(dge_md)


# Data Preprocessing
# Data pre-processing
# Transformations from the raw-scale
# FIG start
md_cpm <- cpm(dge_md)
md_lcpm <- cpm(dge_md, log=TRUE)

md_L <- mean(dge_md$samples$lib.size) * 1e-6
md_M <- median(dge_md$samples$lib.size) * 1e-6
c(md_L, md_M)
summary(md_lcpm)
# FIG end

# Removing genes that are lowly expressed
# 0% of genes in this dataset have zero counts across all  samples.
table(rowSums(dge_md$counts==0)==144)

# FIG start
# Produce Figure

md_lcpm.cutoff <- log2(10/md_M + 2/md_L)
nsamples <- ncol(dge_md)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(md_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,2.5), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=md_lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(md_lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
md_lcpm <- cpm(dge_md, log=TRUE)
plot(density(md_lcpm[,1]), col=col[1], lwd=2, ylim=c(0,2.5), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=md_lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(md_lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
par(mfrow=c(1,1))
# FIG end

# Normalising gene expression distributions
# density plot or boxplot  showing the per sample expression distributions, 
# is useful in determining whether any samples are dissimilar to others.
#

dge_md <- calcNormFactors(dge_md, method = "TMM")
dge_md$samples$norm.factors

# FIG start
# Plot DGE after normalization
md_lcpm <- cpm(dge_md, log=TRUE)
boxplot(md_lcpm, las=2, col=col, main="")
title(main="Monomer-Dimer - Normalized data (trimmed mean of M-values)",ylab="Log-cpm")
par(mfrow=c(1,1))
# FIG end


#Differential expression analysis
# Creating a design matrix and contrasts
# design
design_md <- model.matrix(~0+group)
design_md
colnames(design_md) <- gsub("group", "", colnames(design_md))
design_md
# contrasts
contr.matrix.md <- makeContrasts(
  BRAF.DMSOvsDim = RAF1.BRAF.DMSO - RAF1.BRAF.Dim, 
  BRAF.DMSOvsSor = RAF1.BRAF.DMSO - RAF1.BRAF.Sor, 
  BRAF.DMSOvsSor.Dim = RAF1.BRAF.DMSO - RAF1.BRAF.Sor.Dim,
  BRAF.DMSOvsVem = RAF1.BRAF.DMSO - RAF1.BRAF.Vem, 
  BRAF.DMSOvsVem.Dim = RAF1.BRAF.DMSO - RAF1.BRAF.Vem.Dim,
  BRAFV600E.DMSOvsDim = RAF1.BRAFV600E.DMSO - RAF1.BRAFV600E.Dim, 
  BRAFV600E.DMSOvsSor = RAF1.BRAFV600E.DMSO - RAF1.BRAFV600E.Sor, 
  BRAFV600E.DMSOvsSor.Dim = RAF1.BRAFV600E.DMSO - RAF1.BRAFV600E.Sor.Dim,
  BRAFV600E.DMSOvsVem = RAF1.BRAFV600E.DMSO - RAF1.BRAFV600E.Vem, 
  BRAFV600E.DMSOvsVem.Dim = RAF1.BRAFV600E.DMSO - RAF1.BRAFV600E.Vem.Dim, 
  levels = colnames(design_md))
contr.matrix.md

# Removing heteroscedascity from count data
# FIG START
par(mfrow=c(1,2))
v <- voom(dge_md, design_md, plot=TRUE)
v

vfit_md <- lmFit(v, design_md)
vfit_md <- contrasts.fit(vfit_md, contrasts=contr.matrix.md)
efit_md <- eBayes(vfit_md)
plotSA(efit_md, main="Final model: Mean-variance trend")
vfit_md
mtext("Monomer-dimer dataset", side = 3, line = -32, outer = TRUE,  font = 2)
par(mfrow=c(1,1))
# FIG end

# Fitting linear models for comparisons of interest

# Examining the number of DE genes
summary(decideTests(efit_md))

# stricter definition on significance
tfit_md <- treat(vfit_md, lfc=0.005)
dt_md <- decideTests(tfit_md)
summary(dt_md)
tfit_md

# Find the upregulated genes in 1st contrast BRAF.DMSOvsDim
de.md.up <- which(dt_md[,1]==1)
length(de.md.up)
tfit_md$genes[de.md.up]
# Find the downregulated genes in 1st contrast BRAF.DMSOvsDim
de.md.down <- which(dt_md[,1]==-1)
length(de.md.down)
tfit_md$genes[de.md.down]


# twenty of DE genes
# common between BRAF.DMSOvsDim and BRAF.DMSOvsSor.Dim
de.common.md <- which(dt_md[,1]!=0 & dt_md[,2]!=0)
length(de.common.md)
head(tfit_md$genes[de.common.md], n=20)
vennDiagram(dt_md[,1:2], circle.col=c("turquoise", "salmon"))

# Common genes across all BRAF conditions
# twenty of DE genes
de.common.md <- which(dt_md[,1]!=0 & dt_md[,2]!=0 & dt_md[,3]!=0 & dt_md[,4]!=0 & dt_md[,5]!=0)
length(de.common.md)
head(tfit_md$genes[de.common.md], n=20)
vennDiagram(dt_md[,1:5], circle.col=c("turquoise", "salmon", "blue", "magenta", "red"))


# Common genes across all BRAFV600E conditions
# twenty of DE genes
de.common.md <- which(dt_md[,6]!=0 & dt_md[,7]!=0 & dt_md[,8]!=0 & dt_md[,9]!=0 & dt_md[,10]!=0)
length(de.common.md)
head(tfit_md$genes[de.common.md], n=20)
vennDiagram(dt_md[,6:10], circle.col=c("turquoise", "salmon", "blue", "magenta", "red"))


#Write out the analysis results
dge_md_results = topTreat(tfit_md, n = Inf)
## Saving output ####
saveRDS(dge_md_results, paste0(results_folder,"Monomer_Dimer_DGE_Analysis_Results.Rds"))
write_xlsx(dge_md_results,paste0(results_folder,"Monomer_Dimer_dge_results.xlsx"))
#write.fit(tfit_md, dt_md, file=paste0(results_folder,"Monomer_Dime_tfit_results.txt"))

# Examining individual DE genes from top to bottom
filename <- paste0(results_folder,"Monomer_Dimer_dge_results.xlsx")
BRAF.DMSOvsDim_results <- topTreat(tfit_md, coef=1, n=Inf)
head(BRAF.DMSOvsDim_results)
write.xlsx(BRAF.DMSOvsDim_results,filename,sheetName="BRAF.DMSOvsDim",append = TRUE)

BRAF.DMSOvsSor_results  <- topTreat(tfit_md, coef=2, n=Inf)
head(BRAF.DMSOvsSor_results )
write.xlsx(BRAF.DMSOvsSor_results,filename,sheetName="BRAF.DMSOvsSor",append = TRUE)

BRAF.DMSOvsSor.Dim_results  <- topTreat(tfit_md, coef=3, n=Inf)
head(BRAF.DMSOvsSor.Dim_results )
write.xlsx(BRAF.DMSOvsSor.Dim_results,filename,sheetName="BRAF.DMSOvsSor.Dim",append = TRUE)

BRAF.DMSOvsVem_results  <- topTreat(tfit_md, coef=4, n=Inf)
head(BRAF.DMSOvsVem_results )
write.xlsx(BRAF.DMSOvsVem_results,filename,sheetName="BRAF.DMSOvsVem",append = TRUE)

BRAF.DMSOvsVem.Dim_results  <- topTreat(tfit_md, coef=5, n=Inf)
head(BRAF.DMSOvsVem.Dim_results )
write.xlsx(BRAF.DMSOvsVem.Dim_results,filename,sheetName="BRAF.DMSOvsVem.Dim",append = TRUE)

BRAFV600E.DMSOvsDim_results  <- topTreat(tfit_md, coef=6, n=Inf)
head(BRAFV600E.DMSOvsDim_results )
write.xlsx(BRAFV600E.DMSOvsDim_results,filename,sheetName="BRAFV600E.DMSOvsDim",append = TRUE)

BRAFV600E.DMSOvsSor_results  <- topTreat(tfit_md, coef=7, n=Inf)
head(BRAFV600E.DMSOvsSor_results )
write.xlsx(BRAFV600E.DMSOvsSor_results,filename,sheetName="BRAFV600E.DMSOvsSor",append = TRUE)

BRAFV600E.DMSOvsSor.Dim_results  <- topTreat(tfit_md, coef=8, n=Inf)
head(BRAFV600E.DMSOvsSor.Dim_results )
write.xlsx(BRAFV600E.DMSOvsSor.Dim_results,filename,sheetName="BRAFV600E.Sor.Dim",append = TRUE)

BRAFV600E.DMSOvsVem_results  <- topTreat(tfit_md, coef=9, n=Inf)
head(BRAFV600E.DMSOvsVem_results )
write.xlsx(BRAFV600E.DMSOvsVem_results,filename,sheetName="BRAFV600E.DMSOvsVem",append = TRUE)

BRAFV600E.DMSOvsVem.Dim_results  <- topTreat(tfit_md, coef=10, n=Inf)
head(BRAFV600E.DMSOvsVem.Dim_results )
write.xlsx(BRAFV600E.DMSOvsVem.Dim_results,filename,sheetName="BRAFV600E.DMSOvsVem.Dim",append = TRUE)


# Useful graphical representations of differential expression results
plotMD(tfit_md, column=1, status=dt_md[,1], main=colnames(tfit_md)[1], 
       xlim=c(8.6,9.6))
plotMD(tfit_md, column=1, status=dt_md[,2], main=colnames(tfit_md)[2], 
       xlim=c(8.6,9.6))
# Take snapshot from glMDPlot !!!
glMDPlot(tfit_md, coef=1, status=dt_md, main=colnames(tfit_md)[1],
         side.main="ENTREZID", counts=md_lcpm, groups=group, launch=TRUE)
glMDPlot(tfit_md, coef=2, status=dt_md, main=colnames(tfit_md)[2],
         side.main="ENTREZID", counts=md_lcpm, groups=group, launch=TRUE)
# FIG end


test.topgenes <- BRAFV600E.DMSOvsVem.Dim_results$ID[1:100]
i <- which(v$genes %in% test.topgenes)
