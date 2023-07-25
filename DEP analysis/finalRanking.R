
library(tibble)
#set directory
working_dir <- "C:\\Users\\mtheo\\Desktop\\RAF substrates 2023\\"
working_dir <- "D:\\Documents\\AMaria\\Jimmy\\2nd Trimester\\Project\\RAF substrates 2023\\"

results_folder <- paste0(working_dir,"Results\\")

# Find the downregulated genes
down_regulated_genes <- function(dt, tfit, contrast, results) {
  de <- which(dt[,contrast]==-1)
  length(de)
  de.genes <- tfit$genes[de]
  dim(de.genes)
  
  df <- results %>%
    filter(ID %in% de.genes)
  dim(df)
  # find dubplicated
  #df_md_1[duplicated(df_md_1$ID),]
  df_filtered <- df[c("ID","logFC")]
  return (df_filtered)
}

# Find the upregulated genes
up_regulated_genes <- function(dt, tfit, contrast, results) {
  de <- which(dt[,contrast]==1)
  length(de)
  de.genes <- tfit$genes[de]
  dim(de.genes)
  
  df <- results %>%
    filter(ID %in% de.genes)
  dim(df)
  # find dubplicated
  #df_md_1[duplicated(df_md_1$ID),]
  df_filtered <- df[c("ID","logFC")]
  return (df_filtered)
}

# Find the up/down regulated genes 
dis_regulated_genes <- function(dt, tfit, contrast, results) {
  de <- which(dt[,contrast]!=0)
  length(de)
  
  de.genes <- tfit$genes[de]
  dim(de.genes)
  
  df <- results %>%
    filter(ID %in% de.genes)
  dim(df)
  # find dubplicated
  #df_md_1[duplicated(df_md_1$ID),]
  df_filtered <- df[c("ID","logFC")]
  return (df_filtered)
}
# Find up/down regulated genes in Monomer_Dimer
# RAF1_BRAF
md1 <- dis_regulated_genes(dt_md, tfit_md, 1, BRAF.DMSOvsDim_results )
dim(md1)
colnames(md1) <- c("Gene.names", "BRAF.DMSOvsDim.LFC")
md2 <- down_regulated_genes(dt_md, tfit_md, 2, BRAF.DMSOvsSor_results )
dim(md2)
colnames(md2) <- c("Gene.names", "BRAF.DMSOvsSor.LFC")
md3 <- dis_regulated_genes(dt_md, tfit_md, 3, BRAF.DMSOvsSor.Dim_results )
dim(md3)
colnames(md3) <- c("Gene.names", "BRAF.DMSOvsSor.Dim.LFC")
md4 <- down_regulated_genes(dt_md, tfit_md, 4, BRAF.DMSOvsVem_results )
dim(md4)
colnames(md4) <- c("Gene.names", "BRAF.DMSOvsVem.LFC")
md5 <- dis_regulated_genes(dt_md, tfit_md, 5, BRAF.DMSOvsVem.Dim_results )
dim(md5)
colnames(md5) <- c("Gene.names", "BRAF.DMSOvsVem.Dim.LFC")

# RAF1_BRAFV600
md6 <- dis_regulated_genes(dt_md, tfit_md, 6, BRAFV600E.DMSOvsDim_results )
dim(md6)
colnames(md6) <- c("Gene.names", "BRAFV600E.DMSOvsDim.LFC")
md7 <- down_regulated_genes(dt_md, tfit_md, 7, BRAFV600E.DMSOvsSor_results )
dim(md7)
colnames(md7) <- c("Gene.names", "BRAFV600E.DMSOvsSor.LFC")
md8 <- dis_regulated_genes(dt_md, tfit_md, 8, BRAFV600E.DMSOvsSor.Dim_results )
dim(md8)
colnames(md8) <- c("Gene.names", "BRAFV600E.DMSOvsSor.Dim.LFC")
md9 <- down_regulated_genes(dt_md, tfit_md, 9, BRAFV600E.DMSOvsVem_results )
dim(md9)
colnames(md9) <- c("Gene.names", "BRAFV600E.DMSOvsVem.LFC")
md10 <- dis_regulated_genes(dt_md, tfit_md, 10, BRAFV600E.DMSOvsVem.Dim_results )
dim(md10)
colnames(md10) <- c("Gene.names", "BRAFV600E.DMSOvsVem.Dim.LFC")

#### join ####
# monomer_dimer BRAF
run.seq <- function(x) as.numeric(ave(paste(x), x, FUN = seq_along))
L <- list(md1, md2, md3, md4, md5)
L2 <- lapply(L, function(x) cbind(x, run.seq = run.seq(x$Gene.names)))
monomer_dimer_braf <- Reduce(function(...) merge(..., all = TRUE), L2)[-2]
dim(monomer_dimer_braf)
# Monomer_dimer BRAFV600
L <- list(md6, md7, md8, md9, md10)
L2 <- lapply(L, function(x) cbind(x, run.seq = run.seq(x$Gene.names)))
monomer_dimer_brafv600 <- Reduce(function(...) merge(..., all = TRUE), L2)[-2]
dim(monomer_dimer_brafv600)


# Find up/down regulated genes in Abiotin
abiotin1 <- dis_regulated_genes(dt_ab, tfit_ab, 1, FSBA.EF1.RAF1_results )
dim(abiotin1)
colnames(abiotin1) <- c("Gene.names", "FSBA.EF1.RAF1.LFC")
abiotin2 <- dis_regulated_genes(dt_ab, tfit_ab, 2, FSBA.RAF1.kinase_results )
dim(abiotin2)
colnames(abiotin2) <- c("Gene.names", "FSBA.RAF1.kinase.LFC")

#### join ####
L <- list(abiotin1, abiotin2)
L2 <- lapply(L, function(x) cbind(x, run.seq = run.seq(x$Gene.names)))
abiotin <- Reduce(function(...) merge(..., all = TRUE), L2)[-2]
dim(abiotin)

#####################
#### join monomer_dimer & abiotin ####
L <- list(monomer_dimer_braf, monomer_dimer_brafv600, abiotin)
L2 <- lapply(L, function(x) cbind(x, run.seq = run.seq(x$Gene.names)))
md_ab <- Reduce(function(...) merge(..., all = TRUE), L2)[-2]
dim(md_ab)
summary(md_ab)

# rank function by frequency of occurrence in different datasets & by the biggest logFC
# score each row: #of_occurence*Sum(|logFC|)
md_ab$Num_of_occurence <- rowSums(!is.na(md_ab[,2:13]))
md_ab$Total_logFC <- rowSums(abs(md_ab[,2:13]),na.rm=TRUE)
md_ab <- add_column(md_ab, score = md_ab$Num_of_occurence * md_ab$Total_logFC, .after = 1)
#Rank of “score” column in descending order, i.e. max score gets small rank
#in case of same values, average rank of both the values are assigned.
md_ab <- add_column(md_ab, rank = rank(desc(md_ab$score)), .after = 1)
ordered_md_ab <- md_ab[order(md_ab$rank), ]

write_xlsx(ordered_md_ab, paste0(results_folder, "Monomer_Dimer_vs_Abiotin_ranking.xlsx"))

# Find up/down regulated genes in KCLIP
kclip1 <- dis_regulated_genes(dt, tfit, 1, BRAFWT_results )
dim(kclip1)
colnames(kclip1) <- c("Gene.names", "BRAFWT.LFC")
kclip2 <- dis_regulated_genes(dt, tfit, 2, BRAFTK483M_results )
dim(kclip2)
colnames(kclip2) <- c("Gene.names", "BRAFTK483M.LFC")
kclip3 <- dis_regulated_genes(dt, tfit, 3, BRAFV600E_results )
dim(kclip3)
colnames(kclip3) <- c("Gene.names", "BRAFV600E.LFC")
kclip4 <- dis_regulated_genes(dt, tfit, 4, BRAFV600EK483M_results )
dim(kclip4)
colnames(kclip4) <- c("Gene.names", "BRAFV600EK483M.LFC")

#### join ####
L <- list(kclip1, kclip2, kclip3, kclip4)
L2 <- lapply(L, function(x) cbind(x, run.seq = run.seq(x$Gene.names)))
kclip <- Reduce(function(...) merge(..., all = TRUE), L2)[-2]
dim(kclip)

#####################
#### join monomer_dimer & kclip ####
L <- list(monomer_dimer_braf, monomer_dimer_brafv600, kclip)
L2 <- lapply(L, function(x) cbind(x, run.seq = run.seq(x$Gene.names)))
md_kc <- Reduce(function(...) merge(..., all = TRUE), L2)[-2]
dim(md_kc)
summary(md_kc)

# rank function by frequency of occurrence in different datasets & by the biggest logFC
# score each row: #of_occurence*Sum(|logFC|)
md_kc$Num_of_occurence <- rowSums(!is.na(md_kc[,2:13]))
md_kc$Total_logFC <- rowSums(abs(md_kc[,2:13]),na.rm=TRUE)
md_kc <- add_column(md_kc, score = md_kc$Num_of_occurence * md_kc$Total_logFC, .after = 1)
#Rank of “score” column in descending order, i.e. max score gets small rank
#in case of same values, average rank of both the values are assigned.
md_kc <- add_column(md_kc, rank = rank(desc(md_kc$score)), .after = 1)
ordered_md_kc <- md_kc[order(md_kc$rank), ]

write_xlsx(ordered_md_kc,paste0(results_folder,"Monomer_Dimer_vs_KCLIP_Ranking.xlsx"))
