#BiocManager::install("Go.db")
#BiocManager::install("limma")

suppressPackageStartupMessages({
  library("readxl")
  library("writexl")
  library(tidyverse)
  library(dplyr)
  library("GO.db")
  library("limma")
  library("org.Hs.eg.db") # for annotations
})
#set directory
setwd("C:\\Users\\mtheo\\Desktop\\RAF substrates 2023")
 
# -------------- Quality Check of RAF1_ATPbiotin_Gavin ----------------

abiotin_filename = "RAF1_ATPbiotin_Gavin\\RAF1 kinase assay_protein groups.xlsx"
abiotin_data <- read_excel(abiotin_filename, sheet = "proteinGroups RAF kinase assay")
names(abiotin_data)<-str_replace_all(names(abiotin_data), c(" " = "." , "," = "", "_" = ".", "-" = "."   ))


abiotin_intensity <- abiotin_data[c("Gene.names",
                                         "FSBA.RAF1.kinase.02",
                                         "Intensity.DMSO.EF1.RAF1.05",
                                         "Intensity.DMSO.No.RAF1.03",
                                         "Intensity.DMSO.RAF1.kinase.01",
                                         "Intensity.FSBA.EF1.RAF1.06",
                                         "Intensity.FSBA.No.RAF1.04")]

colnames(abiotin_intensity) <- c("Gene.names",
                                 "Intensity.FSBA.RAF1.kinase",
                                 "Intensity.DMSO.EF1.RAF1",
                                 "Intensity.DMSO.No.RAF1",
                                 "Intensity.DMSO.RAF1.kinase",
                                 "Intensity.FSBA.EF1.RAF1",
                                 "Intensity.FSBA.No.RAF1"
                                 )
# Drop rows with Gene names empty
no_of_nulls <- sum(is.na(abiotin_intensity$Gene.names))
abiotin_intensity<-abiotin_intensity[!(is.na(abiotin_intensity$Gene.names)),]

# FSBA columns!!!
# Drop them rows with Zero in Intensity_FSBA_EF1_RAF1 and Intensity_FSBA_RAF1_kinase columns
fsba_intensity<-abiotin_intensity[!(abiotin_intensity$Intensity.FSBA.EF1.RAF1==0 &
             abiotin_intensity$Intensity.FSBA.RAF1.kinase==0),]

# Get those rows that have only RAF1 counts (thus EF1-RAF1 = 0)
fsba_only_RAF1 <- fsba_intensity[(fsba_intensity['Intensity.FSBA.EF1.RAF1'] == 0),]
# Get those rows that have only EF1-RAF1 counts (thus RAF1 = 0)
fsba_only_EF1_RAF1 <- fsba_intensity[(fsba_intensity['Intensity.FSBA.RAF1.kinase'] == 0),]
# Get those rows that have both RAF1 counts and EF1-RAF1 counts
fsba_both <- fsba_intensity[(fsba_intensity['Intensity.FSBA.RAF1.kinase'] != 0 &
                             fsba_intensity['Intensity.FSBA.EF1.RAF1'] != 0),]

# ----------- nice to be computed -----------------
# Separate rows with multiple genes
fsba_only_RAF1_expanded <- fsba_only_RAF1 %>%
        separate_longer_delim(c(Gene.names), delim = ";")
fsba_only_EF1_RAF1_expanded <- fsba_only_EF1_RAF1 %>%
  separate_longer_delim(c(Gene.names), delim = ";")

# convert gene symbols to entrez ids
symbols_only_RAF1 <- fsba_only_RAF1_expanded$Gene.names
Genes_only_RAF1 <- mapIds(org.Hs.eg.db, symbols_only_RAF1, 'ENTREZID', 'SYMBOL')
symbols_only_EF1_RAF1 <- fsba_only_EF1_RAF1_expanded$Gene.names
Genes_only_EF1_RAF1 <- mapIds(org.Hs.eg.db, symbols_only_EF1_RAF1, 'ENTREZID', 'SYMBOL')

# get GO pathways
g_only_RAF1 <- goana(Genes_only_RAF1, species="Hs")
g_only_RAF1_list <- topGO(g_only_RAF1)

g_only_EF1_RAF1 <- goana(Genes_only_EF1_RAF1, species="Hs")
g_only_EF1_RAF1_list <- topGO(g_only_EF1_RAF1)

common <- intersect(g_only_RAF1_list$Term, g_only_EF1_RAF1_list$Term)

#keg <- kegga(Genes)
#topKEGG(keg, n=15, truncate=34)

# DMSO columns
dmso_intensity<-abiotin_intensity[!(abiotin_intensity$Intensity.DMSO.EF1.RAF1==0 &
                                      abiotin_intensity$Intensity.DMSO.RAF1.kinase==0),]

# Get those rows that have only RAF1 counts (thus EF1-RAF1 = 0)
dmso_only_RAF1 <- dmso_intensity[(dmso_intensity['Intensity.DMSO.EF1.RAF1'] == 0),]
# Get those rows that have only EF1-RAF1 counts (thus RAF1 = 0)
dmso_only_EF1_RAF1 <- dmso_intensity[(dmso_intensity['Intensity.DMSO.RAF1.kinase'] == 0),]
# Get those rows that have both RAF1 counts and EF1-RAF1 counts
dmso_both <- dmso_intensity[(dmso_intensity['Intensity.DMSO.RAF1.kinase'] != 0 &
                               dmso_intensity['Intensity.DMSO.EF1.RAF1'] != 0),]


# -------------- Quality Check of momomer dimer ----------------

### Load All interactors.xlsx
all_interactors_filename = "BRAF monomer-dimer screens\\All interactors.xlsx"
all_interactors_data <- read_excel(all_interactors_filename, sheet ="proteinGroups")
names(all_interactors_data)<-str_replace_all(names(all_interactors_data), c(" " = "." , "," = "", "_" = ".", "-" = "."   ))

# Common between All_Interactors and Abiotin
common <- intersect(fsba_intensity$Gene.names, all_interactors_data$Gene.names)


### Load "Differential interactors identified and quantified by log fold change.xlsx"
monomer_dimer_filename = "BRAF monomer-dimer screens\\Differential interactors identified and quantified by log fold change.xlsx"
jointed_monomer_dimer_data = data.frame()
sheet_names <- excel_sheets(monomer_dimer_filename) 
sheet_names <- sheet_names[sheet_names != "Legend"]; # Remove 'Legend' sheet
for(p in sheet_names){
  #p = sheet_names[1] # only for testing!!!
  monomer_dimer_data <- read_excel(monomer_dimer_filename, sheet = p)
  monomer_dimer_data <- monomer_dimer_data[c("Gene", "lfc", "FDR")]
  s <- sum(is.na(monomer_dimer_data$Gene))
  print(paste("in sheet: ",p," found: ",s," rows with null protein name"))
  monomer_dimer_data<-monomer_dimer_data[!(is.na(monomer_dimer_data$Gene)),]
  

  number_of_dublicates <- dim(monomer_dimer_data[duplicated(monomer_dimer_data$Gene),])[1] 
  dubs <- monomer_dimer_data[duplicated(monomer_dimer_data$Gene),]
  for(g in dubs$Gene){
    df <- monomer_dimer_data[monomer_dimer_data$Gene == g,]
    if(nrow(df)>2){
      print("????")
    }else{
      if(df$lfc[1]>0 && df$lfc[2]>0){
        m <- max(df$lfc)
        monomer_dimer_data <- monomer_dimer_data %>%  filter(!(Gene==g & lfc != m))
      }else if(df$lfc[1]<0 && df$lfc[2]<0){
        m <- min(df$lfc)
        monomer_dimer_data <- monomer_dimer_data %>%  filter(!(Gene==g & lfc != m))
      }else{
        print(paste(g," !!!!"))
        if(df$lfc[1]+df$lfc[2]>0){
          # get the positive 
          monomer_dimer_data <- monomer_dimer_data %>%  filter(!(Gene==g & lfc < 0))
        }else{
          # get the negative
          monomer_dimer_data <- monomer_dimer_data %>%  filter(!(Gene==g & lfc > 0))
        }
      }
      
    }
  }
  colnames(monomer_dimer_data)[2] = p #rename lfc to sheet_name
  colnames(monomer_dimer_data)[3] = paste(p,".FDR") #rename FDR to "sheet_name.FDR"
  names(monomer_dimer_data)<-str_replace_all(names(monomer_dimer_data), c('\\+' = "."))
  if(p == sheet_names[1]){
    jointed_monomer_dimer_data <- monomer_dimer_data
  }else{
    jointed_monomer_dimer_data <- merge(x = jointed_monomer_dimer_data, y = monomer_dimer_data, 
                                        by = "Gene",no.dups=TRUE)
  }
}

# Common between Monomer-Dimer and Abiotin
common <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data$Gene)


# BrafV600E DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO < 0, ])
which(jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO == 0)
common <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO > 0, ]$Gene)

# BrafV600E DMSO vs Sor
jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor <- jointed_monomer_dimer_data$RAF1.BRAFV600E.Sor - jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor < 0, ])
common_sor <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor < 0, ]$Gene)

# BrafV600E DMSO vs Vem
jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem <- jointed_monomer_dimer_data$RAF1.BRAFV600E.Vem - jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem < 0, ])
common_vem <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem > 0, ]$Gene)
common <- intersect(common_sor,common_vem)

# BrafV600E DMSO vs Dim
jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Dim <- jointed_monomer_dimer_data$RAF1.BRAFV600E.Dim - jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Dim > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Dim < 0, ])
common_Dim_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Dim > 0, ]$Gene)
common_Dim_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Dim < 0, ]$Gene)

# BrafV600E DMSO vs Sor_Dim
jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor_Dim <- jointed_monomer_dimer_data$RAF1.BRAFV600E.Sor.Dim - jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor_Dim > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor_Dim < 0, ])
common_Sor_Dim_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor_Dim > 0, ]$Gene)
common_Sor_Dim_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Sor_Dim < 0, ]$Gene)

# BrafV600E DMSO vs Vem_Dim
jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem_Dim <- jointed_monomer_dimer_data$RAF1.BRAFV600E.Vem.Dim - jointed_monomer_dimer_data$RAF1.BRAFV600E.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem_Dim > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem_Dim < 0, ])
common_Vem_Dim_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem_Dim > 0, ]$Gene)
common_Vem_Dim_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Brafv600e_DMSO_vs_Vem_Dim < 0, ]$Gene)

# Braf DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$RAF1.BRAF.DMSO > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$RAF1.BRAF.DMSO < 0, ])
which(jointed_monomer_dimer_data$RAF1.BRAF.DMSO == 0)
common_braf_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$RAF1.BRAF.DMSO > 0, ]$Gene)
common_braf_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$RAF1.BRAF.DMSO < 0, ]$Gene)

# Braf DMSO vs Sor
jointed_monomer_dimer_data$Braf_DMSO_vs_Sor <- jointed_monomer_dimer_data$RAF1.BRAF.Sor - jointed_monomer_dimer_data$RAF1.BRAF.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Sor > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Sor < 0, ])
common_braf_Sor_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Sor > 0, ]$Gene)
common_braf_Sor_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Sor < 0, ]$Gene)

# Braf DMSO vs Vem
jointed_monomer_dimer_data$Braf_DMSO_vs_Vem <- jointed_monomer_dimer_data$RAF1.BRAF.Vem - jointed_monomer_dimer_data$RAF1.BRAF.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Vem > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Vem < 0, ])
common_braf_Vem_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Vem > 0, ]$Gene)
common_braf_Vem_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Vem < 0, ]$Gene)

# Braf DMSO vs Dim
jointed_monomer_dimer_data$Braf_DMSO_vs_Dim <- jointed_monomer_dimer_data$RAF1.BRAF.Dim - jointed_monomer_dimer_data$RAF1.BRAF.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Dim > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Dim < 0, ])
common_braf_Dim_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Dim > 0, ]$Gene)
common_braf_Dim_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Dim < 0, ]$Gene)

# Braf DMSO vs Sor_Dim
jointed_monomer_dimer_data$Braf_DMSO_vs_Sor_Dim <- jointed_monomer_dimer_data$RAF1.BRAF.Sor.Dim - jointed_monomer_dimer_data$RAF1.BRAF.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Sor_Dim > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Sor_Dim < 0, ])
common_braf_Sor_Dim_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Sor_Dim > 0, ]$Gene)
common_braf_Sor_Dim_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Sor_Dim < 0, ]$Gene)

# Braf DMSO vs Vem_Dim
jointed_monomer_dimer_data$Braf_DMSO_vs_Vem_Dim <- jointed_monomer_dimer_data$RAF1.BRAF.Vem.Dim - jointed_monomer_dimer_data$RAF1.BRAF.DMSO
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Vem_Dim > 0, ])
nrow(jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Vem_Dim < 0, ])
common_braf_Vem_Dim_g <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Vem_Dim > 0, ]$Gene)
common_braf_Vem_Dim_l <- intersect(fsba_intensity$Gene.names, jointed_monomer_dimer_data[jointed_monomer_dimer_data$Braf_DMSO_vs_Vem_Dim < 0, ]$Gene)



# -------------TO BE DELETED ----------------------------------------

# Control = Intensity.FSBA.No.RAF1
# Check control's quality
which(abiotin_intensity$Intensity.FSBA.No.RAF1 >0)
# Only 2 proteins have Control > 0
df[ which(df$TR == 4), "RE"]
abiotin_intensity[which(abiotin_intensity$Intensity.FSBA.No.RAF1 >0), "Gene.names"]
abiotin_intensity[which(abiotin_intensity$Intensity.DMSO.No.RAF1 >0), "Gene.names"]



abiotin_intensity$logFC.FSBA.RAF1 <- log2(abiotin_intensity$Intensity.FSBA.RAF1.kinase/abiotin_intensity$Intensity.FSBA.No.RAF1)
abiotin_intensity$logFC.FSBA.EF1.RAF1 <- log2(abiotin_intensity$Intensity.FSBA.EF1.RAF1/abiotin_intensity$Intensity.FSBA.No.RAF1)

abiotin_intensity <- abiotin_intensity %>% 
  mutate(FSBA.RAF1.greaterthan.FSBA.No.RAF1 = if_else(Intensity.FSBA.RAF1.kinase == Intensity.FSBA.No.RAF1, TRUE, FALSE))

ggplot(abiotin_intensity, aes(Intensity.FSBA.RAF1.kinase, logFC.FSBA.RAF1, colour = FSBA.RAF1.greaterthan.FSBA.No.RAF1), size = 8) + geom_point()

