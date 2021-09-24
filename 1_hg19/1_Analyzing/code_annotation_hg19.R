## =============================================================================
## Code for comparing the annotation files which was downloaded from public databases
## See the GitHub Page for more information:
## https://github.com/theresetruong/0_DATA_REF/blob/main/2_annot_data_ref/Gene_Annotation_Strategy/hg19/1_Analyzing
## Author : Yazdan Asgari, Last Edited:14/06/2021
## https://cesp.inserm.fr/en/equipe/exposome-and-heredity
## =============================================================================

# PATH for the annotation data (extracted version)
#======================================================
setwd("PATH_OF_DOWNLOADED_ANNOTATION_DATA")

#======================================================
# Packages needed for running all parts of the code
#======================================================
# BiocManager::install("rtracklayer")
# BiocManager::install("ggplot")
# BiocManager::install("ggvenn")

#======================================================
# import annotation data from different databases
#======================================================
# Needed for importing gtf & gff file formats
library(rtracklayer)

gencode <- rtracklayer::import('gencode.v19.annotation.gtf')
gencode=as.data.frame(gencode)

gencode_superset <- rtracklayer::import('gencode.v19.chr_patch_hapl_scaff.annotation.gtf')
gencode_superset=as.data.frame(gencode_superset)

ENSEMBL <- rtracklayer::import('ENSEMBL_Homo_sapiens.GRCh37.87.gtf')
ENSEMBL=as.data.frame(ENSEMBL)

Biomart <- read.delim('Biomart_GRCh37_martquery_0503134431_814.txt',header = TRUE,sep = "\t", quote = "\"",
                      dec = ".", fill = TRUE)

Biomart_unique <- read.delim('Biomart_GRCh37_unique_results_only_martquery_0503134950_260.txt',header = TRUE,sep = "\t", quote = "\"",
                             dec = ".", fill = TRUE)

UCSC <- read.delim('ucsc_GRCh37',header = TRUE,sep = "\t", quote = "\"",
                   dec = ".", fill = TRUE)

REMAP <- rtracklayer::import('REMAP_gencode.v28lift37.annotation.gtf')
REMAP=as.data.frame(REMAP)

FANTOM <- read.delim('FANTOM_hg19.txt',header = TRUE,sep = "\t", quote = "\"",
                     dec = ".", fill = TRUE)

PE <- rtracklayer::import('PE_GRCh37_latest_genomic.gff')
PE=as.data.frame(PE)

NCBI <- rtracklayer::import('NCBI_GRCh37_latest_genomic.gff')
NCBI=as.data.frame(NCBI)

#======================================================
# Cheking some features for the imported data
#======================================================
length(Biomart$Gene.type)
u <- unique(Biomart$Gene.type)

table(unlist(Biomart$Gene.type))
table(unlist(ENSEMBL$gene_biotype))
table(unlist(gencode_superset$gene_type))
table(unlist(REMAP$gene_type))

#==================================================================================
# 1- Compare PE & NCBI (protein_coding Genes)
#==================================================================================
coding_PE <- subset(PE$gene, PE$gene_biotype=="protein_coding")
coding_NCBI <- subset(NCBI$gene, NCBI$gene_biotype=="protein_coding")

# coding_PE <- subset(PE$Name, PE$gene_biotype=="protein_coding")
# coding_NCBI <- subset(NCBI$Name, NCBI$gene_biotype=="protein_coding")
# coding_PE <- subset(PE$transcript_id, PE$gene_biotype=="protein_coding")
# coding_NCBI <- subset(NCBI$transcript_id, NCBI$gene_biotype=="protein_coding")

x <- list(
  coding_PE = coding_PE,
  coding_NCBI = coding_NCBI
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 3
)

U <- setdiff(coding_PE,coding_NCBI)
length(U)
U
V <- setdiff(coding_NCBI,coding_PE)
length(V)
V

write.table(U, file = "coding_PE.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)
write.table(V, file = "coding_NCBI.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)

#=========================================================
# 2- Compare gencode & gencode_superset (All & coding genes)
#=========================================================
x <- list(
  gencode = gencode$gene_name,
  gencode_superset = gencode_superset$gene_name
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 3
)

coding_gencode <- subset(gencode$gene_name, gencode$gene_type=="protein_coding")
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")

x <- list(
  coding_gencode = coding_gencode,
  coding_gencode_superset = coding_gencode_superset
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 3
)

A <- setdiff(coding_gencode_superset,coding_gencode)
length(A)

#=========================================================
# 3- Compare Biomart & Biomart_unique  (All & coding genes)
#=========================================================
x <- list(
  Biomart = Biomart$Gene.name,
  Biomart_unique = Biomart_unique$Gene.name
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 3
)

coding_Biomart <- subset(Biomart$Gene.name, Biomart$Gene.type=="protein_coding")
coding_Biomart_unique <- subset(Biomart_unique$Gene.name, Biomart_unique$Gene.type=="protein_coding")

x <- list(
  coding_Biomart = coding_Biomart,
  coding_Biomart_unique = coding_Biomart_unique
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 3
)

#============================================================================================
# 4- Compare gencode_superset & ENSEMBL & Biomart & REMAP_gencode_ebi_ftp (Protein_coding Genes)
#============================================================================================
coding_Biomart <- subset(Biomart$Gene.name, Biomart$Gene.type=="protein_coding")
coding_ENSEMBL <- subset(ENSEMBL$gene_name, ENSEMBL$gene_biotype=="protein_coding")
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_name, REMAP$gene_type=="protein_coding")

x <- list(
  coding_Biomart = coding_Biomart,
  coding_ENSEMBL = coding_ENSEMBL,
  coding_gencode_superset = coding_gencode_superset,
  coding_REMAP_gencode_ebi_ftp = coding_REMAP_gencode_ebi_ftp
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

A <- setdiff(coding_gencode_superset,union(union(coding_Biomart,coding_ENSEMBL),coding_REMAP_gencode_ebi_ftp))
length(A)

#=========================================================================================
# 5- Compare gencode_superset & Biomart & REMAP_gencode_ebi_ftp & PE (Protein_coding Genes)
#=========================================================================================
coding_Biomart <- subset(Biomart$Gene.name, Biomart$Gene.type=="protein_coding")
coding_PE <- subset(PE$gene, PE$gene_biotype=="protein_coding")
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_name, REMAP$gene_type=="protein_coding")

# library(ggplot2)
# library(ggvenn)
x <- list(
  coding_Biomart = coding_Biomart,
  coding_PE = coding_PE,
  coding_gencode_superset = coding_gencode_superset,
  coding_REMAP_gencode_ebi_ftp = coding_REMAP_gencode_ebi_ftp
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

U <- setdiff(UCSC$hg19.kgXref.geneSymbol,union(union(Biomart$Gene.name,gencode_superset$gene_name),REMAP$gene_name))
length(U)

#============================================================================================
# 6- Compare gencode_superset & Biomart & REMAP_gencode_ebi_ftp & NCBI (Protein_coding Genes)
#============================================================================================
coding_Biomart <- subset(Biomart$Gene.name, Biomart$Gene.type=="protein_coding")
coding_NCBI <- subset(NCBI$gene, NCBI$gene_biotype=="protein_coding")
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_name, REMAP$gene_type=="protein_coding")

x <- list(
  coding_Biomart = coding_Biomart,
  coding_NCBI = coding_NCBI,
  coding_gencode_superset = coding_gencode_superset,
  coding_REMAP_gencode_ebi_ftp = coding_REMAP_gencode_ebi_ftp
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

U <- setdiff(coding_FANTOM,union(union(coding_REMAP_gencode_ebi_ftp,coding_gencode_superset),coding_Biomart))
length(U)

#============================================================================================
# 7- Compare gencode_superset & Biomart & REMAP_gencode_ebi_ftp & FANTOM (Protein_coding Genes)
#============================================================================================
coding_Biomart <- subset(Biomart$Gene.name, Biomart$Gene.type=="protein_coding")
coding_FANTOM <- subset(FANTOM$gene_symbol, FANTOM$Gene_class=="coding_mRNA")
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_name, REMAP$gene_type=="protein_coding")

x <- list(
  coding_Biomart = coding_Biomart,
  coding_FANTOM = coding_FANTOM,
  coding_gencode_superset = coding_gencode_superset,
  coding_REMAP_gencode_ebi_ftp = coding_REMAP_gencode_ebi_ftp
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

U <- setdiff(coding_FANTOM,union(union(coding_REMAP_gencode_ebi_ftp,coding_gencode_superset),coding_Biomart))
length(U)

#============================================================================================
# 8- Consider union of 5 files (gencode_superset, Biomart, REMAP_gencode_ebi_ftp, PE, NCBI)
# (Protein_coding Genes)
#============================================================================================
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")
coding_Biomart <- subset(Biomart$Gene.name, Biomart$Gene.type=="protein_coding")
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_name, REMAP$gene_type=="protein_coding")
coding_PE <- subset(PE$gene, PE$gene_biotype=="protein_coding")
coding_NCBI <- subset(NCBI$gene, NCBI$gene_biotype=="protein_coding")

gene_annot_coding <- union(coding_gencode_superset,
                           union(coding_Biomart,
                                 union(coding_REMAP_gencode_ebi_ftp,
                                             union(coding_PE,coding_NCBI)
                                       )
                                 )
                           )
length(gene_annot_coding)

#============================================================================================
# 9- Consider union of 3 files (gencode_superset, Biomart, REMAP_gencode_ebi_ftp)
# (Protein_coding Genes) / based on Gene Symbol 
#============================================================================================
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")
coding_Biomart <- subset(Biomart$Gene.name, Biomart$Gene.type=="protein_coding")
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_name, REMAP$gene_type=="protein_coding")

gene_annot_coding <- union(coding_gencode_superset,
                           union(coding_Biomart,coding_REMAP_gencode_ebi_ftp)
                           )
length(gene_annot_coding)

#============================================================================================
# 10- Consider union of 3 files (gencode_superset, Biomart, REMAP_gencode_ebi_ftp)
# (Protein_coding Genes) / based on ENSEMBL ID
#============================================================================================
coding_gencode_superset <- subset(gencode_superset$gene_id, gencode_superset$gene_type=="protein_coding")
coding_Biomart <- subset(Biomart$Gene.stable.ID, Biomart$Gene.type=="protein_coding")
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_id, REMAP$gene_type=="protein_coding")

coding_gencode_superset <- substr(coding_gencode_superset,1,15)
coding_REMAP_gencode_ebi_ftp <- substr(coding_REMAP_gencode_ebi_ftp,1,15)

gene_annot_coding <- union(coding_gencode_superset,
                           union(coding_Biomart,coding_REMAP_gencode_ebi_ftp)
                           )

length(gene_annot_coding)

#============================================================================================
# 11- Consider union of 2 files (PE, NCBI)
# (Protein_coding Genes) / based on Gene Symbol 
#============================================================================================
coding_PE <- subset(PE$gene, PE$gene_biotype=="protein_coding")
coding_NCBI <- subset(NCBI$gene, NCBI$gene_biotype=="protein_coding")

gene_annot_coding <- union(coding_PE,coding_NCBI)
length(gene_annot_coding)

#============================================================================================
# 12- Consider union of 6 files (gencode_superset, Biomart, REMAP_gencode_ebi_ftp, PE, NCBI, FANTOM)
# (Protein_coding Genes)
#============================================================================================
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")
coding_Biomart <- subset(Biomart$Gene.name, Biomart$Gene.type=="protein_coding")
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_name, REMAP$gene_type=="protein_coding")
coding_PE <- subset(PE$gene, PE$gene_biotype=="protein_coding")
coding_NCBI <- subset(NCBI$gene, NCBI$gene_biotype=="protein_coding")
coding_FANTOM <- subset(FANTOM$gene_symbol, FANTOM$Gene_class=="coding_mRNA")

gene_annot_coding <- union(coding_gencode_superset,
                           union(coding_Biomart,
                                 union(coding_REMAP_gencode_ebi_ftp,
                                       union(coding_PE,
                                             union(coding_NCBI,coding_FANTOM)
                                             )
                                       )
                                 )
                           )
length(gene_annot_coding)

#============================================================================================
write.table(gene_annot_coding, file = "gene_annot_coding.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)
#============================================================================================

#==================================================================================
# 13- Compare gencode_superset & ENSEMBL & Biomart & REMAP_gencode_ebi_ftp (All Genes)
#==================================================================================
x <- list(
  Biomart = Biomart$Gene.name,
  ENSEMBL = ENSEMBL$gene_name,
  gencode_superset = gencode_superset$gene_name,
  REMAP_gencode_ebi_ftp = REMAP$gene_name
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

#=============================================================================
# 14- Compare gencode_superset & Biomart & REMAP_gencode_ebi_ftp & PE (All Genes)
#=============================================================================
x <- list(
  Biomart = Biomart$Gene.name,
  PE = PE$gene,
  gencode_superset = gencode_superset$gene_name,
  REMAP_gencode_ebi_ftp = REMAP$gene_name
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

#==================================================================================
# 15- Compare gencode_superset & Biomart & REMAP_gencode_ebi_ftp & NCBI (All Genes)
#==================================================================================
x <- list(
  Biomart = Biomart$Gene.name,
  NCBI = NCBI$gene,
  gencode_superset = gencode_superset$gene_name,
  REMAP_gencode_ebi_ftp = REMAP$gene_name
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

#==================================================================================
# 16- Compare gencode_superset & Biomart & REMAP_gencode_ebi_ftp & FANTOM (All Genes)
#==================================================================================
x <- list(
  Biomart = Biomart$Gene.name,
  FANTOM = FANTOM$gene_symbol,
  gencode_superset = gencode_superset$gene_name,
  REMAP_gencode_ebi_ftp = REMAP$gene_name
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

#===============================================================================
# 17- Compare gencode_superset & Biomart & REMAP_gencode_ebi_ftp & UCSC (All Genes)
#===============================================================================
y <- list(
  Biomart = Biomart$Gene.name,
  UCSC = UCSC$hg19.kgXref.geneSymbol,
  gencode_superset = gencode_superset$gene_name,
  REMAP_gencode_ebi_ftp = REMAP$gene_name
)

# library(ggplot2)
# library(ggvenn)
ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #fill_color = c("red", "green", "yellow","purple", "green","blue"),
  stroke_size = 0.5, set_name_size = 3
)

U <- setdiff(UCSC$hg19.kgXref.geneSymbol,union(union(Biomart$Gene.name,gencode_superset$gene_name),REMAP$gene_name))
length(U)

#======================================================
# Draw histogram for Gene_type of the annotation data
#======================================================
library(ggplot2)
qplot(Biomart$Gene.type)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

qplot(ENSEMBL$gene_biotype)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

qplot(gencode_superset$gene_type)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

qplot(REMAP$gene_type)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
