## =============================================================================
## Code for performing some analyses on the downloaded annotation files
## See the GitHub Page for more information:
## https://github.com/theresetruong/0_DATA_REF/blob/main/2_annot_data_ref/Gene_Annotation_Strategy/hg19/1_Analyzing
## Author : Yazdan Asgari, Last Edited:14/06/2021
## https://cesp.inserm.fr/en/equipe/exposome-and-heredity
## =============================================================================

# The folder that the file of Inflammation-Related Genes exists
setwd("PATH_OF_THE_DATA")

# Reading the file
inflame_genes_hgnc <- read.csv('Pathways_genes_hgnc_input.csv',header = TRUE,sep = ",", quote = "\"",
                                         fill = TRUE)
View(inflame_genes_hgnc)
# NOTE:
# inflame_genes_hgnc$geneID refers to the gene IDs appear in the supplementary file of the paper
# inflame_genes_hgnc$hgnc_Approved_symbol refers to the gene symbole extracted from HGNC


# for plotting Venn diagrams
library(ggplot2)
library(ggvenn)

# =============================================================================
# NOTE: To run the following codes, you need to import gene annotation files first
# For more info, see the README file in the "Downloading the annotation files" in the GitHub page
# =============================================================================

# 1- coding_gencode
# Importing the data
library(rtracklayer)
setwd("PATH_OF_THE_DATA")
gencode <- rtracklayer::import('gencode.v19.annotation.gtf')
gencode=as.data.frame(gencode)
coding_gencode <- subset(gencode$gene_name, gencode$gene_type=="protein_coding")
inflame_genes_paper_id <- inflame_genes_hgnc$geneID
inflame_genes_hgnc_id <- inflame_genes_hgnc$hgnc_Approved_symbol
x <- list(
  inflame_genes_paper_id = inflame_genes_paper_id,
  inflame_genes_hgnc_id = inflame_genes_hgnc_id,
  coding_gencode = coding_gencode
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)
n_coding_gencode_paper <- length(intersect(coding_gencode,inflame_genes_paper_id))
n_coding_gencode_paper
n_coding_gencode_paper_not_found = length(inflame_genes_paper_id) - n_coding_gencode_paper
n_coding_gencode_paper_not_found
n_coding_gencode_paper_per = (n_coding_gencode_paper/length(inflame_genes_paper_id))*100
n_coding_gencode_paper_per

n_coding_gencode_hgnc <- length(intersect(coding_gencode,inflame_genes_hgnc_id))
n_coding_gencode_hgnc
n_coding_gencode_hgnc_not_found = length(inflame_genes_hgnc_id) - n_coding_gencode_hgnc
n_coding_gencode_hgnc_not_found
n_coding_gencode_hgnc_per = (n_coding_gencode_hgnc/length(inflame_genes_hgnc_id))*100
n_coding_gencode_hgnc_per

# 2- coding_gencode_superset
# Importing the data
library(rtracklayer)
setwd("PATH_OF_THE_DATA")
gencode_superset <- rtracklayer::import('gencode.v19.chr_patch_hapl_scaff.annotation.gtf')
gencode_superset=as.data.frame(gencode_superset)
coding_gencode_superset <- subset(gencode_superset$gene_name, gencode_superset$gene_type=="protein_coding")
inflame_genes_paper_id <- inflame_genes_hgnc$geneID
inflame_genes_hgnc_id <- inflame_genes_hgnc$hgnc_Approved_symbol
x <- list(
  inflame_genes_paper_id = inflame_genes_paper_id,
  inflame_genes_hgnc_id = inflame_genes_hgnc_id,
  coding_gencode_superset = coding_gencode_superset
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)
n_coding_gencode_paper <- length(intersect(coding_gencode_superset,inflame_genes_paper_id))
n_coding_gencode_paper
n_coding_gencode_paper_not_found = length(inflame_genes_paper_id) - n_coding_gencode_paper
n_coding_gencode_paper_not_found
n_coding_gencode_paper_per = (n_coding_gencode_paper/length(inflame_genes_paper_id))*100
n_coding_gencode_paper_per

n_coding_gencode_hgnc <- length(intersect(coding_gencode_superset,inflame_genes_hgnc_id))
n_coding_gencode_hgnc
n_coding_gencode_hgnc_not_found = length(inflame_genes_hgnc_id) - n_coding_gencode_hgnc
n_coding_gencode_hgnc_not_found
n_coding_gencode_hgnc_per = (n_coding_gencode_hgnc/length(inflame_genes_hgnc_id))*100
n_coding_gencode_hgnc_per

# 3- coding_REMAP_gencode_ebi_ftp
# Importing the data
library(rtracklayer)
setwd("PATH_OF_THE_DATA")
REMAP <- rtracklayer::import('REMAP_gencode.v28lift37.annotation.gtf')
REMAP=as.data.frame(REMAP)
coding_REMAP_gencode_ebi_ftp <- subset(REMAP$gene_name, REMAP$gene_type=="protein_coding")
inflame_genes_paper_id <- inflame_genes_hgnc$geneID
inflame_genes_hgnc_id <- inflame_genes_hgnc$hgnc_Approved_symbol
y <- list(
  inflame_genes_paper_id = inflame_genes_paper_id,
  inflame_genes_hgnc_id = inflame_genes_hgnc_id,
  coding_REMAP_gencode_ebi_ftp = coding_REMAP_gencode_ebi_ftp
)
ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)
n_coding_gencode_paper <- length(intersect(coding_REMAP_gencode_ebi_ftp,inflame_genes_paper_id))
n_coding_gencode_paper
n_coding_gencode_paper_not_found = length(inflame_genes_paper_id) - n_coding_gencode_paper
n_coding_gencode_paper_not_found
n_coding_gencode_paper_per = (n_coding_gencode_paper/length(inflame_genes_paper_id))*100
n_coding_gencode_paper_per

n_coding_gencode_hgnc <- length(intersect(coding_REMAP_gencode_ebi_ftp,inflame_genes_hgnc_id))
n_coding_gencode_hgnc
n_coding_gencode_hgnc_not_found = length(inflame_genes_hgnc_id) - n_coding_gencode_hgnc
n_coding_gencode_hgnc_not_found
n_coding_gencode_hgnc_per = (n_coding_gencode_hgnc/length(inflame_genes_hgnc_id))*100
n_coding_gencode_hgnc_per

# 4- coding_REMAP_latest
# Importing the data
library(rtracklayer)
setwd("PATH_OF_THE_DATA")
REMAP_latest <- rtracklayer::import('REMAP_latest_gencode.v38lift37.annotation.gtf')
REMAP_latest=as.data.frame(REMAP_latest)
coding_REMAP_latest <- subset(REMAP_latest$gene_name, REMAP_latest$gene_type=="protein_coding")
inflame_genes_paper_id <- inflame_genes_hgnc$geneID
inflame_genes_hgnc_id <- inflame_genes_hgnc$hgnc_Approved_symbol
y <- list(
  inflame_genes_paper_id = inflame_genes_paper_id,
  inflame_genes_hgnc_id = inflame_genes_hgnc_id,
  coding_REMAP_latest = coding_REMAP_latest
)
ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)
n_coding_gencode_paper <- length(intersect(coding_REMAP_latest,inflame_genes_paper_id))
n_coding_gencode_paper
n_coding_gencode_paper_not_found = length(inflame_genes_paper_id) - n_coding_gencode_paper
n_coding_gencode_paper_not_found
n_coding_gencode_paper_per = (n_coding_gencode_paper/length(inflame_genes_paper_id))*100
n_coding_gencode_paper_per

n_coding_gencode_hgnc <- length(intersect(coding_REMAP_latest,inflame_genes_hgnc_id))
n_coding_gencode_hgnc
n_coding_gencode_hgnc_not_found = length(inflame_genes_hgnc_id) - n_coding_gencode_hgnc
n_coding_gencode_hgnc_not_found
n_coding_gencode_hgnc_per = (n_coding_gencode_hgnc/length(inflame_genes_hgnc_id))*100
n_coding_gencode_hgnc_per

# 5- gencode_all
# Note that we have importing the annotation data before (in the coding section of the current code)
inflame_genes_paper_id <- inflame_genes_hgnc$geneID
inflame_genes_hgnc_id <- inflame_genes_hgnc$hgnc_Approved_symbol
x <- list(
  inflame_genes_paper_id = inflame_genes_paper_id,
  inflame_genes_hgnc_id = inflame_genes_hgnc_id,
  gencode_all = gencode$gene_name
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)
n_gencode_paper <- length(intersect(gencode$gene_name,inflame_genes_paper_id))
n_gencode_paper
n_gencode_paper_not_found = length(inflame_genes_paper_id) - n_gencode_paper
n_gencode_paper_not_found
n_gencode_paper_per = (n_gencode_paper/length(inflame_genes_paper_id))*100
n_gencode_paper_per

n_gencode_hgnc <- length(intersect(gencode$gene_name,inflame_genes_hgnc_id))
n_gencode_hgnc
n_gencode_hgnc_not_found = length(inflame_genes_hgnc_id) - n_gencode_hgnc
n_gencode_hgnc_not_found
n_gencode_hgnc_per = (n_gencode_hgnc/length(inflame_genes_hgnc_id))*100
n_gencode_hgnc_per

gencode_superset <- rtracklayer::import('gencode.v19.chr_patch_hapl_scaff.annotation.gtf')
gencode_superset=as.data.frame(gencode_superset)

# 6- gencode_superset_all
# Note that we have importing the annotation data before (in the coding section of the current code)
inflame_genes_paper_id <- inflame_genes_hgnc$geneID
inflame_genes_hgnc_id <- inflame_genes_hgnc$hgnc_Approved_symbol
x <- list(
  inflame_genes_paper_id = inflame_genes_paper_id,
  inflame_genes_hgnc_id = inflame_genes_hgnc_id,
  gencode_superset_all = gencode_superset$gene_name
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)
n_gencode_paper <- length(intersect(gencode_superset$gene_name,inflame_genes_paper_id))
n_gencode_paper
n_gencode_paper_not_found = length(inflame_genes_paper_id) - n_gencode_paper
n_gencode_paper_not_found
n_gencode_paper_per = (n_gencode_paper/length(inflame_genes_paper_id))*100
n_gencode_paper_per

n_gencode_hgnc <- length(intersect(gencode_superset$gene_name,inflame_genes_hgnc_id))
n_gencode_hgnc
n_gencode_hgnc_not_found = length(inflame_genes_hgnc_id) - n_gencode_hgnc
n_gencode_hgnc_not_found
n_gencode_hgnc_per = (n_gencode_hgnc/length(inflame_genes_hgnc_id))*100
n_gencode_hgnc_per


# 7- REMAP_gencode_ebi_ftp_all
# Note that we have importing the annotation data before (in the coding section of the current code)
inflame_genes_paper_id <- inflame_genes_hgnc$geneID
inflame_genes_hgnc_id <- inflame_genes_hgnc$hgnc_Approved_symbol
y <- list(
  inflame_genes_paper_id = inflame_genes_paper_id,
  inflame_genes_hgnc_id = inflame_genes_hgnc_id,
  REMAP_gencode_ebi_ftp_all = REMAP$gene_name
)
ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)
n_gencode_paper <- length(intersect(REMAP$gene_name,inflame_genes_paper_id))
n_gencode_paper
n_gencode_paper_not_found = length(inflame_genes_paper_id) - n_gencode_paper
n_gencode_paper_not_found
n_gencode_paper_per = (n_gencode_paper/length(inflame_genes_paper_id))*100
n_gencode_paper_per

n_gencode_hgnc <- length(intersect(REMAP$gene_name,inflame_genes_hgnc_id))
n_gencode_hgnc
n_gencode_hgnc_not_found = length(inflame_genes_hgnc_id) - n_gencode_hgnc
n_gencode_hgnc_not_found
n_gencode_hgnc_per = (n_gencode_hgnc/length(inflame_genes_hgnc_id))*100
n_gencode_hgnc_per


# 8- REMAP_latest_all
# Note that we have importing the annotation data before (in the coding section of the current code)
inflame_genes_paper_id <- inflame_genes_hgnc$geneID
inflame_genes_hgnc_id <- inflame_genes_hgnc$hgnc_Approved_symbol
y <- list(
  inflame_genes_paper_id = inflame_genes_paper_id,
  inflame_genes_hgnc_id = inflame_genes_hgnc_id,
  REMAP_latest_all = REMAP_latest$gene_name
)
ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)
n_gencode_paper <- length(intersect(REMAP_latest$gene_name,inflame_genes_paper_id))
n_gencode_paper
n_gencode_paper_not_found = length(inflame_genes_paper_id) - n_gencode_paper
n_gencode_paper_not_found
n_gencode_paper_per = (n_gencode_paper/length(inflame_genes_paper_id))*100
n_gencode_paper_per

n_gencode_hgnc <- length(intersect(REMAP_latest$gene_name,inflame_genes_hgnc_id))
n_gencode_hgnc
n_gencode_hgnc_not_found = length(inflame_genes_hgnc_id) - n_gencode_hgnc
n_gencode_hgnc_not_found
n_gencode_hgnc_per = (n_gencode_hgnc/length(inflame_genes_hgnc_id))*100
n_gencode_hgnc_per

# =============================================================================
# Find and Save Infammation Gene symbols that could NOT found in the database
# =============================================================================
U <- setdiff(inflame_genes_hgnc$hgnc_Approved_symbol,
             intersect(inflame_genes_hgnc$hgnc_Approved_symbol,REMAP_latest$gene_name)
             )
length(U)
U
write.table(U, file = "Not_found_inflame_genes.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)


# =============================================================================
# Find the gene class of Infammation Genes that found in the REMAP_latest_all database
# =============================================================================
inf_genes_found <- intersect(inflame_genes_hgnc$hgnc_Approved_symbol,REMAP_latest$gene_name)
class(inf_genes_found)

A <- REMAP_latest[unlist(lapply(inf_genes_found, function(x) match(x, REMAP_latest$gene_name))),]

qplot(A$gene_type)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

write.table(A, file = "found_inflame_genes.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)


