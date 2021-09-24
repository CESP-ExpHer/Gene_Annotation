## =============================================================================
## Some Analyses on the "REMAP_latest_gencode.v38lift37.annotation.gtf" file
## See the GitHub Page for more information:
## https://github.com/theresetruong/0_DATA_REF/blob/main/2_annot_data_ref/Gene_Annotation_Strategy/hg19/1_Analyzing
## Author : Yazdan Asgari, Last Edited:14/06/2021
## https://cesp.inserm.fr/en/equipe/exposome-and-heredity
## =============================================================================

# Path where the "REMAP_latest_gencode.v38lift37.annotation.gtf" file exists
setwd("PATH_OF_DOWNLOADED_ANNOTATION_DATA")

# =============================================================================
# Total number of genes 
genes_REMAP_latest <- subset(REMAP_latest$gene_name, REMAP_latest$type=="gene") 
genes_REMAP_latest <- REMAP_latest[which(REMAP_latest$type=="gene"),]
length(genes_REMAP_latest$gene_name)

# Total number of unique genes (JUST by looking at gene name) 
length(unique(genes_REMAP_latest$gene_name))

# Total number of unique genes (by looking at gene name & chromosome) 
genes_uniques_REMAP_latest <- REMAP_latest[!duplicated(REMAP_latest[c("gene_name", "seqnames")]),] 
length(genes_uniques_REMAP_latest$gene_name)

# =============================================================================
# Total number of protein_coding genes
genes_coding_REMAP_latest <- subset(REMAP_latest$gene_name, ( REMAP_latest$type=="gene"  & 
                                                                REMAP_latest$gene_type=="protein_coding" 
)
)
length(genes_coding_REMAP_latest)

# Total number of unique genes (JUST by looking at gene name) 
length(unique(genes_coding_REMAP_latest$gene_name))

# Total number of unique genes (by looking at gene name & chromosome) 
genes_coding_REMAP_latest <- REMAP_latest[which(REMAP_latest$type=="gene"  & 
                                                  REMAP_latest$gene_type=="protein_coding" 
                                                ),
                                          ]
genes_coding_rem_dup_REMAP_latest <- genes_coding_REMAP_latest[!duplicated(genes_coding_REMAP_latest[c("gene_name", "seqnames")]),] 
length(genes_coding_rem_dup_REMAP_latest$gene_name)

# =============================================================================
# Show the duplicate found based on gene_name column (show ALL rows including duplicates as well)
A1 <- genes_coding_REMAP_latest[duplicated(genes_coding_REMAP_latest["gene_name"]) | duplicated(genes_coding_REMAP_latest[,"gene_name"], fromLast=TRUE),]
A1 <- A1[order(A1$gene_name),]

# Show the duplicate found based on gene_name column (JUST show one duplicate)
A2 <- genes_coding_REMAP_latest[duplicated(genes_coding_REMAP_latest["gene_name"]),]
A2 <- A2[order(A2$gene_name),]

# Show the duplicate found based on gene_name & chromosome columns (show ALL rows including duplicates as well)
A3 <- genes_coding_REMAP_latest[duplicated(genes_coding_REMAP_latest[c("gene_name", "seqnames")]) | duplicated(genes_coding_REMAP_latest[c("gene_name", "seqnames")], fromLast=TRUE),]
# Show the duplicate found based on gene_name & chromosome columns
A4 <- genes_coding_REMAP_latest[duplicated(genes_coding_REMAP_latest[c("gene_name", "seqnames")]),] 

# Remove duplicates based on gene_name & chromosome columns
A5 <- genes_coding_REMAP_latest[!duplicated(genes_coding_REMAP_latest[c("gene_name", "seqnames")]),] 
# =============================================================================

