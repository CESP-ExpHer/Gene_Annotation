# GRCh37 (hg19) Annotation 
Created by: Yazdan Asgari<br>
Creation date: 7 May 2021<br>
Update: 14 Apr 2022<br><br>

## Here are some information and tips that you can use before working with the annotation files
## Downloading the annotation files
For understanding how to download annotation files related to hg19 genome assembly from various public databases, you could see the [Readme file](0_Download) in the 0_Download folder.
## Analyses of the annotation files
Then, for Comparing the downloaded annotation files, we performed some analyses on them. You could see the [Readme file](1_Analyzing) in the 1_Analyzing folder.
## Creation of the annotation file
Finally, we performesd some cleaning on a selected annotation file in order to be prepared for further usage. We did it for both "Protein-coding" genes and "All" genes in this section. You could see them in the [Readme file](2_Creation_annot_file) in the 2_Creation_annot_file folder.

## How to annotate your SNPs list with an annotation file
There are various scenarios to do that, but here we provide three of them in summary.
### Writing programming script
Suppose you have a list of SNPs (which contains chromosomes and positions columns as well). NOw, you could write a script (for example in R or Python, etc.) to compare every SNP with any gene in the annotation file (its chr, start, and end base pair position) to check if the SNP lie inside the gene or not.<br>
**NOTE 1:** You could consider a specific distance in upstream or downstream of a gene (for instance +/-10Kb) by taking into account this distance into start/end position of the gene.
**NOTE 2:** Because there is biologically possible the genes to have overlaps, you should expect some SNPs might mapped to more than one gene and you should consider in your script how they would be stored in the output file in order to work with them easily in the future.
### Using PLINK

### Using ANNOVAR

