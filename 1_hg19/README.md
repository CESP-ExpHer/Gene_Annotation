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
### 1- Writing a programming script
Suppose you have a list of SNPs (which contains chromosomes and positions columns as well). NOw, you could write a script (for example in R or Python, etc.) to compare every SNP with any gene in the annotation file (its chr, start, and end base pair position) to check if the SNP lie inside the gene or not.<br><br>
**NOTE 1:** You could consider a specific distance in upstream or downstream of a gene (for instance +/-10Kb) by taking into account this distance into start/end position of the gene.<br><br>
**NOTE 2:** Because there is biologically possible the genes to have overlaps, you should expect some SNPs might mapped to more than one gene and you should consider in your script how they would be stored in the output file in order to work with them easily in the future.

### 2- Using PLINK
You could use **PLINK** tool in order to annotate your SNPs list [Plink web page](https://zzz.bwh.harvard.edu/plink/annot.shtml)
Here, we provide you a summary. Assume you have a SNPs list file called ***GWAS.txt***.<br>
This file could be used as an input for Plink annotation process.<br>
So, prepare a file called ***GWAS.assoc*** with the following four columns:<br>
 | CHR | SNP | BP | P |
 | --- | ---------- | ---------- | ---------- |
 |  11  | rs2212434  | 76281593 | 1E-10 |
 |  11  | rs61893460 | 76291154 | 1E-9 |

It is also needed to have an annotation list file that you prepared before (as we talked about it earlier in this page) (suppose we called it ***annotation_hg19.txt***). <br>
Another file called ***snp129.attrib*** is needed for the annotation. <br>
You could download it directly from the [Plink web page](https://zzz.bwh.harvard.edu/plink/res.shtml#attrib).<br>
Then, all you need is to run the following command for annotation of the GWAS file:<br>
```
plink --annotate GWAS.assoc attrib=snp129.attrib.gz ranges=annotation_hg19.txt
```
**Note 1:** All files should be in the same folder.<br><br>
**Note 2:** To specify a particular distance for genes/ranges (i.e. upstream/downstream of a gene), 
use the following command, e.g. for for 1000 bp (1Kb) upstream and downstream
```
plink --annotate GWAS.assoc attrib=snp129.attrib.gz ranges=annotation_hg19.txt --border 1
```
After running, an output called ***plink.annot*** would be created. <br>
Here is some lines of the annotation output:<br>
| CHR | SNP | BP | P | ANNOT |
| --- | ---------- | ---------- | ---------- | ---------- |
|  11  | rs2212434  | 76281593 | 1E-9 | . |
|  11  | rs7931483 | 76302067 | 1E-9 | 'RP11-672A2.7'(0) |
|  12  | rs3024971 | 57493727 | 1E-9 | 'CATG00000009161.1'(0)\|'STAT6'(0) |

**Note 3:** Always check that the line numbers of this file should be the same as the ***GWAS.assoc*** line numbers.<br><br>
**Note 4:** The value in the parenthesis show if a SNP lie inside a gene (0) or outside of it (non-zero).<br><br>
**Note 5:** If a SNP overlapped with more than one gene, they separated by "|" character.

### 3- Using ANNOVAR
You could use "ANNOVAR" tool to annotate a list of SNPs. It gives you multiple output files including annotation for each SNP, detailed information for exonic SNPs, if a SNP got error during the annotation process, and a log file. <br>
You could find how to use this tool in our GitHub page [ANNOVAR GitHUb Page](/../../Tutorial_ANNOVAR)

