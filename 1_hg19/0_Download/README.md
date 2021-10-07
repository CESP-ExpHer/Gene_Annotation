# Downloading annotation files related to hg19 genome assembly from public databases
Created by: Yazdan Asgari<br>
First creation date: 7 May 2021<br>
Last Updated: 7 Oct 2021<br><br>
**NOTE:** You could find all downloaded annotation files (explained in this section) in the **"2_Gene_annotation"** folder on the server.

- [Downloading annotation files related to hg19 genome assembly from public databases](#downloading-annotation-files-related-to-hg19-genome-assembly-from-public-databases)
  * [Overall Procedure](#overall-procedure)
    + [GENCODE Database](#gencode-database)
    + [ENSEMBL Database](#ensembl-database)
    + [BioMart Tool](#biomart-tool)
    + [UCSC Database](#ucsc-database)
    + [Gene Annotation used in one of our papers](#gene-annotation-used-in-one-of-our-papers)
    + [FANTOM Database](#fantom-database)
    + [Gene Annotation in the Server](#gene-annotation-in-the-server)
    + [NCBI Database](#ncbi-database)
    + [Summary of 10 downloaded gene annotation files](#summary-of-10-downloaded-gene-annotation-files)

## Overall Procedure
As you can see in the following image, we wanted to download 10 different annotation data for hg19 genome assembly from public databases.
<br></br>
<kbd> <img src="0_Images/Slide2.PNG"/> </kbd>
### GENCODE Database
Based on the GENCODE webpage, *"The goal of the GENCODE project is to identify and classify all gene features in the human and mouse genomes with high accuracy based on biological evidence, and to release these annotations for the benefit of biomedical research and genome interpretation."*<br>
If you click on the **"Human"** Tab in the website, it is possible to use **"Current release"** for the Human.
<br></br>
<kbd> <img src="0_Images/Slide5.PNG"/> </kbd>
<br></br>
Here, there is files for the lates release. However, they are for **hg38** assembly version. So, you could click on **"Go to GRCh37 version of this release"**.
<br></br>
<kbd> <img src="0_Images/Slide6.PNG"/> </kbd>
<br></br>
In this page, you will see the **GRCh37** data. But if you note, it is mentioned that *"gene annotation originally created on the GRCh38 reference chromosomes, mapped to the GRCh37 primary assembly with gencode-backmap"*. So, it is not what we are looking for.
<br></br>
<kbd> <img src="0_Images/Slide7.PNG"/> </kbd>
<br></br>
Instead, from the **"Human"** tab menu, you should click on **"Release history"**:
<br></br>
<kbd> <img src="0_Images/Slide8.PNG"/> </kbd>
<br></br>
Then, you could select the **GRCh37** assembly (which is GENCODE release 19):
<br></br>
<kbd> <img src="0_Images/Slide9.PNG"/> </kbd>
<br></br>
So, in this page, you could download gene annotation file in a *GTF* file format.<br>
As you can see, there are two files: the *main annotation* file AND the *superset of the main annotation* file. We downloaded both files for further analyses.
<br></br>
<kbd> <img src="0_Images/Slide10.PNG"/> </kbd>
<br></br>
Now, we have two annotation files downloaded:
<br></br>
<kbd> <img src="0_Images/Slide11.PNG"/> </kbd>
<br></br>
### ENSEMBL Database
As you know, Ensembl is a data resource for the European Bioinformatics Institute (EMBL-EBI).
<br></br>
<kbd> <img src="0_Images/Slide12.PNG"/> </kbd>
<br></br>
In the Ensemble webpage, you need to click on the **"Human"** option.
<br></br>
<kbd> <img src="0_Images/Slide13.PNG"/> </kbd>
<br></br>
Then, in the **"Gene annotation"** box, you could click on **"Download GTF or GFF3"** to download the annotation files. But we know that this data is for the **GRCh38** assembly.
<br></br>
<kbd> <img src="0_Images/Slide14.PNG"/> </kbd>
<br></br>
So, in the **"Genome assembly"** box, use the **"Other assemblies"** drop-down menu to select the **hg19** assembly:
<br></br>
<kbd> <img src="0_Images/Slide15.PNG"/> </kbd>
<br></br>
Similar to the previous assembly page, we expected to obtain the annotation file by clicking on the **"Gene annotation"** box. But the links in this page does NOT work:
<br></br>
<kbd> <img src="0_Images/Slide16.PNG"/> </kbd>
<br></br>
So, you need to use another approach:<br>
This time, start from the **CRCh38** assembly address and click on **"../"** part of the page to go to the upper levels (called parent folder) of the ftp pages:
<br></br>
<kbd> <img src="0_Images/Slide17.PNG"/> </kbd>
<br></br>
When you reached the **"/pub"** folder, click on the **"grch37"** text, then **"release-104"**, then **"gtf"**, then **"homo_sapiens"**. Now, you can download the annotation file for **hg19** assembly.
<br></br>
<kbd> <img src="0_Images/Slide18.PNG"/> </kbd>
<br></br>
So, we have downloaded 3 annotation files so far:
<br></br>
<kbd> <img src="0_Images/Slide19.PNG"/> </kbd>
<br></br>
### BioMart Tool
If you notice to the Ensembl webpage, you could see **"BioMart"** at the top of the screen.
<br></br>
<kbd> <img src="0_Images/Slide20.PNG"/> </kbd>
<br></br>
In the Wikipedia, you see that *"BioMart is a community-driven project to provide a single point of access to distributed research data."*<br>
You could select different datasets through **"CHOOSE DATABASE"** drop-down menu in the BioMart page to download different data.
<br></br>
<kbd> <img src="0_Images/Slide21.PNG"/> </kbd>
<br></br>
But as you know, this data is for the **GRCh38** assembly. For access to the **GRCh37** assembly, click on the following link:<br>
[BioMaer GRCh37 assembly link](https://grch37.ensembl.org/biomart/martview/4f4cea56e8f6f9bbe161ba1f6c97d923)
<br></br>
<kbd> <img src="0_Images/Slide22.PNG"/> </kbd>
<br></br>
Now, in this page, select **"Ensembl Genes 104"**, then **"Human genes (GRCh37.p13)"**.<br>
If you click on the **"Count"** button, it tells you the number of genes in the Dataset.
<br></br>
<kbd> <img src="0_Images/Slide23.PNG"/> </kbd>
<br></br>
You could set **"Filters"** and adds **"Attributes"** to the dataset to choose what you want to be downloaded.<br>
Then, click on the **"Results"** button to see 10 lines of the result.<br>
For download the data, after choosing the output file formar, click on the **"Go"** button.
<br></br>
<kbd> <img src="0_Images/Slide24.PNG"/> </kbd>
<br></br>
This is the fourth annotation file:
<br></br>
<kbd> <img src="0_Images/Slide25.PNG"/> </kbd>
<br></br>
But if you see next to the **"Go"** button, there is an option called **"Unique results only"**. So, selection of that option, we downloaded the next annotation data:
<br></br>
<kbd> <img src="0_Images/Slide26.PNG"/> </kbd>
<br></br>
<br></br>
<kbd> <img src="0_Images/Slide27.PNG"/> </kbd>
<br></br>
### UCSC Database
In the UCSC database, there is a **"Table Browser"** tool used to download the data:
<br></br>
<kbd> <img src="0_Images/Slide28.PNG"/> </kbd>
<br></br>
In the Table browser, by default, there is not a Gene symbol column for downloaded data. For adding Gene Symbol, you should link the official gene symbol via a linkage between the knownGene table and the *kgXref table* in the human assembly as follows:<br>
Select the following parameters:<br>
**group:** Genes and Gene Predictions<br>
**track:** UCSC Genes<br>
**table:** knownGene<br>
Then, using the **"paste list"** button, you should upload list of UCSC IDs into the Table Browser.<br> 
**output format:** selected fields from primary and related tables<br>
Then press the **“get output”** button.<br>
On this page, you can select **“geneSymbol”** from **“hg19.kgXref fields”** (and other fields that you like)<br>
Then, press the **“get output”** button again.
<br></br>
<kbd> <img src="0_Images/Slide29.PNG"/> </kbd>
<br></br>
So, we downloaded 6 annotation files so far:
<br></br>
<kbd> <img src="0_Images/Slide30.PNG"/> </kbd>
<br></br>
### Gene Annotation used in one of our papers
The next file is mentioned in the following paper:<br>
*Gene network and biological pathways associated with susceptibility to differentiated thyroid carcinoma*
<br>
*Kulkami et al., Sci Rep., 2021 Apr 26;11(1):8932. doi: 10.1038/s41598-021-88253-0*
<br>
<br></br>
<kbd> <img src="0_Images/paper1.png"/> </kbd>
<br></br>
And here is the address for the annotation file:<br>
[Annotation file Link](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/)
<br></br>
<kbd> <img src="0_Images/Slide31.PNG"/> </kbd>
<br></br>
You need to download **"gencode.v28lift37.annotation.gtf.gz"** file.<br>
**NOTE:** if you see the *"README_GRCh37_mapping.txt"*, you find that this file contain GENCODE annotation originally built on GRCh38 that has been *mapped back* to the GRCh37 assembly.
So, we have the 7th annotation file:
<br></br>
<kbd> <img src="0_Images/Slide32.PNG"/> </kbd>
<br></br>
### FANTOM Database
For FANTOM annotation file, you need to go to the FANTOM project webpage:
<br></br>
<kbd> <img src="0_Images/fantom0.png"/> </kbd>
<br></br>
Then go to the **"Data"**: 
<br></br>
<kbd> <img src="0_Images/fantom1.png"/> </kbd>
<br></br>
In this page, click on the **"Hon_et_al_2016"** data, then **"assembly"**, then **"lv3_robust"**:
<br></br>
<kbd> <img src="0_Images/fantom2.png"/> </kbd>
<br></br>
<br></br>
<kbd> <img src="0_Images/fantom3.png"/> </kbd>
<br></br>
Then, you need to download **"FANTOM_CAT.lv3_robust.info_table.gene.tsv.gz"** file.
<br></br>
<kbd> <img src="0_Images/fantom4.png"/> </kbd>
<br></br>
So, we have the 8th annotation file:
<br></br>
<kbd> <img src="0_Images/Slide34.PNG"/> </kbd>
<br></br>
### Gene Annotation in the Server
Also, there is a file in our server called **"GRCh37_latest_genomic.gff"**:
<br></br>
<kbd> <img src="0_Images/35_new.png"/> </kbd>
<br></br>
So, we have the 9th annotation file:
<br></br>
<kbd> <img src="0_Images/Slide36.PNG"/> </kbd>
<br></br>
### NCBI Database
You could find the annotation file from the NCBI database based on the following path:<br>
From the NCBI Main page **/ Genome / Human Genome / GRCh37 / gff3**
<br></br>
<kbd> <img src="0_Images/Slide37.PNG"/> </kbd>
<br></br>
<br></br>
<kbd> <img src="0_Images/Slide38.PNG"/> </kbd>
<br></br>
<br></br>
<kbd> <img src="0_Images/Slide39.PNG"/> </kbd>
<br></br>
Downloading the gff3 file, we have the 10th annotation file:
<br></br>
<kbd> <img src="0_Images/Slide40.PNG"/> </kbd>
<br></br>
### Summary of 10 downloaded gene annotation files
So, we have 10 different annotation files, all for **GRCh37 (hg19)** assembly of the Human Genome from the most important public databases.<br>
For more convenience, we add the name of the database as the starting string for each downloaded file:
<br></br>
<kbd> <img src="0_Images/Slide83.PNG"/> </kbd>
<br></br>
Now you can go to the [Analysis Section](../1_Analyzing) to see the similarity and differences of these annotation files.
