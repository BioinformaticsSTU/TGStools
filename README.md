# TGStools
TGStools is a bioinformatics suit to facilitate transcriptome analysis of long reads from third generation sequencing platform.


----------------------------
# Table of Contents
----------------------------

   * [Overview](#overview)
   * [Installation](#installation)
   * [Command and subcommand structure](#command-and-subcommand-structure)
      * [geneDisplay](#geneDisplay)
      * [staDist](#staDist)
      * [INCP](#INCP)
      * [extract_lncRNA_gtf](#extract_lncrna_gtf)
      * [tiss_specific](#tiss_specific)
      * [staAS](#staas)
      * [calScoreD](#calscored)
      * [GOenrich](#goenrich)

   * [License](#license)




----------------------------
# Overview
----------------------------

The third generation sequencing can de novo detect long reads of several thousand base pairs, thus provides a global view of the full length transcriptome. But due to less sequencing accurate rate, it often yields many spurious transcripts. It's important to prioritize the results by a visualization framework that automatically integrates rich annotation information with known genomic features. Therefore, we developed TGStools, a bioinformatics suit to facilitate routine tasks such as characterizing the full-length transcripts and detecting the shifted types of alternative splicing in post transcriptome analysis.


----------------------------
# Installation
----------------------------

TGStools has been developed in Python 3.5. 

If necessary, to install python3 we recommend to download from the official site https://www.python.org/downloads/ the corresponding version for your OS.


----------------------------
# Command and subcommand structure
----------------------------

TGStools works with a command/subcommand structure:

```
python3 TGStools.py subcommand options

```
where the subcommand can be one of these:

- **geneDisplay**    :
- **staDist**    :
- **INCP**    : an integration classification tool of CNCI and PLEK for identify coding or non-coding transcripts (fasta file and gtf file).
- **extract_lncRNA_gtf**       : A tool that extract lncRNA information of GTF format based on the tanscript ID of the candidate lncRNA.
- **tiss_specific**       : A tool that extract cancer-specific lncRNA information of GTF format.
- **staAS**        : calculate the proportion of each alternative splicing event in different samples and create graphs.
- **calScoreD**     : calculate score_D of each gene.
- **GOenrich**     : select top genes and make GO enrichment analysis.


----------------------------
## geneDisplay
----------------------------

input GTF files, gene id and trans_quant files(quantity of transcripts), then the macroscopic image showing the expression of transcripts of queried gene can be gotten. Similar to UCSC's gene query, users can choose to input multiple epigenetic data in the same folder, and make statistical drawings of the transcripts and epigenetic data under the gene.

### Input files

#### gtf file
An annotation file in GTF format is like:

```
chr14 Ensembl exon  73741918  73744001  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73749067  73749213  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";  
chr14 Ensembl exon  73750789  73751082  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73753818  73754022  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
```

#### path of histone files or fatom5 files

#### trans_quant
```
i2_LQ_K5103rd_c89674/f1p10/3040	ENSG00000230183_novel01	ENSG00000230183
i3_LQ_K5103rd_c15439/f1p0/3034	ENSG00000230183_novel01	ENSG00000230183
i3_LQ_K5103rd_c13337/f1p5/3578	ENST00000489294	ENSG00000152332
i2_LQ_K5103rd_c68909/f4p3/2516	ENST00000489294	ENSG00000152332
i2_LQ_K5103rd_c52655/f1p22/2851	ENST00000489294	ENSG00000152332
```

### Usage
```
python3 TGStools.py geneDisplay -g <gtf> -i <gene_id> -q <trans_quant> -p <path>
```

- **-g**  | **--gtf**: gtf file

- **-i**  | **--id**: gene id thao you want to query

- **-q**  | **--quant**: quantity of transcript

- **-p**  | **--path**: path of histone files or fatom5 files

### Example
```
python3 TGStools.py geneDisplay -g K510_3rd.gtf  -i ENSG00000035141 -q K510_3rd.id2id.xls  -p histone
```

### Output files


#### display of transcripts
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/ENSG00000035141.png" width = "550" height = "400"  />
ENSG00000035141.png


----------------------------
## staDist
----------------------------

statistical tables and histone/fantom5 histone files

### Input files

#### gtf file

#### path of histone files or fatom5 files

### Usage
```
python3 TGStools.py staDist -g <gtf> -p <path> -f <flag>
```

- **-g**  | **--gtf**: input file of fasta file or gtf file, if the input is fasta file,the file format must be the twolineFasta

- **-p**  | **--path**: path of histone files or fatom5 files

- **-f**  | **--flag**: flag that tells programme the type of files in path


### Example
```
python3 TGStools.py staDist -g K510_3rd.gtf  -p histone  -f histone
or
python3 TGStools.py staDist -g K510_3rd.gtf  -p fantom5  -f fantom5
```



----------------------------
## INCP
----------------------------

an integration classification tool of CNCI and PLEK for identify coding or non-coding transcripts (fasta file and gtf file)
**This step will cost much time and resource, we strongly recommend user to run this step at night while the computer is free.**

### Input files

#### hg19.2bit and hg38.2bit
hg19.2bit and hg38.2bit are 2bit format of human genomes which can be downloaded on UCSC.

#### gtf file
An annotation file in GTF format is like:

```
chr14 Ensembl exon  73741918  73744001  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73749067  73749213  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";  
chr14 Ensembl exon  73750789  73751082  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73753818  73754022  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
```
Note: chromosome id must have prefix chr.

#### fasta file
An annotation file in FASTA format is like:
```
> SEQ1
AGCACCAGCCGCCCCATCGCCACCGCCGCCGCCGCCGCCCGGATCCTGGCGCGCTGAATGCAGACTAACA
> SEQ2
CCATAAAGACTTCAAGGAACTAAGGTACAATAAGTGTCTTATGAACTTCAGCTGCAATGGAAAGAATGGAAGCT
```
Note: fasta file format must be the twolineFasta

### Usage
```
python3 TGStools.py INCP -i <file> -p <parallel> -r <reference> -g
or 
python3 TGStools.py INCP -i <file> -p <parallel>
```

- **-i**  | **--input**: input file of fasta file or gtf file, if the input is fasta file,the file format must be the twolineFasta

- **-p**  | **--parallel**: assign the running CUP numbers

- **-g**  | **--gtf**: if your input file is gtf format please use this parameter

- **-r**  | **--reference**: if you use the -g this parameter must be assigned, within this parameter please assign the path of your reference genome. Some reference files which has been prepared could be download at [hg38](hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit), [hg19](hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit).

### Example
```
python3 TGStools.py INCP -i candidate.gtf -p 6 -g -d hg38.2bit
or 
python3 TGStools.py INCP -i candidate.fasta -p 6
```

### Output files

mainly contains 4 files and 1 directory

#### input_no_suffix directory
directory of the CNCI output result.

#### input_no_suffix_PLEK
directory of the PLEK output result.

#### input_no_suffix_union_plek_cnci.txt
output of union of the software CNCI and PLEK, in which the first column is the tanscript ID

#### input_no_suffix_intersect_plek_cnci.txt
output of intersect of the software CNCI and PLEK, in which the first column is the tanscript ID

#### venny_plek_cnci.pdf 
the summary of the venny between the output of the CNCI and PLEK

<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/venn.png" width = "550" height = "400"  />

venny image

----------------------------
## extract_lncRNA_gtf
----------------------------

A tool that extract lncRNA information of GTF format based on the tanscript ID of the candidate lncRNA. This tool is only for result of gtf mode.

### Input files

#### index index
input file of the candidate lncRNA, in which the first column is the tanscript ID of the candidate lncRNA.
```
ENSG00000174151_novel02 noncoding       -0.047104       2124    2604    2758
ENSG00000118482_novel01 coding  0.084   0       1512    2792
ENSG00000118482_novel02 coding  0.073   0       1455    2120
ENSG00000118482_novel04 coding  0.064   0       1695    2979
ENSG00000118482_novel06 noncoding       -0.047104       2718    2922    3078
ENSG00000118482_novel07 noncoding       -0.0217088      2166    2355    2510
```
#### gtf file

### Usage
```
python3 TGStools.py extract_lncRNA_gtf -i <file> -g <gtf> -o <out>
```

- **-i**  | **--input**: input file of the candidate lncRNA, in which the first column is the tanscript ID of the candidate lncRNA. This file also can be the output file of INCP

- **-g**  | **--gtf**: GTF file corresponding to fasta in the incp, where the last column contain the tanscript ID

- **-o**  | **--out**: output name extracted lncRNA information of GTF format

### Example
```
python3 TGStools.py extract_lncRNA_gtf -i intersect_plek_cnci.txt -g unannotation.gtf -o out
```

### Output files
output file extract lncRNA information of GTF format

----------------------------
## tiss_specific
----------------------------

A tool that extract cancer-specific lncRNA information of GTF format.

**If you want to run this step faster, skip the parameter -t.**

### Input files

#### FILE 
files of the candidate lncRNA gtf format.

#### TISSUE FILE


### Usage
```
python3 TGStools.py tiss_specific -i <file>[,<file>] [-t <tss>] -r <reference> -o <out>
```

- **-i**  | **--input**: input files of the candidate lncRNA gtf format, if the input files have two splited by ',', the first set control sample, the other set cancer sample. If the input files have only one, there are two situations. The one is based on a background control tissue. The other have not a background control tissue. Related knowledge can refer to the  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3185964/

- **-t**  | **--tissue**: set tissue name if input files have the only one. If there are two input files, this parameter can be ignored.The background control tissue : 'adipose', 'adrenal', 'brain', 'breast', 'colon', 'heart', 'kidney', 'liver', 'lung', 'lymphNode', 'ovary', 'prostate', 'skeltalMuscle', 'whiteBloodCell', 'testes', 'thyroid', 'placenta', 'foreskin', 'hLF'.

- **-r**  | **--reference**: set refgene name,'hg38' and 'hg19' can be chosen.

- **-o**  | **--out**: output name extracted lncRNA-specific information of GTF format

### Example
```
python3 TGStools.py tiss_specific -f sample.gtf[,control.gtf] -t breast -o out.gtf
```

### Output files
GTF file which extract lncRNA-specific information

----------------------------
## staAS
----------------------------
calculate the proportion of each alternative splicing event in different samples and create graphs.
### Input files

#### gtf file
An annotation file in GTF format is required:

```
chr14 Ensembl exon  73741918  73744001  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73749067  73749213  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";  
chr14 Ensembl exon  73750789  73751082  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73753818  73754022  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
```

### Usage
```
python3 TGStools.py staAS -n <names> -p <prefix>
```
List of options available:

- **-n**  | **--names**: names of GTF format file/files containing at least "exon" lines

- **-p**  | **--prefix**: prefix of output files

The command line to generate local AS events will be of the form:

### Example
```
python3 TGStools.py staAS -n K140.gtf,K510.gtf,SEC.gtf,SHEE.gtf,TE5.gtf -p TEST
```

### Output files

#### statistics of alternative splicing events
TEST_AS.txt
```
	A3	A5	AF	AL	MX	RI	SE
K140	2131	2224	4850	1050	330	2652	4608
K510	2638	2737	5562	1170	462	3130	5397
SEC	1689	1628	2911	431	261	2104	3267
SHEE	2849	2666	6876	954	389	3851	4870
TE5	1922    1870    4182    611     328     2553    4288
```

#### bar image
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/TEST_bar.1.png" width = "500" height = "400"  />
TEST_bar.png

#### barh_align image
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/TEST_barh_align.1.png" width = "400" height = "400"  />
TEST_barh_align.png

----------------------------
## calScoreD
----------------------------

calculate score_D of each gene

### Input files

#### ioi file
ioi file contains the transcript "events" of gene. produced by staAS, cantained in *_Event dietctory.

### Usage
```
python3 TGStools.py calScoreD -c <control> -t <treated>
```
List of options available:

- **-c**  | **--control**: name of control sample

- **-t**  | **--treated**: names of treated samples

- **-p**  | **--prefix**: prefix of output files

The command line to generate local AS events will be of the form:

### Example
```
python3 TGStools.py calScoreD -c K140 -t K510,SEC -p TEST
```

### Output files

TEST_score_D.txt
```
	K140_transcripts	K510_transcripts	K510_score	SEC_transcripts	SEC_score	
ENSG00000151466	ENST00000281142,ENST00000506368	ENST00000281142,ENST00000511426	0.667	ENST00000506368,ENST00000511426	0.667	1.334
ENSG00000128059	ENST00000264220	ENST00000264220	0	ENST00000264220	0	0
ENSG00000033178	ENST00000322244,ENST00000429659	ENST00000322244,ENST00000429659	0	ENST00000322244,ENST00000429659	0	0
ENSG00000145414	ENST00000274054	ENST00000274054	0	ENST00000274054	0	0
ENSG00000138767	ENST00000504123,ENST00000512485	ENST00000504123,ENST00000512485	0	ENST00000504804	1	1
```

To quantify the differential isoform usage between cells, we defined the score D of each gene as follows:

<img src="https://github.com/BioinformaticsSTU/SMRCanaToolkits/blob/master/CDZ/formula.png"  />

where gene j has isoform set a , and set b respectively in cell line X and Y ; c is the number of isoform intersection for set a and set b; d is the number of isoform union for set a and set b. Thus D sums up scores when comparing the control sample and treated samples.

----------------------------
## GOenrich
----------------------------

select top genes and make GO enrichment analysis

### Input files

#### score D file

### Usage
```
python3 TGStools.py GOenrich -i <input> -t <threshold> -n <number> -f <type> -p <prefix>
```
List of options available:

- **-i**  | **--input**: score D result

- **-t**  | **--threshold**: threshold for adjusted p-value

- **-n**  | **--number**: number of top genes for analysis

- **-f**  | **--type**: type of image, 'bar' and 'scatter' can be chosen

- **-p**  | **--prefix**: prefix of output files

### Example
```
python3 TGStools.py GOenrich -i TEST_score_D.txt -t 0.05 -n 500 -f all -p TEST
```

### Output files

#### TEST_GO_reports.txt
```
Gene_set	Term	Overlap	P-value	Adjusted P-value	Old P-value	Old Adjusted P-value	Z-score	Combined Score	Genes
GO_Biological_Process_2018	cellular response to DNA damage stimulus (GO:0006974)	26/330	9.13787228190808e-09	1.6576100319381273e-05	1.6286367911914452e-08	2.9543471392212826e-05	-1.3493239020889152	24.977116526080287	RPAIN;MUS81;BOD1L1;NUDT1;RECQL4;HERC2;RECQL5;CHEK2;MACROD1;NEK1;POLK;FNIP2;ZNF385A;VAV3;SHPRH;SLF1;DDX11;FANCA;RAD52;RNF168;RAD51B;PSME4;ATM;MMS19;CEP63;TP73
GO_Biological_Process_2018	regulation of striated muscle cell differentiation (GO:0051153)	3/9	0.0007078338095233947	0.18343007578220555	0.0015745180149574901	0.3173528532369876	-2.903223500208859	21.057954568005567	HDAC5;NRG1;HDAC9
GO_Biological_Process_2018	intraciliary retrograde transport (GO:0035721)	3/11	0.0013474224910073086	0.22220221806247806	0.0025281461092817154	0.36497538704896487	-2.67387109608216	17.67311619442282	DYNC2H1;ICK;DYNC2LI1
GO_Biological_Process_2018	double-strand break repair (GO:0006302)	11/142	0.00021322248264889278	0.09669639588127289	0.00028269216039523195	0.12820089473923768	-1.9459585449063224	16.4495269901104	RAD52;RECQL4;RNF168;RAD51B;HERC2;RECQL5;KDM4D;EME2;CHEK2;ATM;POLK
```

#### barh image
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/TEST_GO_barh.png"  width="600" height="450" />
TEST_GO_enrichment_barh.png

#### scatter image
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/TEST_GO_scatter.png" width="600" height="400" />
TEST_GO_enrichment_scatter.png





----------------------------
# License
----------------------------

TGStools is released under the MIT license.


