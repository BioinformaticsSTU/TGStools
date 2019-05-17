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

Third generation sequencing can de novo detect long reads of several thousand base pairs, thus provides a global view of the full length transcriptome. But due to less sequencing accurate rate, it often yields many spurious transcripts. It's important to prioritize the results by a visualization framework that automatically integrates rich annotation information with known genomic features. Therefore, we developed TGStools, a bioinformatics suit to facilitate routine tasks such as display transcripts of gene, characterizing the full-length transcripts and detecting the shifted types of alternative splicing in post transcriptome analysis.


----------------------------
# Installation
----------------------------

TGStools can run under linux and Python 3.5. In order to install TGStools successfully, it should be installed in ***conda*** environment and dependencies of TGStools are also designed to be installed through ***conda command***.
Conda should be installed and activate.
```
conda create -n python3 python=3.5
source activate python3
```
Because github has data upload limit, we have upload ***supplement data (gtfAnnotation.gtf)*** in OneDrive and you can download supplement data 
[here](https://stumail-my.sharepoint.cn/:u:/g/personal/d_z_chen_stu_edu_cn/ERG1zRvBkVFAn7mCyLeNvVoBGVbslQZQJIy-FUhF3LuGtA?e=bcgdYa
) and then put this file in ***source*** directory.
After conda environment has been activated, TGStools can be installed through command below.
By running setup.sh, TGStools and depedencies(PLEK, CNCI, libsvm, matplotlib, matplotlib_venn, pandas and gseapy) will be installed automatically.
```
source setup.sh
```


----------------------------
# Command and subcommand structure
----------------------------

TGStools works with a command/subcommand structure:

```
python3 TGStools.py subcommand options

```
where the subcommand can be one of these:

- **geneDisplay**    : create a macroscopic image showing transcripts of queried gene.
- **staDist**    : distances distribution of transcript-start-site (TSS ) in each full-length transcript to the closest epigenetic marks and CAGE tags.
- **INCP**    : an integration classification tool of CNCI and PLEK for identifying coding or non-coding transcripts (fasta file and gtf file).
- **extract_lncRNA_gtf**       : extract lncRNA information of GTF format based on the tanscript ID of the candidate lncRNA.
- **tiss_specific**       : extract tissue-specific lncRNA information of GTF format.
- **staAS**        : calculate the proportion of each alternative splicing event in different samples and create graphs.
- **calScoreD**     : calculate score_D(its formula can be seen below) of each gene.
- **GOenrich**     : select top ranked genes and conduct GO enrichment analysis.


----------------------------
## geneDisplay
----------------------------

By providing GTF files, gene id and optionally trans_quant files(quantity of transcripts) to get the macroscopic image which display the isoforms of queried gene. Users could provide multiple annotation data(bed format) such as known transcripts, epigenetic marks and CAGE tags in the same folder. These auxiliary annotations will be used to evaluate the isoforms detected from long reads sequencing.

For your convenience, we have upload epigenetic marks data in OneDrive. Uploaded data contains 3 types of histone marks from 15 tissues. You can download data [here](https://stumail-my.sharepoint.cn/:f:/g/personal/d_z_chen_stu_edu_cn/Enfeh4BW0vFJhi7cCsFaTUEBWimU5c5BH0ndF5SSw2TyLw?e=TLflsU).

### Input files

#### gtf file
An annotation file in GTF format is like:

```
chr14 Ensembl exon  73741918  73744001  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73749067  73749213  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";  
chr14 Ensembl exon  73750789  73751082  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73753818  73754022  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
```

#### path of auxiliary annotations files

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

- **-i**  | **--id**: gene id that you want to query

- **-q**  | **--quant**: quantity of transcript

- **-p**  | **--path**: path of auxiliary annotations.

### Example
```
python3 TGStools.py geneDisplay -g K510_3rd.gtf  -i ENSG00000035141 -q trans_quant.txt  -p histone
or
python3 TGStools.py geneDisplay -g K510_3rd.gtf  -i ENSG00000035141 -p histone
```

### Output files


#### isoforms comparison
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/img/ENSG00000035141.png" width = "550" height = "400"  />
ENSG00000035141.pdf


----------------------------
## staDist
----------------------------

Distances distribution of transcript-start-site(TSS ) in each full-length transcript to the closest epigenetic marks and CAGE tags.

**This step will cost much time and resource, we strongly recommend user to run this step at night while the computer is free.**

### Input files

#### gtf file

#### directory which contains auxiliary annotations files

### Usage
```
python3 TGStools.py staDist -g <gtf> -p <path> -f <flag>
```

- **-g**  | **--gtf**: gtf file

- **-p**  | **--path**: directory which contains auxiliary annotations files

- **-f**  | **--flag**: flag that tells programme the type of files in path

- **-r**  | **--prefix**: prefix of output files

### Example
```
python3 TGStools.py staDist -g K510_3rd.gtf  -p histone  -f histone -r TEST
or
python3 TGStools.py staDist -g K510_3rd.gtf  -p fantom5  -f fantom5 -r TEST
```

### Output files

<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/img/staDist.png" width = "550" height = "400"  />
staDist.pdf

----------------------------
## INCP
----------------------------

INCP(identify non-coding transcript from CNCI and PLEK) is an integration classification tool of CNCI and PLEK for identifying coding or non-coding transcripts(fasta file and gtf file).  
**This step will cost much time and resource, we strongly recommend user to run this step at night while the computer is free.**

### Input files

#### hg19.2bit and hg38.2bit
[hg19.2bit](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit) and [hg38.2bit](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit) are 2bit format of human genomes which can be downloaded on UCSC.

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

- **-p**  | **--parallel**: assign the running CPU numbers which should not larger than the CPU number your computer has

- **-g**  | **--gtf**: if your input file is gtf format please use this parameter

- **-r**  | **--reference**: if you use the -g this parameter must be assigned, within this parameter please assign the path of your reference genome. Some reference files which has been prepared could be download at [hg38](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit), [hg19](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit).

### Example
```
python3 TGStools.py INCP -i candidate.gtf -p 6 -g -r hg38.2bit
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

<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/img/venn.png" width = "550" height = "400"  />

venny image

----------------------------
## extract_lncRNA_gtf
----------------------------

Extract lncRNA information from GTF file based on the tanscript ID of the candidate lncRNA. This tool is only for result of gtf mode.

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

extract tissue-specific lncRNA information of GTF format.  
**This step will cost much time and resource, we strongly recommend user to run this step at night while the computer is free.**  
**If you want to run this step faster, input control and cancer sample and skip the parameter -t.**

### Input files

#### FILE 
files of the candidate lncRNA gtf format.

### Usage
```
python3 TGStools.py tiss_specific -i <file>[,<file>] [-t <tissue>] -r <reference> -o <out>
```

- **-i**  | **--input**: input files of the candidate lncRNA. If there are two input files which was splited by ',', the first file is of control sample, the other is of cancer sample. This step get cancer-specific transctipts from these file. If there is only one input file, background control tissue must be added. Related knowledge can refer to the  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3185964/

- **-t**  | **--tissue**: set tissue name if there is only one input file. If there are two input files, this parameter should be ignored. The background control tissue : 'adipose', 'adrenal', 'brain', 'breast', 'colon', 'heart', 'kidney', 'liver', 'lung', 'lymphNode', 'ovary', 'prostate', 'skeltalMuscle', 'whiteBloodCell', 'testes', 'thyroid', 'placenta', 'foreskin', 'hLF'.

- **-r**  | **--reference**: refgene name, 'hg38' and 'hg19' can be chosen.

- **-o**  | **--out**: output name extracted lncRNA-specific information of GTF format

### Example
```
python3 TGStools.py tiss_specific -i sample.gtf[,control.gtf] -t breast -r hg38 -o out.gtf
```

### Output files
GTF file which extract lncRNA-specific information

----------------------------
## staAS
----------------------------
Analyzing alternative events by SUPPA and calculate the proportion of each alternative splicing event in different samples and produce graphs.
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

#### Chi-square test result of samples in pair
TEST_chi2_result.txt
```
chi-square test result of K510 and K140:
chisq-statistic=14.5507, p-value=0.0241, df=6 expected_frep=[[ 2583.5706325   2687.58521866  5640.62432911  1202.66865258
    429.06016795  3132.35592306  5420.13507614]
 [ 2185.4293675   2273.41478134  4771.37567089  1017.33134742
    362.93983205  2649.64407694  4584.86492386]].
Proportion of alternative splicing in samples has difference.

chi-square test result of K510 and SHEEC:
chisq-statistic=127.8865, p-value=0.0000, df=6 expected_frep=[[ 2734.06990745  2758.08069009  5353.77266601  1011.61218438
    456.83673286  3307.1693773   5474.45844191]
 [ 1592.93009255  1606.91930991  3119.22733399   589.38781562
    266.16326714  1926.8306227   3189.54155809]].
Proportion of alternative splicing in samples has difference.
```

#### bar image
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/img/TEST_bar.1.png" width = "500" height = "400"  />
TEST_bar.png

#### barh_align image
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/img/TEST_barh_align.1.png" width = "400" height = "400"  />
TEST_barh_align.png

----------------------------
## calScoreD
----------------------------

calculate score_D of each gene.  
In order to quantify the differential isoform usage between cells, we defined the score D of each gene as follows:
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/img/formula.png"  />
where gene j has isoform set a , and set b respectively in cell line X and Y ; c is the number of isoform intersection for set a and set b; d is the number of isoform union for set a and set b. Thus D sums up scores when comparing the control sample and treated samples.
Genes with a higher D value are more diversely spliced.

### Input files

#### ioi file
ioi file which contains transcript "events" of each gene in *_Event dietctory and produced at last step staAS.

### Usage
```
python3 TGStools.py calScoreD -c <control> -t <treated> -p <prefix>
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

----------------------------
## GOenrich
----------------------------

select top ranked genes from score_D result and conduct GO enrichment analysis

### Input files

#### score D file

### Usage
```
python3 TGStools.py GOenrich -i <input> -t <threshold> -n <number> -f <type> -p <prefix>
```
List of options available:

- **-i**  | **--input**: score D result

- **-t**  | **--threshold**: threshold for adjusted p-value of GO enrichment analysis result

- **-n**  | **--number**: number of score_D top ranked genes 

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
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/img/TEST_GO_barh.png"  width="600" height="450" />
TEST_GO_enrichment_barh.png

#### scatter image
<img src="https://github.com/BioinformaticsSTU/TGStools/blob/master/img/TEST_GO_scatter.png" width="600" height="400" />
TEST_GO_enrichment_scatter.png


----------------------------
# License
----------------------------

TGStools is released under the MIT license.


