# ChIP-seq Practical 2019


1. [Introduction](#introduction)  
2. [Downloading ChIP-seq reads from NCBI](#download)
3. [Connect to the server and set up your environment](#setup)
4. [Quality control of the reads and statistics](#qc)
5. [Mapping the reads with Bowtie](#mapping)
6. [Estimating the number of duplicated reads](#dup)
7. [ChIP quality controls](#cqc)
8. [Visualizing the data in a genome browser](#visualize)
9. [Peak calling with MACS](#macs)
10. [Motif analysis](#motif)
11. [Peak annotation using R](#peakr)
12. [FAQ](#faq)
13. [References](#ref)



## Introduction <a name="introduction"></a>
### Goal
The aim is to :
  * have an understanding of the nature of ChIP-Seq data
  * perform a complete analysis workflow including quality check (QC), read mapping, visualization in a genome browser and peak-calling. Use command line and open source software for each step of the workflow and feel the complexity of the task
  * have an overview of possible downstream analyses


### Summary
This training gives an introduction to ChIP-seq data analysis, covering the processing steps starting from the reads to the peaks. Among all possible downstream analyses, the practical aspect will focus on motif analyses. This training does not cover all methods available today. It does not aim at bringing users to a professional NGS analyst level but provides enough information to allow biologists understand what DNA sequencing practically is and to communicate with NGS experts for more in-depth needs.

### Dataset description
For this training, we will use this dataset:
* a dataset produced by Myers et al [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/23818864) involved in the regulation of gene expression under anaerobic conditions in bacteria. We will focus on one factor: **FNR**. The advantage of this dataset is its small size, allowing real time execution of all steps of the dataset

## Downloading ChIP-seq reads from NCBI <a name="download"></a>
**Goal**: Identify the dataset corresponding to the studied article and retrieve the data (reads as FASTQ file) corresponding to one experiment and the control.  

### 1 - Obtaining an identifier for a chosen dataset
Within an article of interest, search for a sentence mentioning the deposition of the data in a database. Here, the following sentence can be found at the end of the Materials and Methods section:
*"All genome-wide data from this publication have been deposited in NCBI’s Gene Expression Omnibus (**GSE41195**)."*
We will thus use the **GSE41195** identifier to retrieve the dataset from the **NCBI GEO** (Gene Expression Omnibus) database.

NGS datasets are (usually) made freely accessible for other scientists, by depositing these datasets into specialized databanks. [Sequence Read Archive (SRA)](http://www.ncbi.nlm.nih.gov/sra) located in USA hosted by NCBI, and its European equivalent [European Nucleotide Archive (ENA)](http://www.ebi.ac.uk/ena) located in England hosted by EBI both contains **raw reads**.

Functional genomic datasets (transcriptomics, genome-wide binding such as ChIP-seq,...) are deposited in the databases [Gene Expression Omnibus (GEO)](http://www.ncbi.nlm.nih.gov/geo/) or its European equivalent [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/).

### 2 - Accessing GSE41195 from GEO
1.  The GEO database hosts processed data files and many details related to the experiments. The SRA (Sequence Read Archive) stores the actual raw sequence data.
2. Search in Google **GSE41195**. Click on the first link to directly access the correct page on the GEO database.
![alt text][geo]
3. This GEO entry is a mixture of expression analysis and chip-seq. At the bottom of the page, click on the subseries related to the chip-seq datasets. (this subseries has its own identifier: **GSE41187**).
![alt text][geo2]
4. From this page, we will focus on the experiment **FNR IP ChIP-seq Anaerobic A**. At the bottom of the page, click on the link "**GSM1010219** - FNR IP ChIP-seq Anaerobic A".
5. In the new page, go to the bottom to find the SRA identifier. This is the identifier of the raw dataset stored in the SRA database.  
![alt text][geo3]
6. Copy the identifier **SRX189773** (do not click on the link that would take you to the SRA database, see below why)

### 3 - Downloading FASTQ file from the ENA database
Although direct access to the SRA database at the NCBI is doable, SRA does not store the sequence FASTQ format. In practice, it's simpler (and quicker!!) to download datasets from the ENA database (European Nucleotide Archive) hosted by EBI (European Bioinformatics Institute) in UK. ENA encompasses the data from SRA.

1. Go to the [EBI](http://www.ebi.ac.uk/) website. Paste your SRA identifier (SRX189773) and click on the button "search".
![alt text][ebi4]
2. Click on the first result. On the next page, there is a link to the FASTQ file. For efficiency, this file has already been downloaded and will be available on the server. (SRR576933.fastq.gz)  
![alt text][ebi5]

**tip**: To download the control dataset, we should redo the same steps starting from the GEO web page specific to the chip-seq datasets (see step 2.4) and choose **anaerobic INPUT DNA**.  
The downloaded FASTQ file will also be available on the server (SRR576938.fastq.gz)

**At this point, you have two FASTQ files, one for the experiment, one for the control.**

## Connect to the server and set up your environment <a name="setup"></a>
### 1 - Sign in on the server
  * On MobaXterm
> Session : ssh  
> Host : core.cluster.france-bioinformatique.fr
> Specify username : ticked and filled in  
> Advanced SSH settings : X11-Forwarding  
  * On MacOS and Linux
```bash
ssh -XY <login>@core.cluster.france-bioinformatique.fr
```

### 2 - Set up your working environment
1. Go to your working directory
```bash
cd /shared/projects/ens_m2_2019/<login>
```
2. Load the conda virtual environment which contains all bioinformatics tools used to analyze ChIP-seq data.
```bash
module load conda
source activate eba2018_chipseq
```
3. Start an interactive session
```bash
sinteractive
```
4. Create a directory that will contain all results of the upcoming analyses.
```bash
mkdir cours_chipseq
```
5. Go to the newly created directory
```bash
cd cours_chipseq
```
6. <mark>Start an interactive session</mark>
```bash
sinteractive
```
7. Copy the directory containing data

```bash
cp -r /shared/projects/ens_m2_2019/data .
```

8. Your directory structure should be like this
 ```
/shared/projects/ens_m2_2019/<login>/cours_chipseq
│
└───data
```

## Quality control of the reads and statistics <a name="qc"></a>
**Goal**: Get some basic information on the data (read length, number of reads, global quality of dataset)  

### 1 - Getting the FASTQC report
Before you analyze the data, it is crucial to check the quality of the data. We will use the standard tool for checking the quality of data generated on the Illumina platform: [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

1. Create a directory named **01-QualityControl** in which to output results from fastqc
```bash
mkdir 01-QualityControl
```
2. Go to the directory you've just created
```bash
cd 01-QualityControl
```

Your directory structure should be like this
 ```
/shared/projects/ens_m2_2019/<login>/cours_chipseq
│
└───data
│   
└───01-QualityControl <- you should be in this folder
```
3. Check the help page of the program to see its usage and parameters.

```bash
fastqc --help
```
4. Launch the FASTQC program on the experiment file (SRR576933.fastq.gz)
  * -o: creates all output files in the specified output directory. '.' means current directory.
```bash
fastqc ../data/SRR576933.fastq.gz -o .
```  
5. Wait until the analysis is finished. Check the files output by the program.
```bash
ls
```
> SRR576933_fastqc.html  SRR576933_fastqc.zip

6. Download the HTML file SRR576933_fastqc.html on your local machine (either with ssh or the program you used to upload your data on the server). Using a bash command it would look like this.
```bash
### OPEN A NEW TERMINAL
## Create a directory where to put generated files on your computer
mkdir ~/Desktop/cours_chipseq

## Go to the location on your computer, where you want to put the data, for example:
cd ~/Desktop/cours_chipseq

## Download the file
scp <login>@core.cluster.france-bioinformatique.fr:/shared/projects/ens_m2_2019/<login>/cours_chipseq/01-QualityControl/SRR576933_fastqc.html .
# Enter your password
```
7. On your machine, open this file in Firefox.  
8. Launch the FASTQC program on the control file (SRR576938.fastq)

**Analyze the result of the FASTQC program:  
How many reads are present in the file ?  
What is the read length ?  
Is the overall quality good ?  
Are there any concerns raised by the report ? If so, can you tell where the problem might come from ?**  


### 2 - Organism length
Knowing your organism size is important to evaluate if your dataset has sufficient coverage to continue your analyses. For the human genome (3 Gb), we usually aim at least 10 Million reads.

1. Go to the [NCBI Genome](http://www.ncbi.nlm.nih.gov/genome) website, and search for the organism **Escherichia coli**
2. Scroll down up to the **Escherichia coli str. K-12 substr. MG1655** to access statistics on this genome.
![alt text][genome6]

**How long is the genome ?  
Do both FASTQ files contain enough reads for a proper analysis ?**

**At this point, you should be confident about the quality of the datasets, and whether it's worth following with analyzing the datasets.**

## Mapping the reads with Bowtie <a name="mapping"></a>
**Goal**: Obtain the coordinates of each read on the reference genome.  

### 1 - Choosing a mapping program
There are multiple programs to perform the mapping step. For reads produced by an Illumina machine for ChIP-seq, the currently "standard" programs are BWA and Bowtie (versions 1 and 2), and STAR is getting popular. We will use **Bowtie version 1.2.1.1** (Langmead B et al., Genome Biol, 2009) for this exercise, as this program remains effective for short reads (< 50bp).

### 2 - Bowtie
1. Try out bowtie
```bash
bowtie
```
This prints the help of the program. However, this is a bit difficult to read ! If you need to know more about the program, it's easier to directly check the manual on the [website](http://bowtie-bio.sourceforge.net/manual.shtml).

2. bowtie needs the reference genome to align each read on it. This genome needs to be in a specific format (=index) for bowtie to be able to use it. Several pre-built indexes are available for download on the bowtie webpage, but our genome is not there. You will need to make this index file.

3. Create a directory named **02-Mapping** in which to output mapping results
```bash
cd ..
mkdir 02-Mapping
```
4. Go to the directory you've just created
```bash
cd 02-Mapping
```


### 3 - Prepare the index file
1. To make the index file, you will need the complete genome, in FASTA format. It has already been downloaded to gain time (Escherichia_coli_K12.fasta in the course folder) (The genome was downloaded from the NCBI). Note that we will not work with the latest version (NC_000913.3) but the previous one (NC_000913.2), because the available tools for visualization have not been updated yet to the latest version. This will not affect our results.
2. Create a directory named **index** in which to output bowtie indexes
```bash
mkdir index
```
3. Go to the newly created directory
```bash
cd index
```
4. Try out bowtie-build
```bash
bowtie-build
```
5. Build the index for bowtie
```bash
## Creating genome index : provide the path to the genome file and the name to give to the index (Escherichia_coli_K12)
bowtie-build ../../data/Escherichia_coli_K12.fasta Escherichia_coli_K12
```
6. Go back to upper directory i.e 02-Mapping
```bash
cd ..
```

### 4 - Mapping the experiment
1. Create a directory named **IP** in which to put mapping results for IP
```bash
mkdir IP
```
2. Go to the newly created directory
```bash
cd IP
```
Your directory structure should be like this:
```
/shared/projects/ens_m2_2019/<login>/cours_chipseq
│
└───data
│   
└───01-QualityControl
│   
└───02-Mapping
|    └───index
|    └───IP
```


3. Let's see the parameters of bowtie before launching the mapping:
  * Escherichia_coli_K12 is the name of our genome index file
  * Number of mismatches for SOAP-like alignment policy (-v): to 2, which will allow two mismatches anywhere in the read, when aligning the read to the genome sequence.
  * Suppress all alignments for a read if more than n reportable alignments exist (-m): to 1, which will exclude the reads that do not map uniquely to the genome.
  * -q indicates the input file is in FASTQ format. SRR576933.fastq is the name of our FASTQ file.
  * -3 will trim x base from the end of the read. As our last position is of low quality, we'll trim 1 base.
  * -S will output the result in SAM format
  * 2> SRR576933.out will output some statistics about the mapping in the file SRR576933.out
```bash  
## Run alignment
bowtie ../index/Escherichia_coli_K12 ../../data/SRR576933.fastq.gz -v 2 -m 1 -3 1 -S 2> SRR576933.out > SRR576933.sam
```  
This should take few minutes as we work with a small genome. For the human genome, we would need either more time and more resources.

Bowtie output is a [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file. The SAM format correspond to large text files, that can be compressed ("zipped") into BAM format. The BAM files are usually sorted and indexed for fast access to the data it contains. The index of a given bam file is names .bam.bai or .bai file. Some tools require to have the index of the bam file to process it.

4. Sort the sam file and create a bam file using samtools
  * -b: output BAM
```bash
samtools sort SRR576933.sam | samtools view -b > SRR576933.bam
```

5. Create an index for the bam file
```bash
samtools index SRR576933.bam
```

6. Compress the .sam file (you could also delete the file)
```bash
gzip SRR576933.sam
```

**Analyze the result of the mapped reads:  
Open the file SRR576933.out (for example using the ` less ` command), which contains some statistics about the mapping. How many reads were mapped? How many multi-mapped reads were originally present in the sample? To quit less press 'q'**

### 5 - Mapping the control
1. Repeat the steps above (in 4 - Mapping the experiment) for the file SRR576938.fastq.gz in a directory named "**Control**" within the directory 02-Mapping.

**Analyze the result of the mapped reads:  
Open the file SRR576938.out. How many reads were mapped?**

## Estimating the number of duplicated reads <a name="dup"></a>
**Goal**: Duplicated reads i.e reads mapped at the same positions in the genome are present in ChIP-seq results. They can arise for several reasons including a biased amplification during the PCR step of the library prep, DNA fragment coming from repetitive elements of the genome, sequencing saturation or the same clusters read several times on the flowcell (i.e optical duplicates). As analyzing ChIP-Seq data consist at some point in detecting signal enrichment, we can not keep duplicated reads for subsequent analysis. So let's detect them using [Picard](http://broadinstitute.github.io/picard/)   

1. Go to the directory with alignment file of treatment (IP)
```bash
cd /shared/projects/ens_m2_2019/<login>/cours_chipseq/02-Mapping/IP
```
2. Run Picard markDuplicates to mark duplicated reads (= reads mapping at the exact same location on the genome)
  * CREATE_INDEX: Create .bai file for the result bam file with marked duplicate reads
  * INPUT: input file name to mark for duplicate reads
  * OUTPUT: output file name
  * METRICS: file with duplicates marking statistics
  * VALIDATION_STRINGENCY: Validation stringency for all SAM files read by picard.
```bash
srun picard MarkDuplicates \
CREATE_INDEX=true \
INPUT=SRR576933.bam \
OUTPUT=Marked_SRR576933.bam \
METRICS_FILE=metrics \
VALIDATION_STRINGENCY=STRICT
```

To determine the number of duplicated reads marked by Picard, we can run the `samtools flagstat` command:

```bash
samtools flagstat Marked_SRR576933.bam
```


Go back to working home directory (i.e /shared/projects/training/<login>/EBA2018_chipseq/)
```bash
## If you are in 02-Mapping/IP
cd ../..
```

## Visualizing the data in a genome browser <a name="visualize"></a>
**Goal**: View the data in their genomic context, to check whether the IP worked  

### 1 - Choosing a genome browser
There are several options for genome browsers, divided between the local browsers (need to install the program, eg. IGV) and the online web browsers (eg. UCSC genome browser, Ensembl). We often use both types, depending on the aim and the localization of the data.
If the data are on your computer, to prevent data transfer, it's easier to visualize the data locally (IGV). Note that if you're working on a non-model organism, the local viewer will be the only choice. If the aim is to share the results with your collaborators, view many tracks in the context of many existing annotations, then the online genome browsers are more suitable.

### 2 - Viewing the raw alignment data in IGV
1. Download the following files from the server onto your computer
  * data/Escherichia_coli_K12.fasta
  * data/Escherichia_coli_K_12_MG1655.annotation.fixed.bed
  * 02-Mapping/IP/SRR576933.bam
  * 02-Mapping/IP/SRR576933.bam.bai
  * 02-Mapping/Control/SRR576938.bam
  * 02-Mapping/Control/SRR576938.bam.bai
2. Open IGV on your computer
3. Load the genome
  * Genomes / Load Genome from File...
  * Select the fasta file Escherichia_coli_K12.fasta located into the data directory
4. Load an annotation file named Escherichia_coli_K_12_MG1655.annotation.fixed.bed into IGV
  * File / Load from File...
  * Select the annotation file. The positions of the genes are now loaded.
5. Load the two bam files (SRR576933.bam and SRR576938.bam) in IGV.

**Browse around in the genome. Do you see peaks?**  
**Browse into IGV. Go to the following genes: b1127, b1108**

However, looking at BAM file as such does not allow to directly compare the two samples as data are not normalized. Let's generate normalized data for visualization.

### 3 - Viewing scaled data
[bamCoverage](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html) from deepTools generates BigWigs from BAM files
1. Try it out
```bash
bamCoverage --help
```
2. Create a directory named **04-Visualization** to store bamCoverage outputs
```bash
mkdir 03-Visualization
```
3. Go to the newly created directory
```bash
cd 03-Visualization
```

Your directory structure should be like this:
```
/shared/projects/ens_m2_2019/<login>/cours_chipseq
│
└───data
│   
└───01-QualityControl
│   
└───02-Mapping
|    └───index
|    └───IP
|    └───Control
│   
└───03-Visualization <- you should be in this folder
```

4. Generate a scaled bigwig file on the IP with bamCoverage
  * --bam: BAM file to process
  * --outFileName: output file name
  * --outFileFormat: Output file type
  * --effectiveGenomeSize : size of the mappable genome
  * --normalizeUsing : different overall normalization methods; we will use RPGC method corresponding to 1x average coverage
  * --skipNonCoveredRegions: skip non-covered regions
  * --extendReads 200: Extend reads to fragment size
  * --ignoreDuplicates: reads that have the same orientation and start position will be considered only once
```bash
bamCoverage --bam ../02-Mapping/IP/Marked_SRR576933.bam \
--outFileName SRR576933_nodup.bw --outFileFormat bigwig --effectiveGenomeSize 4639675 \
--normalizeUsing RPGC --skipNonCoveredRegions --extendReads 200 --ignoreDuplicates
```
5. Do it for the control (be careful for the control you will need **5G** of memory to process the file)
6. Download the two bigwig files you have just generated
  * 04-Visualization/SRR576933_nodup.bw
  * 04-Visualization/SRR576938_nodup.bw
7. Load the two bigwig files in IGV
  * File / Load from File...
  * Select the two bigwig files.
8. Set the visualization of the two bigwig files to be autoscaled
  * Click right on the name of the tracks and select **Autoscale**

**Go back to the genes we looked at earlier: b1127, b1108. Look at the shape of the signal.**  
**Keep IGV opened.**

Go back to working home directory (i.e /shared/projects/ens_m2_2019/<login>/cours_chipseq)
```bash
## If you are in 03-Visualization
cd ..
```

## Peak calling with MACS <a name="macs"></a>
**Goal**: Define the peaks, i.e. the region with a high density of reads, where the studied factor was bound

### 1 - Choosing a peak-calling program
There are multiple programs to perform the peak-calling step. Some are more directed towards histone marks (broad peaks) while others are specific to narrow peaks (transcription factors). Here we will use MACS version 1.4.2 because it's known to produce generally good results, and it is well-maintained by the developer. A new version (MACS2) is available.

### 2 - Calling the peaks
1. Create a directory named **04-PeakCalling** to store annotatePeaks outputs
```bash
mkdir 04-PeakCalling
```
2. Go to the newly created directory
```bash
cd 04-PeakCalling
```
3. Try out MACS
```bash
macs
```
This prints the help of the program.

4. Let's see the parameters of MACS before launching the mapping:
  * ChIP-seq tag file (-t) is the name of our experiment (treatment) mapped read file SRR576933.bam
  * ChIP-seq control file (-c) is the name of our input (control) mapped read file SRR576938.bam
  * --format BAM indicates the input file are in BAM format. Other formats can be specified (SAM,BED...)
  * --gsize Effective genome size: this is the size of the genome considered "usable" for peak calling. This value is given by the MACS developers on their website. It is smaller than the complete genome because many regions are excluded (telomeres, highly repeated regions...). The default value is for human (2700000000.0), so we need to change it. As the value for E. coli is not provided, we will take the complete genome size 4639675.
  * --name provides a prefix for the output files. We set this to FNR_Anaerobic_A, but it could be any name.
  * --bw The bandwidth is the size of the fragment extracted from the gel electrophoresis or expected from sonication. By default, this value is 300bp. Usually, this value is indicated in the Methods section of publications. In the studied publication, a sentence mentions "400bp fragments (FNR libraries)". We thus set this value to 400.
  * --keep-dup specifies how MACS should treat the reads that are located at the exact same location (duplicates). The manual specifies that keeping only 1 representative of these "stacks" of reads is giving the best results. We doesn't mention it as 1 is the default value.
  <!-- * --bdg --single-profile will output a file in BEDGRAPH format to visualize the peak profiles in a genome browser. There will be one file for the treatment, and one for the control. -->
  * --diag is optional and increases the running time. It tests the saturation of the dataset, and gives an idea of how many peaks are found with subsets of the initial dataset.
  * &> MACS.out will output the verbosity (=information) in the file MACS.out
```bash
macs -t ../02-Mapping/IP/SRR576933.bam -c ../02-Mapping/Control/SRR576938.bam --format BAM  --gsize 4639675 \
--name "FNR_Anaerobic_A" --bw 400 --diag &> MACS.out
```
3. This should take a few minutes, mainly because of the --diag option. Without, the program runs faster.

### 3 - Analyzing the MACS results
**Look at the files that were created by MACS. Which files contains which information ?**  
**How many peaks were detected by MACS ?**

**At this point, you should have a BED file containing the peak coordinates.**

Go back to working home directory (i.e /shared/projects/ens_m2_2019/<login>/cours_chipseq)
```bash
## If you are in 04-PeakCalling
cd ..
```

### 4 - Visualize peaks into IGV

1. Download the BED file 04-PeakCalling/FNR_Anaerobic_A_peaks.bed to visualise in IGV.

**Go back again to the genes we looked at earlier: b1127, b1108. Do you see peaks?**


## References <a name="ref"></a>
This practical is jointly prepared with the [EBAI course](https://www.france-bioinformatique.fr/fr/evenements/EBAI2019)

[geo]: https://github.com/slegras/EBAI2017/blob/master/images/1_GEO.png "GEO"
[geo2]: https://github.com/slegras/EBAI2017/blob/master/images/2_GEO.png "GEO2"
[geo3]: https://github.com/slegras/EBAI2017/blob/master/images/3_GEO.png "GEO3"
[ebi4]: https://github.com/slegras/EBAI2017/blob/master/images/4_EBI.png "EBI"
[ebi5]: https://github.com/slegras/EBAI2017/blob/master/images/5_EBI.png "EBI"
[genome6]: https://github.com/slegras/EBAI2017/blob/master/images/6_Genomes.png "E. Coli K-12"
