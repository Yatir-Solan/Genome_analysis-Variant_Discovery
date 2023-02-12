# Variant Discovery - GitHub repository - 2nd edition

# Introduction

Genome analysis is a fundamental field within bioinformatics. Nevertheless, because my research subject as part of my computational biology MSc studies, was in other bioinformatics realms, I did not have the opportunity to dive deeply into it.

Therefore, I decided to gain knowledge in this field by conducting a small personal project - Whole Exome Sequencing (WES) analysis and **Variant discovery**. In this project, I decided to use **GATK.** GATK stands for Genome Analysis Toolkit, it was developed by the Broad institute and it is widely common for genome analysis and variant discovery. Besides developing individual tools, GATK is maintaining recommended research workflows, so scientists could follow them. These workflows, which are named GATK’s **best practices workflows ([link](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices)),** involve both tools from GATK, and other software, and they are comprehensively depicted and documented on GATK’s website. Here, I shall follow the workflow for ****Germline short variant discovery ([link](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)).** 

The main *.md file (the one you currently read), will follow the conducted research pipeline.  Here, code lines will appear beside a comprehensive text to elaborate on the tools, methods, and concepts behind it. The raw pipeline, which is mainly a bash script (research_pipeline.sh) can be found within the code section of the repository. In the final step of the pipeline, I use a Python script in a jupyter notebook, which is also provided (wannovar_analyse.ipynb). 

Hopefully, this repository, together with the documents in it, will help others who are new to this field. 

---

# Case study and objectives

- **Case study** - A newborn male baby, was born to two healthy parents. Several days after his birth he started having diarrhea, throwing up, and his skin turned somewhat yellow. The infant’s exome was obtained (SRR5439568**),** and was sent for genetic analysis.
- **Objectives** - To analyze the given exome and generate genetic diagnosis for his health condition.

---

# Technical subjects

- The research pipeline was run on a remote Linux machine. On it, all the programs and modules used were installed previously.
- Scientific writing - To keep things fast and flexible, I’m not going to follow strict academic writing rules. But, the tools I’m going to use will be clear and obvious, so that anyone interested in using them, could find them.

---

# Research Overview

1. **Obtaining the data**.
2. **Data Pre-processing** :
    1. Reads quality control.
    2. Mark duplicated reads.
    3. Base Quality Scores Recalibration (BQSR).
3. **Variant discovery**.
    1. Variant calling.
    2. Variant filtration.
4. **Variant** **filtration.**
5. **Variant annotation**.

### Germline short variant discovery

**Germline** mutations are mutations that can be found already in the germ cells, and therefore they are present within all other body cells and considered inheritable. In opposed to that, **Somatic** mutations are mutations that appear spontaneously within specific cells, and many times they are causative of cancer. 

---

# Research Pipeline Walkthrough

## 1. Obtaining the data

### Initialization

Loading the relevant modules.

```bash
module load AWScli/aws-cli-2.1.0
module load sratoolkit.3.0.1
module load bwa # version 0.7.17
module load samtools # version 1.9
module load bcftools/bcftools-1.6
module load gatk # version 4.1.7.0
```

The work will be arranged in 4 directories - reads, alignment, variants, and annotation. In addition, I shall arrange all paths as bash parameters, so all are commands will be tighter and more readable. 

```bash
# Creating directories :
mkdir /davidb/yatirsolan/NGS_variant_discovery # the analysis main directory.
cd /davidb/yatirsolan/NGS_variant_discovery
mkdir reads alignment variants annotationgatk MarkDuplicates \
    -I ${alignment_dir}/SRR5439568_srtd.bam \
    -O ${alignment_dir}/SRR5439568_srtd_mrkdpl.bam \
    -M ${alignment_dir}/SRR5439568_srtd_mrkdpl.mtrx \
    -ASO coordinate # Assume sort order - coordinate, because the BAM file is already sorted in coordinate order.

# Paths & Parameters initialization:
reads_dir="/davidb/yatirsolan/NGS_analysis/reads"
alignment_dir="/davidb/yatirsolan/NGS_analysis/alignment"
variants_dir="/davidb/yatirsolan/NGS_analysis/variants"
annotation_dir='/davidb/yatirsolan/NGS_variant_discovery/2nd_pipeline/annotation'
```

### Getting a reference genome

- **Reference genome** **download** - Because the assignment is a human genome analysis a human reference genome was needed. The latest available reference genome is the GRCh38. I decided to use the GATK version of it. It is a version that suits all GATK commands. The genome was downloaded from [here](https://ewels.github.io/AWS-iGenomes/) by the following bash lines.
    
    ```bash
    # Downloading GRCh38 reference genome
    aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ ./reference/Homo_sapiens/GATK/GRCh38/
    
    -->> ls -1 
    alignment
    reads
    reference
    variants
    annotation
    
    reference_genome_dir="/davidb/yatirsolan/NGS_analysis/reference/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta" # ref. genome directory path
    reference_genome=${reference_genome_dir}"/Homo_sapiens_assembly38.fasta" # reference genome path
    
    -->> ls ${reference_genome_dir} -1 
    Homo_sapiens_assembly38.fasta # fasta file
    Homo_sapiens_assembly38.fasta.fai # fasta index file 
    Homo_sapiens_assembly38.dict # sequence dictionary file 
    ```
    
- **Reference genome preparation** - In order to work with the reference genome, one should verify the presence of three components :
    - **Sequence** (*.faa) - the genome sequences as a fasta file.
    - **Sequence random access index files** (*.fai, *.dict) - the result of ***samtools faidx*** and ***samtools dict*** respectively. These files are crucial for tools such as GATK, but not mandatory for all softwares. The files are already placed beside the main fasta file. So, in the case of using GATK’s reference genome, there is no need to create them.
        - ***.fai** - In the file, each line describes the contigs within the fasta file - one contig per line. The values represented by columns are described in the code below.
            
            ```bash
            -->> samtools faidx ${reference_genome} # was not run eventually, because it was already present.
            -->> head -1 ${reference_genome_dir}/Homo_sapiens_assembly38.fasta.fai
            chr1    248956422       112     100     101  # contig, size (base pairs), location, bases per line, bytes per line 
            ```
            
        - ***.dict** - The file depicts the fasta file content.
            
            ```bash
            -->> samtools dict ${reference_genome} # was not run eventually.
            -->> head -2 ${reference_genome_dir}/Homo_sapiens_assembly38.dict
            @HD     VN:1.5  SO:unsorted
            @SQ     SN:chr1 LN:248956422    M5:6aef897c3d6ff0c78aff06ac189178dd     AS:38   UR:/seq/references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta        SP:Homo sapiens
            ```
            
    - **Alignment index files** (*.amb, *.ann, *.bwt, *.pac, *.sa) - Indexing the reference genome is a mandatory step for all alignment tools.
        - **Relocating index files** - Since I have used a GATK ‘ready to use’ product, the reference genome was already indexed. Nevertheless, in order for the pipeline to work, index files had to be relocated to be with the fasta files.
            
            ```bash
            # Examine the index files that cane in with the reference genome
            -->> ls /davidb/yatirsolan/NGS_analysis/reference/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/ -1
            Homo_sapiens_assembly38.fasta.64.alt **# what exactly is this file?**
            Homo_sapiens_assembly38.fasta.64.amb # BWA index file 
            Homo_sapiens_assembly38.fasta.64.ann # BWA index file 
            Homo_sapiens_assembly38.fasta.64.bwt # BWA index file 
            Homo_sapiens_assembly38.fasta.64.pac # BWA index file 
            Homo_sapiens_assembly38.fasta.64.sa # BWA index file 
            
            # copying the reference genome index files, so they will be in the same directory as the genome itself.
            -->> cp /davidb/yatirsolan/NGS_analysis/reference/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/* ${reference_genome_dir}
            ```
            
        - **Un-indexed genome** - In the case of using an un-indexed reference genome, we can use ***bwa index*** in order to index it. That can be done by the following line.
            
            ```bash
            -->> bwa index ${reference_genome} # was not run eventually.
            ```
            

### Getting the raw reads

- **Reads download** - The raw reads were obtained using ***fasterq-dmp*** of the sratoolkit, from the given [link](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5439568&display=metadata), using the accession SRR5439568.
    - *-O → output directory.*
    - *--split-files → creates two separated files.*
    - *accession → our sample (SRR5439568)*
    
    ```bash
    fasterq-dump -O ${reads_dir} --split-files SRR5439568
    
    -->> ls ${reads_dir} -1
    SRR5439568_1.fastq 
    SRR5439568_2.fastq
    
    raw_reads1=${reads_dir}"/SRR5439568_1.fastq" # row reads 1 fastq files path
    raw_reads2=${reads_dir}"/SRR5439568_2.fastq" # row reads 2 fastq files path
    ```
    

---

## 2. **Data Pre-processing**

### **Reads inspection**

- The reads are the result of Illumina 1.9 sequencing.
- Paired-End, 4,869,926 reads each.
- A rather small file (1.6 GB each) - which is consistent with an individual sample.
- Reads length - 101.

### Reads q**uality control (QC)**

To assess the quality of the reads I have used FastQC program using ******************fastqc****************** in the command line.

```bash
fastqc ${raw_reads1} ${raw_reads2}
```

**FastQC results -** 

Print screens of one of the files are given when needed. Through all the measurements - the QC results of both files look generally the same.

- **Base quality scores** - above 30 (green zone) for all reads (both files), and there is no significant quality decay toward the margins of the reads.
    
    ![Untitled](Variant%20Discovery%20-%20GitHub%20repository%20-%202nd%20editio%20c1c43825599c4bcaa255c7c6c364817d/Untitled.png)
    
- **Mean sequence quality scores** - The vast majority of the reads maintain mean quality scores that are above 30.
- **Per base sequence content** - Both files show a rather unbiased distribution. (to my knowledge the slight fluctuation within the beginning is fairly normal).
    
    ![Untitled](Variant%20Discovery%20-%20GitHub%20repository%20-%202nd%20editio%20c1c43825599c4bcaa255c7c6c364817d/Untitled%201.png)
    
- **Per sequence GC content** - **Warning**.
    
    ![Untitled](Variant%20Discovery%20-%20GitHub%20repository%20-%202nd%20editio%20c1c43825599c4bcaa255c7c6c364817d/Untitled%202.png)
    
- **Per base N content** - None for both.
- **Sequence length distribution** - 101 for all, in both.
- **Sequence duplication levels - Failure.**
    
    ![Untitled](Variant%20Discovery%20-%20GitHub%20repository%20-%202nd%20editio%20c1c43825599c4bcaa255c7c6c364817d/Untitled%203.png)
    
- **Overrepresented sequences - Warning.**
- **Adaptor content** - None for both.

**Quality Control conclusion** - The sequencing is in a good shape, but it needs to undergo duplicate reads removal (will be applied after the alignment).

### Reads Mapping (alignment)

We shall combine 3 commands into a single piped command line.

- **BWA mem** - Both BWA and STAR were considered as aligners. BWA was chosen because it is faster. Efficiency was also considered when choosing the aligner algorithm within BWA. Therefore, the ***bwa mem*** algorithm was chosen. Anyway, using ***bwa mem*** follows the recommended workflow by GATK.
    - *-R → Read group data. In case of a multiple samples alignment, a tag is needed to distinguish between them within the downstream analysis. In this example, even though we deals with a single individual, it is a good practice to maintain this convention. A proper read group tag header, is a tab-delimited string, that starts with ‘@’. It is common to reference the flowcell and the lane in the ID, the sequencing platform in the PL, and to give an informative SM (Sample Name). **Let’s come back here, and look [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups).***
    - *-t → Threads number.*
- **Coordinate Sort** - Sorting the alignment file using ***samtools sort*** (more about this step in the following).
- **Converting the SAM → BAM** - SAM (Sequence Alignment Map) is usually a large file, therefore it is more convenient to work with a binary file - BAM file (Binary sequence Alignment Map). We can do this with ***samtools view***.
    - *-h → (headers) Include the headers in the output.*
    - *-b → BAM output.*
    - *-o → output file name*

```bash
read_group="@RG\tID:SRR5439568\tPL:ILLUMINA\tSM:KID"
bwa mem -R ${read_group} -t 4 ${reference_genome} ${raw_reads1} ${raw_reads2} | \ # mapping the reads to the reference genome.
samtools sort | \                                                         # sorting the file in coordination with the reference genome.
samtools view -b -h -o ${alignment_dir}/SRR5439568_srtd.bam      # converting the final output to be in a BAM format.

-->> ls ${alignment_dir} -1
SRR5439568_srtd.bam # the output of the alignment
```

### Preparation for variant calling

- **The meaning of sorting -**
    
    For the sake of understanding this issue I shall create an unsorted alignment file. i.e., without piping it through a sorting command.
    
    ```bash
    -->> bwa mem -R ${read_group} -t 4 ${reference_genome} ${raw_reads1} ${raw_reads2} | \ # mapping the reads to the reference genome.
    		 samtools view -b -h -o ${alignment_dir}/SRR5439568.bam      # converting the final output to be in a BAM format.
    -->> ls ${alignment_dir} -1
    SRR5439568.bam
    SRR5439568_srtd.bam
    ```
    
    The reads in the unsorted BAM file (SRR5439568.bam), are located in the same order they were in the fastq files (the raw reads). We can observe this using that ***samtools view*** as shown in the code below. e.g., the pair of SRR5439568.1 comes first, then SRR5439568.2, SRR5439568.3, and so on. It makes more sense, to sort the reads in the order that they fall with respect to the reference genome. This procedure, which can be done by the ***samtools sort*** command, is also called coordinate sorting. In addition, software such as IGV requires a sorted alignment for them to work. Therefore, observing the sorted BAM file (SRR5439568_srtd.bam) shows a chromosomal-order file. i.e., the order of the reads is now coordinated with the genomic location on which they were mapped (chr 1).
    
    ```bash
    -->> samtools view ${alignment_dir}/SRR5439568.bam | head -6
    SRR5439568.1    99      chr18   9698878 60      101M    =       9699017 240     GGGGTTCTCTAGAGA...
    SRR5439568.1    147     chr18   9699017 60      101M    =       9698878 -240    CCAAAATCTGCAGGG...
    SRR5439568.2    83      chr12   102844371       60      101M    =       102844317       -155   ...
    SRR5439568.2    163     chr12   102844317       60      101M    =       102844371       155    ...
    SRR5439568.3    99      chr19   49125973        60      101M    =       49126098        226    ...
    SRR5439568.3    147     chr19   49126098        60      101M    =       49125973        -226   ...
    
    -->> samtools view ${alignment_dir}/SRR5439568_srtd.bam | head -6
    SRR5439568.2273567      163     chr1    10029   11      101M    =       10046   118     CCTAA...
    SRR5439568.2273567      83      chr1    10046   11      101M    =       10029   -118    CCCTA...
    SRR5439568.623404       163     chr1    10316   40      37M1I63M        =       10470   255  ...
    SRR5439568.623404       83      chr1    10470   40      101M    =       10316   -255    GCGGT...
    SRR5439568.698818       163     chr1    12249   0       101M    =       12354   206     CTCCT...
    SRR5439568.698818       83      chr1    12354   0       101M    =       12249   -206    CCGGC...
    ```
    
    In addition, the sorted file is smaller in size. That happens by virtue of the compression of the BAM file, which due to the sorting works more effectively.
    
    ```bash
    -->> du -h ${alignment_dir}/*
    876M    SRR5439568.bam
    519M    SRR5439568_srtd.bam
    -->> rm -f ${alignment_dir}/SRR5439568.bam # after using for demonstration we can delete the unsorted BAM file.
    ```
    
- **Alignment quality control** -
    
    ***samtools flagstat*** provides counts (statistics) of different SAM/BAM bitwise flag categories.  which are found in the SAM/BAM format.
    
    ```bash
    -->> samtools flagstat ${alignment_dir}/SRR5439568_srtd.bam
    9927085 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    187233 + 0 supplementary
    0 + 0 duplicates
    9779360 + 0 mapped (98.51% : N/A)
    9739852 + 0 paired in sequencing
    4869926 + 0 read1
    4869926 + 0 read2
    9268404 + 0 properly paired (95.16% : N/A)
    9566154 + 0 with itself and mate mapped
    25973 + 0 singletons (0.27% : N/A)
    170982 + 0 with mate mapped to a different chr
    146009 + 0 with mate mapped to a different chr (mapQ>=5)
    ```
    
- **Marking duplicates -**
    
    Duplicated reads produced by NGS are not rare. The duplication is being produced first, in the PCR amplification (PCR duplicates) step, on which duplicates are created on purpose, and therefore, by chance, some of the reads will be duplicated. In addition, duplicated reads can arise in the sequencing step, caused by optical errors of neighboring flowcells (Optical duplicates). Reads duplication can reduce the reliability of the analysis. As can be seen in the aforementioned statistics results, to this point, none of the reads are considered duplicates. Therefore we should locate them and tag them, so later they will be ignored. Here, we shall use the GATK **MarkDuplicates** command. A command that adds tags to the duplicated reads. Later on in the pipeline, most downstream tools ignore the tagged duplicates. In fact, **MarkDuplicates** is GATK’s wrapper for the **Picard** program algorithm (with the same name), and it could be run directly using Picard.
    
    - *-I → (input) → the raw BAM file.*
    - *-O → (output) → the marked BAM file.*
    - *-M → (metrics file) → a file to write duplication metrics to.*
    - -*ASO → (Assume sort order) → coordinate, because the BAM file is already sorted in coordinant order.*
    
    After we have marked duplicate reads, ***samtools flagstat*** command counts them as well (6,403,892). 
    
    ```bash
    gatk MarkDuplicates \
        -I ${alignment_dir}/SRR5439568_srtd.bam \
        -O ${alignment_dir}/SRR5439568_srtd_mrkdpl.bam \
        -M ${alignment_dir}/SRR5439568_srtd_mrkdpl.mtrx \
        -ASO coordinate
    
    -->> samtools flagstat ${alignment_dir}/SRR5439568_srtd_fltrd_mrkdpl.bam
    9927085 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    187233 + 0 supplementary
    6403892 + 0 duplicates
    9779360 + 0 mapped (98.51% : N/A)
    9739852 + 0 paired in sequencing
    4869926 + 0 read1
    4869926 + 0 read2
    9268404 + 0 properly paired (95.16% : N/A)
    9566154 + 0 with itself and mate mapped
    25973 + 0 singletons (0.27% : N/A)
    170982 + 0 with mate mapped to a different chr
    146009 + 0 with mate mapped to a different chr (mapQ>=5)
    ```
    
- **Base quality recalibration (BQSR) -**
    
    Quality scores of the bases (which were given to us by the sequencer) have a significant impact on a variant discovery analysis. Basically, the base scores are taken into account within variant calling algorithms, so that bases with lower scores will be treated with less confidence as opposed to others with higher quality scores. But, it has been shown, that quality scores can be significantly biased, mostly as a result of nucleotide context. Therefore, this step is aimed to locate systematically biased quality scores and to recalibrate them. BQSR procedure involves a ML process in which a set of known variants is given as input, then mismatches of the sample to the reference which are present in the known set are skipped so that the quality scores of unknown mismatches could be re-inferred. The procedure involves two GATK functions, ***BaseRecalibrator*** which creates the recalibrated scores and outputs these into a recalibration table, and ***ApplyBQSR*** applies the recalibrated scores to the BAM file. 
    
    - ***BaseRecalibrator -***
        - *-R → (reference genome)*
        - *-I → (input) → The BAM file the model should use for learning.*
        - *--known-sites → Database of known polymorphic sites used to exclude regions around known polymorphisms from analysis.*
        - *O → (output) → The output recalibration table file to create.*
    - ***ApplyBQSR -***
        - *-R → reference genome*
        - *-I → (input) → The BAM file that should be recalibrated.*
        - *-bqsr → (recalibration table) → the output of **BaseRecalibrator.***
        - *O → (output) → The calibrated BAM file.*
    
    Interestingly, ***ApplyBQSR*** also provides us with a new *.bai index file (will be elaborated soon).
    
    ```bash
    known_sites_dir='/davidb/yatirsolan/NGS_variant_discovery/2nd_pipeline/known_sites'
    
    # Downloading the files crucial for the known set in the BQSR stage.
    wget -P ${known_sites_dir} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
    wget -P ${known_sites_dir} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
    
    known_sites=${known_sites_dir}/Homo_sapiens_assembly38.dbsnp138.vcf
    
    # BQSR - Base quality scores adjustments
     # 1) creating the model
    gatk BaseRecalibrator \
    	  -R ${reference_genome} \
    	  -I ${alignment_dir}/SRR5439568_srtd_mrkdpl.bam \
    	  --known-sites ${known_sites} \
    	  -O ${alignment_dir}/recalibration.table
    
    # 2) applying the model
    gatk ApplyBQSR \
        -R ${reference_genome} \
        -I ${alignment_dir}/SRR5439568_srtd_mrkdpl.bam \
        *-bqsr* ${alignment_dir}/recalibration.table \
        -O ${alignment_dir}/SRR5439568_srtd_mrkdpl_bqsr.bam
    
    -- >> ls ${alignment_dir}/*_bqsr* -1
    SRR5439568_srtd_mrkdpl_bqsr.bam
    SRR5439568_srtd_mrkdpl_bqsr.bai # gatk ApplyBQSR creates an index file.
    ```
    
- **Indexing the BAM file** -
    
    A BAM index file is a file that allows efficient access to specific regions within the alignment. This functionality is also crucial to visualization tools, such as IGV, which for it *.bai file (the index file), is mandatory. In our pipeline, we have benefited from GATK’s ***ApplyBQSR*** command, which by default also creates a *.bai file (as can be seen above). When using equivalent tools which do not produce an index file, we can utilize ***samtools index*** to do so***.*** Besides visualization tools, whenever an index file is present, we can manually address reads within the alignment using ***samtools view***. In the code below, I show how can ***samtools index*** can be used, even though it was not in use within our pipeline. Additionally, I show how can we use the indexed file with ***samtools view**,* and we lack this functionality when a BAM file is not indexed.
    
    ```bash
    -->> samtools index ${alignment_dir}/SRR5439568_srtd_mrkdpl_bqsr.bam # was not run eventually. 
    -->> ls  ${alignment_dir} -1
    SRR5439568_srtd.bam
    SRR5439568_srtd_mrkdpl.bam
    SRR5439568_srtd_mrkdpl_bqsr.bai # gatk ApplyBQSR creates an index file.
    SRR5439568_srtd_mrkdpl_bqsr.bam
    SRR5439568_srtd_mrkdpl.mtrx
    
    -->> samtools view ${alignment_dir}/SRR5439568_srtd_mrkdpl_bqsr.bam chr1:1000-10500 # calling reads that were aligned in this range
    SRR5439568.2273567      163     chr1    10029   11      101M    =       10046   118     CCTAAC...
    SRR5439568.2273567      83      chr1    10046   11      101M    =       10029   -118    CCCTAA...
    SRR5439568.623404       163     chr1    10316   40      37M1I63M        =       10470   255   ...
    SRR5439568.623404       83      chr1    10470   40      101M    =       10316   -255    GCGGTA...
    
    -->> samtools view ${alignment_dir}/SRR5439568_srtd_mrkdpl.bam chr1:1000-10500 # trying to do the same call with an unindexed file
    [main_samview] random alignment retrieval only works for indexed BAM or CRAM files.
    ```
    
- **Analysis ready BAM file -**
    
    The data preprocessing step objective is to provide an alignment file, which is filtered and clean, on which the next step (variant discovery), could conduct further analysis. Hence, marking the duplicates was the final step in the data preprocessing stage. 
    
    ```bash
    analaysis_ready_bam=${alignment_dir}/SRR5439568_srtd_mrkdpl_bqsr.bam
    ```
    

---

## 3. Variants Discovery

### Variants Calling

This step is performed by GATK’s ***HaplotypeCaller*** function. This function is the most essential in our pipeline, and its details are covered in more depth in the [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035531412-HaplotypeCaller-in-a-nutshell). Here, we shall brief the principles. The algorithm’s main objective is to locate systematic mismatches of the aligned reads with respect to the reference genome. Or, in simple words - to identify variants. First, the program utilizes a sliding window to look for locations with significant evidence for variation, which are defined as **active regions**. Then, each active region is **reassembled** using De Bruijn-like to produce possible haplotypes, and re-align them to the reference. Afterward, a **haplotype likelihood** for each genomic position within the active regions is calculated. And eventually, a **genotype likelihood** for each variant position is calculated.

- *-R → (reference genome)*
- *-I → (input) → The BAM file after all the preparation which we did on it.*
- *O → (output) → A VCF file, which holds all the variants found by the algorithm.*
    
    ```bash
    gatk HaplotypeCaller \
    	  -R ${reference_genome} \
    	  -I ${analaysis_ready_bam} \
    	  -O ${variants_dir}/SRR5439568_hptpcl.vcf
    ```
    

---

## 4. Variants Filtration

In this stage, we aim to omit false positives variants, which have emerged from the variant calling step. There are two filtration approaches - Hard filtering and VQSR (Variant Quality Score Recalibration). **VQSR** is the more sophisticated approach. It uses large datasets of known variants to teach ML models the characteristics of true and false variants so that deciphering between variants on the examined data could be made. Nevertheless, VQSR requires data that is present in humans mostly, and for some model organisms. But it does not in many others. On the other hand, **Hard filtering** is a relatively crude and simple procedure, because it basically utilizes thresholds of different measurements to filter the variants. Therefore, Hard filtering is much more robust and fits all organisms. Because of that, even though it is human research, we shall stick with the **Hard filtering** approach. 

- **Dividing the variants -**
    
    GATK holds specific recommendations regarding the proper thresholds for variant filtration. These thresholds differentiate between SNPs and Indels. Therefore, we need to divide the variants, so that we can manipulate each of the variant types with the specification they require. We will utilize GATK’s ***SelectVariants*** function ******to separate the variants into two distinct files.
    
    - *-V → (variants) → A VCF file with the variants.*
    - *-select-type → the type of variants that should be included.*
    - *O → (output) → A VCF file, on which the selected variants should be written.*
        
        ```bash
        # Deviding the variants to SNOs and Indels
         # 1) SNPs
        gatk SelectVariants \
            -V ${variants_dir}/SRR5439568_hptpcl.vcf \
            -select-type SNP \
            -O ${variants_dir}/SRR5439568_hptpcl_snps.vcf 
        
         # 2) Indels
        gatk SelectVariants \
             -V ${variants_dir}/SRR5439568_hptpcl.vcf \
             -select-type INDEL \
             -O ${variants_dir}/SRR5439568_hptpcl_indels.vcf
        ```
        
- **Adding Filtration tags to the variants** -
    
    Now, with the use of GATK’s ***VariantFiltration*** function, variants in each of the files will be tagged with respect to the thresholds they will be run through. Filters should be applied both on site-level properties and on sample-level properties. Regarding site-level filtration, we can specify the corresponding thresholds measurements, so that tags in the FILTER column (within the VCF file), will appear for variants that didn’t pass the filtration, and the tag PASS will appear for the ones that did. The logic is the same regarding sample-level filtration, but the filtration tags will appear in the corresponding column (SAMPLE), and the FORMAT column will be updated for it as well. Site-level measurements on which thresholds will be applied are: *QD*, *FS*, *MQ*, *SOR*, *MQRankSum*, and *ReadPosRankSum*. For the sample-level: *DP, GQ*. To keep things tight, I will not go over each of the thresholds meaning, and the reasons for it to be included. For further elaboration regarding the parameters - [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).  
    
    - *-V → (variants) → Single-type variants VCF file.*
    - *-filter-name → A string to appear in the FILTER column of non-passing threshold variants.*
    - *-filter → The threshold itself.*
    - *-genotype-filter-name → A string to appear in the Sample column of non-passing threshold variants.*
    - *-genotype-filter-expression → The actual sample-level filtration.*
    - *O → (output) → The same VCF file enriched with filtration tags at the FILTER column.*
    
    ```bash
    # Adding filtration tags
     # 1) SNPs
    gatk VariantFiltration \
        -V ${variants_dir}/SRR5439568_hptpcl_snps.vcf \
        -O ${variants_dir}/SRR5439568_hptpcl_snps_tag.vcf \
        -filter-name 'QD_flt' -filter 'QD < 2.0' \
        -filter-name 'FS_flt' -filter 'FS > 60.0' \
        -filter-name 'SOR_flt' -filter 'SOR > 4.0' \
        -filter-name 'MQ_flt' -filter 'MQ < 40.0' \
        -filter-name 'MQRankSum_flt' -filter 'MQRankSum < -12.5' \
        -filter-name 'ReadPosRankSum_flt' -filter 'ReadPosRankSum < -8.0' \
        -genotype-filter-expression 'DP < 10' \
        -genotype-filter-name 'DP_flt' \
        -genotype-filter-expression 'GQ < 10' \
        -genotype-filter-name 'GQ_flt' 
    
     # 2) Indels
    gatk VariantFiltration \
        -V ${variants_dir}/SRR5439568_hptpcl_indels.vcf \
        -O ${variants_dir}/SRR5439568_hptpcl_indels_tag.vcf \
        -filter-name 'QD_flt' -filter 'QD < 2.0' \
        -filter-name 'FS_flt' -filter 'FS > 200.0' \
        -filter-name 'SOR_flt' -filter 'SOR > 10.0' \
        -genotype-filter-expression 'DP < 10' \
        -genotype-filter-name 'DP_flt' \
        -genotype-filter-expression 'GQ < 10' \
        -genotype-filter-name 'GQ_flt'
    ```
    
    To see the results, for each of the SNPs files (tagged and untagged) I grepped the first 3 variant records, so we could see the filtration tags as they appear in the FILTER column. Note the dots (’.’) at the untagged file, which represent records that were not targeted with any filtration procedure. 
    
    ```bash
    -->> grep -v ^## ${variants_dir}/SRR5439568_hptpcl_snps.vcf | head -4        # The first 3 records of the original SNPs VCF 
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  KID
    chr1    16534   .       C       T       30.64   .       AC=1;AF=0.500;AN=2;BaseQRankSum=0.674;DP=5;ExcessHet=3.0103;FS=0.000;MLEAC...
    chr1    183629  .       G       A       32.64   .       AC=1;AF=0.500;AN=2;BaseQRankSum=-0.967;DP=3;ExcessHet=3.0103;FS=0.000;MLEA...
    chr1    187102  .       C       G       34.64   .       AC=1;AF=0.500;AN=2;BaseQRankSum=-0.402;DP=32;ExcessHet=3.0103;FS=1.845;MLE...
    
    -->> grep -v ^## ${variants_dir}/SRR5439568_hptpcl_snps_tag.vcf | head -4    # The first 3 records of the tagged SNPs VCF
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  KID
    chr1    16534   .       C       T       30.64   MQ_flt  AC=1;AF=0.500;AN=2;BaseQRankSum=0.674;DP=5;ExcessHet=3.0103;FS=0.000;MLEAC...
    chr1    183629  .       G       A       32.64   PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=-0.967;DP=3;ExcessHet=3.0103;FS=0.000;MLEA...
    chr1    187102  .       C       G       34.64   QD_flt  AC=1;AF=0.500;AN=2;BaseQRankSum=-0.402;DP=32;ExcessHet=3.0103;FS=1.845;MLE...
    ```
    
    Nevertheless, the FILTER columns hold tags for site-level filtration only. As was demonstrated previously, we also filtered variants concerning the sample-level properties. To see the tags regarding those, we should examine the **KID** column, which is the only sample presence in our case. In the code box below I concluded five specific variant records, all of which passed the site-level filters, therefore they all are marked with PASS in the FILTER column. Nevertheless, concerning sample-level, only one of them passed the filters. That can be seen in the KID column, for which the filtration names (DP_flt, GQ_flt) we applied are found in the ones they failed it. 
    
    ```bash
    -->> grep -v '^##' ${variants_dir}/SRR5439568_hptpcl_snps_tag.vcf | awk 'NR==1||NR>=57&&NR<=61' | cut -f 7,9,10
    FILTER     FORMAT                  KID
    PASS       GT:AD:DP:FT:GQ:PL       0/1:1,2:3:DP_flt:22:65,0,22
    PASS       GT:AD:DP:FT:GQ:PL       1/1:0,4:4:DP_flt:12:143,12,0
    PASS       GT:AD:DP:GQ:PL          1/1:0,12:12:36:377,36,0 # a variant record that passed both site, and sample levels filtration.
    PASS       GT:AD:DP:FT:GQ:PL       1/1:0,2:2:DP_flt;GQ_flt:6:49,6,0
    PASS       GT:AD:DP:FT:GQ:PL       1/1:0,2:2:DP_flt;GQ_flt:6:49,6,0
    ```
    
- **Excluding variants that didn’t pass filtration -**
    
    To exclude records that didn’t pass site-level filtration we use GATK’s ***SelectVariants*** function.
    
    - *-V → (variants) → A VCF file with the variants.*
    - *--exclude-filtered → Don't include filtered sites {true, false}*
    - *O → (output) → A VCF file, without the site-level failed filtration variant records.*
    
    Afterward, I use bash scripting to manually exclude the records that failed the sample-level filtration. Specifically, we exclude (using grep -v) sample-level filtration names we specified.
    
    ```bash
    # Excluding the filtered variants
     # 1.1) SNPs - Exclude variants that failed to pass site-level filtration.
    gatk SelectVariants \
        -V ${variants_dir}/SRR5439568_hptpcl_snps_tag.vcf \
        -O ${variants_dir}/SRR5439568_hptpcl_snps_fltrd.vcf \
        --exclude-filtered true 
     # 1.2) SNPs - Exclude variants that failed to pass sample-level filtration.
    echo "$(grep -v -E 'DP_flt|GQ_flt' ${variants_dir}/SRR5439568_hptpcl_snps_fltrd.vcf)" > \
    			${variants_dir}/SRR5439568_hptpcl_snps_fltrd.vcf
    
    # --------------- for the indels, the same script except for file names --------------- #
    
     # 2.1) Indels - Exclude variants that failed to pass site-level filtration.
    gatk SelectVariants \
        -V ${variants_dir}/SRR5439568_hptpcl_indels_tag.vcf \
        -O ${variants_dir}/SRR5439568_hptpcl_indels_fltrd.vcf \
        --exclude-filtered true
     # 2.2) Indels - Exclude variants that failed to pass sample-level filtration.
    echo "$(grep -v -E 'DP_flt|GQ_flt' ${variants_dir}/SRR5439568_hptpcl_indels_fltrd.vcf)" > \
    			${variants_dir}/SRR5439568_hptpcl_indels_fltrd.vcf
    ```
    

The previous step concluded all the preparations the variants should apply with before being functionally annotated. Therefore, the final VCF files that were produced are considered ready for functional annotation VCF files.

```bash
annotation_ready_SNPs=${variants_dir}/SRR5439568_hptpcl_snps_fltrd.vcf
annotation_ready_indels=${variants_dir}/SRR5439568_hptpcl_indels_fltrd.vcf
```

---

## 5. Variants Annotation (Post VCF stage)

This is the final major step in our analysis. We have obtained a VCF file with filtered variants records. Now, we need to understand which one of them causing the clinical problem. The idea is pretty simple, we benefit from the comprehensive available data that treasure within it details regarding known variants and their medical significance, and genomic positions in general. Therefore, the procedure we aim to do, is to compare the variants we identified with this data, and if the variants we found are known, i.e., have been reported and studied before, we shall extract the available data, and annotate our own with it. In GATK there is a tool oriented exactly for the aforementioned, and its name is ***Funcotator*** (FUNCtional annOTATOR).

### Obtating data to annotate with

First, we need to provide GATK with the data sources (the so-called, known data). To improve the simplicity of ***Funcotator*** usage, GATK built an additional tool named ***FuncotatorDataSourceDownloader**.* Running this tool, as can be seen in the code box below, will download a suitable data source for us to use. 

- *-gremline → Because our analysis involves germline survey, we specify it. One can choose out of the options {true/false}.*
- -*extract-after-download →*  *true. so that extraction will happen immediately after the download. One can choose out of the options {true/false}.*
- *-validate-integrity →  true, so validation will be applied. One can choose out of the options {true/false}.*

```bash
gatk FuncotatorDataSourceDownloader -germline true -validate-integrity true -extract-after-download true 

-->> ls -1 
alignment
annotation
funcotator_dataSources.v1.6.20190124g # the data source folder which was already extracted
funcotator_dataSources.v1.6.20190124g.tar.gz
known_sites
variants

data_sources_dir='/davidb/yatirsolan/NGS_variant_discovery/2nd_pipeline/funcotator_dataSources.v1.6.20190124g'
```

### Functional annotation

Now, we move to the actual part of the annotation. In the code box below the usage of ***Funcotator*** can be seen. ******

- *-R → reference genome.*
- *--ref-version → The reference genome version we work with, which in our case is the hg38 (equivalent to GRCh38).*
- *--data-sources-path → The path to the directory we have just downloaded by **FuncotatorDataSourceDownloader***.
- *-V → (variants) → The VCF file targeted for functional annotation.*
- *O → (output) → The annotations-enriched VCF file.*
- *--output-file-format → VCF (we should choose from {VCF, MAF, SEG}.*

```bash
# Adding a functional annotation to the INFO column of the VCF.
 # 1) SNPs
gatk Funcotator \
    -R ${reference_genome} \
    --ref-version hg38 \
    --data-sources-path ${data_sources_dir} \
    -V ${annotation_ready_SNPs} \
    -O ${annotation_dir}/SRR5439568_hptpcl_snps_fltrd_anttd.vcf \
    --output-file-format VCF

# --------------- for the indels, the same script except for file names --------------- #

 # 2) Indels
gatk Funcotator \
		-R ${reference_genome} \
    --ref-version hg38 \
    --data-sources-path ${data_sources_dir} \
    -V ${annotation_ready_indels} \
    -O ${annotation_dir}/SRR5439568_hptpcl_indels_fltrd_anttd.vcf \
    --output-file-format VCF
```

The output of GATK’s ***Funcotator*** function is the same VCF but enriched with the annotations within it. The annotations are designated to appear in the INFO column, and they are described in the header sections as FUNCOTATION. In the code box below I’ve selected the first 4 variant records and asked to see the INFO column only. As can be seen, after the ‘FUNCOTATION=’ string in brackets ‘[’ / ‘]’ there are several values delimited by vertical bars ‘|’ pipes. These values are the functional annotations, that were enriched by ***Funcotator.*** Each corresponds with the variant record they are found in.

```bash
-->> grep -v '^##' SRR5439568_Hplocall_snps_fltrd_anttd.vcf | cut -f 8 | head
INFO
AC=2;AF=1.00;AN=2;DP=11;ExcessHet=3.0103;FS=0.000;FUNCOTATION=[AGRN|hg38|chr1|1053552|1053552|INTRON||SNP|G|G|C|g.chr1:1053552G>C|...
AC=2;AF=1.00;AN=2;DP=12;ExcessHet=3.0103;FS=0.000;FUNCOTATION=[MXRA8|hg38|chr1|1361438|1361438|FIVE_PRIME_UTR||SNP|C|C|G|g.chr1:13...
AC=2;AF=1.00;AN=2;DP=33;ExcessHet=3.0103;FS=0.000;FUNCOTATION=[MTHFR|hg38|chr1|11791061|11791061|INTRON||SNP|A|A|G|g.chr1:11791061...
AC=2;AF=1.00;AN=2;DP=102;ExcessHet=3.0103;FS=0.000;FUNCOTATION=[MTHFR|hg38|chr1|11794400|11794400|SILENT||SNP|G|G|A|g.chr1:1179440...
```

To interpret the meanings of each of the values, we sholud refer with the correspondent INFO raw in the header section that starts the VCF. In the code box below I grepped only the relevant header. At this raw, following the string “Funcotation fields are: “ the columns describing the annotations values are found.

```bash
-->> grep '^##' ${annotation_dir}/SRR5439568_hptpcl_snps_fltrd_anttd.vcf | grep 'FUNCOTATION'
##INFO=<ID=FUNCOTATION,Number=A,Type=String,Description="Functional annotation from the Funcotator tool.  Funcotation fields are: Gencode_27_hugoSymbol|Gencode_27_ncbiBuild|Gencode_27_chromosome|Gencode_27_start|Gencode_27_end|Gencode_27_variantClassification|Gencode_27_secondaryVariantClassification|Gencode_27_variantType|Gencode_27_refAllele|Gencode_27_tumorSeqAllele1|Gencode_27_tumorSeqAllele2|Gencode_27_genomeChange|Gencode_27_annotationTranscript|Gencode_27_transcriptStrand|Gencode_27_transcriptExon|Gencode_27_transcriptPos|Gencode_27_cDnaChange|Gencode_27_codonChange|Gencode_27_proteinChange|Gencode_27_gcContent|Gencode_27_referenceContext|Gencode_27_otherTranscripts|ACMGLMMLof_LOF_Mechanism|ACMGLMMLof_Mode_of_Inheritance|ACMGLMMLof_Notes|ACMG_recommendation_Disease_Name|ClinVar_VCF_AF_ESP|ClinVar_VCF_AF_EXAC|ClinVar_VCF_AF_TGP|ClinVar_VCF_ALLELEID|ClinVar_VCF_CLNDISDB|ClinVar_VCF_CLNDISDBINCL|ClinVar_VCF_CLNDN|ClinVar_VCF_CLNDNINCL|ClinVar_VCF_CLNHGVS|ClinVar_VCF_CLNREVSTAT|ClinVar_VCF_CLNSIG|ClinVar_VCF_CLNSIGCONF|ClinVar_VCF_CLNSIGINCL|ClinVar_VCF_CLNVC|ClinVar_VCF_CLNVCSO|ClinVar_VCF_CLNVI|ClinVar_VCF_DBVARID|ClinVar_VCF_GENEINFO|ClinVar_VCF_MC|ClinVar_VCF_ORIGIN|ClinVar_VCF_RS|ClinVar_VCF_SSR|ClinVar_VCF_ID|ClinVar_VCF_FILTER|LMMKnown_LMM_FLAGGED|LMMKnown_ID|LMMKnown_FILTER">
```

### Extracting the details into a readable annotated variants table

To understand the annotations, it is much more convenient to have them as a table. This step is a bit hacky and includes some more advanced bash commands. Nevertheless, all the bash script can be easily done in any other suitable tool, such as Python or R for instance.

- **Create a file with column names** - From the FUNCOTATION header, we import the names of the columns, to be the first raw of a new tab-delimited file (tsv), which is the variants annotations table file.

```bash
 # 1) SNPs 
  # 1.1) Extracting the column's names from the header section of the VCF into an a new empty file. 
funcotation_header=$(grep '^##' ${annotation_dir}/SRR5439568_hptpcl_snps_fltrd_anttd.vcf | \
                   grep 'FUNCOTATION' | \
                   sed 's/ //g;s/>//g;s/"//g') # remove spaces, double quotes, and right arrows from the string.
echo ${funcotation_header##*\Funcotationfieldsare:} | \
     sed 's/|/\t/g' > \ # replace pipes '|' with tabs '\t'.
		 ${annotation_dir}/SRR5439568_annotated_variants_snps.tsv # adding the columns to the head of an empty file.

# --------------- for the indels, the same script except for file names --------------- #

 # 2) Indels
  # 2.1) Extracting the column's names from the header section of the VCF into a new empty file. 
funcotation_header=$(grep '^##' ${annotation_dir}/SRR5439568_hptpcl_indels_fltrd_anttd.vcf | \
                   grep 'FUNCOTATION' | \
                   sed 's/ //g;s/>//g;s/"//g') # remove spaces, double quotes, and right arrows from the string.
echo ${funcotation_header##*\Funcotationfieldsare:} | \
     sed 's/|/\t/g' > \ # replace pipes '|' with tabs '\t'.
		 ${annotation_dir}/SRR5439568_annotated_variants_indels.tsv # adding the columns to the head of an empty file.
```

- **Adding annotations values to the variants annotations table file** - ****First, we use GATK’s ***VariantsToTable*** function to withdraw the FUNCOTATION values from the INFO column. We drop the values into a temporary file, that will soon be deleted.
    - *-V → (variants) → a VCF file, enriched with annotations.*
    - *-O → (output)→ a file to output the values to.*
    - *-F → (fields) → The name of a standard VCF field or an INFO field to include in the output table.*
    
    Later, we manipulate the extracted values and add them to the variants annotations table file right below the header line that was established previously. Eventually, we delete the temporary file.
    
    ```bash
      # 1.2) Extracting the INFO column from the VCF, with the FUNCOTATION data only, and dropping it to a temporary file.
    gatk VariantsToTable \
        -V ${annotation_dir}/SRR5439568_hptpcl_snps_fltrd_anttd.vcf \
        -O ${annotation_dir}/SRR5439568_annotated_variants_snps_tmp.tsv \
        -F FUNCOTATION 
      # 1.3) Uniting the column's name and the data itself to form a table (tsv file).
    grep -v 'FUNCOTATION' ${annotation_dir}/SRR5439568_annotated_variants_snps_tmp.tsv | \
         sed 's/[][]//g;s/|/\t/g' \ # removing square brackets from the string, and replace pipes '|' with tabs '\t'.
    		 >> ${annotation_dir}/SRR5439568_annotated_variants_snps.tsv
    rm ${annotation_dir}/SRR5439568_annotated_variants_snps_tmp.tsv
    
    # --------------- for the indels, the same script except for file names --------------- #
    
      # 2.2) Extracting the INFO column from the VCF, with the FUNCOTATION data only, and dropping it to a temporary file.
    gatk VariantsToTable \
        -V ${annotation_dir}/SRR5439568_hptpcl_indels_fltrd_anttd.vcf \
        -O ${annotation_dir}/SRR5439568_annotated_variants_indels_tmp.tsv \
        -F FUNCOTATION 
      # 2.3) Uniting the column's name and the data itself to form a table (tsv file).
    grep -v "FUNCOTATION" ${annotation_dir}/SRR5439568_annotated_variants_indels_tmp.tsv | \
         sed 's/[][]//g;s/|/\t/g' \ # removing square brackets from the string, and replace pipes '|' with tabs '\t'.
    		 >> ${annotation_dir}/SRR5439568_annotated_variants_indels.tsv
    rm ${annotation_dir}/SRR5439568_annotated_variants_indels_tmp.tsv
    ```
    

Finally, we can examine the final file that contains the variants we have found and filtered, together with their functional annotations.

```bash
-->> less ${annotation_dir}/SRR5439568_annotated_variants_indels_tmp.tsv
Gencode_27_hugoSymbol   Gencode_27_ncbiBuild    Gencode_27_chromosome   Gencode_27_start        Gencode_27_end  Gencode_27_variantClassific...
AGRN    hg38    chr1    1053552 1053552 INTRON          SNP     G       G       C       g.chr1:1053552G>C       ENST00000379370.6       +  ...
MXRA8   hg38    chr1    1361438 1361438 FIVE_PRIME_UTR          SNP     C       C       G       g.chr1:1361438C>G       ENST00000477278.3  ...
MTHFR   hg38    chr1    11791061        11791061        INTRON          SNP     A       A       G       g.chr1:11791061A>G      ENST0000037...
MTHFR   hg38    chr1    11794400        11794400        SILENT          SNP     G       G       A       g.chr1:11794400G>A      ENST0000037...
MTHFR   hg38    chr1    11794419        11794419        MISSENSE                SNP     T       T       G       g.chr1:11794419T>G      ENS...
MTHFR   hg38    chr1    11794698        11794698        INTRON          SNP     G       G       A       g.chr1:11794698G>A      ENST0000037...
```

### Diagnosis

The annotated variant table contains columns that described them. For instance, the clinical significance of the variants (benign, pathogenic, etc), variant classification (missense, splice site, etc), and more. We can use these columns to minimize the table, to only variants whose functional annotations suit the clinical disorder which our case deals with. Here, I included the following prerequisites:

- **Column’s names** - Include columns names will appear when we examine the results
- **Clinical significance:**
    - Include only records, which are considered **pathogenic**.
    - Exclude records that are annotated as **benign**.
- **Variant location** - Exclude variants in intronic realms.

```bash
-->> cat ${annotation_dir}/SRR5439568_annotated_variants_indels.tsv | \ # checking the annotated indels file
		 grep -E "Pathogenic|pathogenic|Gencode_27_hugoSymbol" | \
		 grep -v -E "Benign|benign|synonymous|INTRON"

Gencode_27_hugoSymbol   Gencode_27_ncbiBuild    Gencode_27_chromosome   Gencode_27_start        Gencode_27_end  Gencode_27_variantClassi

-->> cat ${annotation_dir}/SRR5439568_annotated_variants_snps.tsv | \ # checking the annotated SNPs file
		 grep -E "Pathogenic|pathogenic|Gencode_27_hugoSymbol" | \
		 grep -v -E "Benign|benign|synonymous|INTRON"

Gencode_27_hugoSymbol   Gencode_27_ncbiBuild    Gencode_27_chromosome   Gencode_27_start        Gencode_27_end  Gencode_27_variantClassification        Gen...
GALT      hg38    chr9    34647227        34647227        MISSENSE        SNP       T    T    C    g.chr9:34647227T>C      ENST00000556278.1       + ...
GALT      hg38    chr9    34648170        34648170        SPLICE_SITE     MISSENSE  SNP  A    A    G       g.chr9:34648170A>G      ENST00000556278.1 ...
```

**GALT** - The analysis reveals two different variants appearing at the same gene. It is a metabolic gene, that specifically holds a major role within glycolysis - **galactose-1-phosphate uridylyltransferase**. Its main job - converting galactose to glucose. The absence of this gene is causative for classic Galactosemia. In newborns who suffer from this disease, the consumption of lactose is life-threatening. **Galactosemia** is an autosomal recessive inherited rare disorder. In the IGV screenshot, both variants can be seen, both are heterozygous. For a disorder to take place, an individual should possess deleterious mutations in both chromosomes. Because the described symptoms of the newborn coincide with some of the classic galactosemia symptoms (vomiting, diarrhea, and jaundice), there is a great possibility that each of his chromosomes consists of a different mutation. That case is reasonable even with healthy parents. The other option, in it, the two variants are linked together on the same chromosome is also reasonable. In the latter option, the health situation of the infant is not related to the GALT variants.

![Untitled](Variant%20Discovery%20-%20GitHub%20repository%20-%202nd%20editio%20c1c43825599c4bcaa255c7c6c364817d/Untitled%204.png)