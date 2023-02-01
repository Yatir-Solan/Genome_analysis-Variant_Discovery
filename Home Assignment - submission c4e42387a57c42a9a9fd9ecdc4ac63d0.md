# Home Assignment - submission

# Home Assignment

The assignment was implemented as a bash script that was run on a remote Linux machine.

The lines given in the code patches are brought in to represent the pipeline. 

### Disclaimers

- Even though I had two courses in that area of research, I have never done this kind of research hands-on. Therefore, after my initial literature survey and reading I decided to use common tools that I was briefly introduced to before, (**bcftools** for example). Mainly, because it was more convenient for me. After observing much of the literature… if I had started all over again, I would probably consider working with other tools. For instance, more modern tools, which also provide more of a full pipeline (ready out of the box). e.g., **nf-core/sarek** or **GATK**. Or at least, consider utilizing parts in the pipeline with more modern tools such as **DeepVariant**, or **FreeBayes**.
- Some of the bash scriptings in the code lines attached, can be more efficient both in code lines (by Linux piping), and in memory resources (for instance, I could have maintained most of the commands with files still zipped). I did not follow the efficient approach - first because the files are relatively small, and second because I had a need to look at the files along the pipeline ‘with my own eyes’ to gain confidence in what I’m doing.
- Time limitation - As I’m describing in the following lines. Some important stages were not done, due to time limitations. Along this paper, whenever a case like this appears, I am noting it.
- Scientific writing - I’m aware of the non-scientific writing I have adopted here. If needed, sources could be delivered.

### Loading the relevant programs

```bash
module load AWScli/aws-cli-2.1.0
module load sratoolkit.3.0.1
module load bwa/bwa-0.7.17
module load samtools/1.9
module load bcftools/bcftools-1.6
```

### Reference Genome

- **Reference genome** **download** - Because the assignment is a human genome analysis a human reference genome was needed. I have found that the latest reference genome is the GRCh38. After trying several versions of the GRCh38, I decided eventually to use the GATK version of it. The genome was downloaded from [here](https://ewels.github.io/AWS-iGenomes/) by the following bash lines.

```bash
# Downloading GRCh38 reference genome
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ ./references/Homo_sapiens/GATK/GRCh38/
```

- **Reference genome indexing** - Since I have used a GATK ‘ready to use’ product, the reference genome was already indexed. Nevertheless, in order for the pipeline to work, index files needed to be relocated to be with the fasta files.

```bash
# copying the reference genome index files, so they will be in the same directory as the genome itself.
cd /references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta
cp ./../BWAIndex/* ./ 
```

- **Hypothetical indexing** - In the case of using an un-indexed reference genome, we can use the ‘bwa index’ command in order to index it. That can be done by the following line.

```bash
# Hypthetical - was not run.
bwa index /references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
```

### Raw reads

- **Reads download** - The raw reads were obtained using ‘sratoolkit’, from the given [link](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5439568&display=metadata), using the accession SRR5439568.

```bash
mkdir reads
cd reads
/powerapps/share/centos7/miniconda/miniconda3-4.7.12-environmentally/envs/sratoolkit.3.0.1/bin/fasterq-dump --split-files SRR5439568
```

- **Reads inspection**
    - **Reads details**
        - The reads are the result of Illumina 1.9 sequencing.
        - Paired-End, 4,869,926 reads each.
        - A rather small file (1.6 GB each) - which is consistent with an individual sample.
        - Reads length - 101.
    - **Quality Control** - To assess the quality of the reads I have used the ‘FastQC’ program in the command line.
    
    ```bash
    fastqc reads/SRR5439568_1.fastq reads/SRR5439568_2.fastq
    ```
    
    - **FastQC results** (Print screens of one of the files are given when needed. Through all the measurements - the QC results of both files look generally the same.).
        - **Base quality scores** - above 30 (green zone) for all reads (both files), and there is no significant quality decay toward the margins of the reads.
        
        ![Untitled](Home%20Assignment%20-%20submission%20c4e42387a57c42a9a9fd9ecdc4ac63d0/Untitled.png)
        
        - **Mean sequence quality scores** - The vast majority of the reads maintain mean quality scores that are above 30.
        - **Per base sequence content** - Both files show a rather unbiased distribution. (to my knowledge the slight fluctuation within the beginning is fairly normal).
        
        ![Untitled](Home%20Assignment%20-%20submission%20c4e42387a57c42a9a9fd9ecdc4ac63d0/Untitled%201.png)
        
        - **Per sequence GC content** - **Warning**. A warning is raised if the sum of the deviations from the normal distribution represents more than 15% of the reads. **Common reasons for warnings -** Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example), which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species. (Taken from [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html)).
        
        ![Untitled](Home%20Assignment%20-%20submission%20c4e42387a57c42a9a9fd9ecdc4ac63d0/Untitled%202.png)
        
        - **Per base N content** - None for both.
        - **Sequence length distribution** - 101 for all, in both.
        - **Sequence duplication levels - Failure.** This module will issue a error if non-unique sequences make up more than 50% of the total. (Taken from [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)).
        
        ![Untitled](Home%20Assignment%20-%20submission%20c4e42387a57c42a9a9fd9ecdc4ac63d0/Untitled%203.png)
        
        - **Overrepresented sequences - Warning.** This module will issue a warning if any sequence is found to represent more than 0.1% of the total. **Common reasons for warnings -** This module will often be triggered when used to analyze small RNA libraries where sequences are not subjected to random fragmentation, and the same sequence may naturally be present in a significant proportion of the library. (Taken from [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html)).
        - **Adaptor content** - None for both.
    - **Quality Control conclusion** - The sequencing is in a good shape, but it needs to undergo duplicate reads removal.
    
    ### **Duplicated removal using SAMBAMBA**
    
    Due to time constraints, I have not managed to learn this program and utilize it. I believe if I had the time to conduct a filter upon duplicate reads, the next FastQC report regarding it, would have been much better concerning both reads in over-representation, duplication levels, and Per sequence GC content).
    
    ### Mapping the reads and preparation for variants calling
    
    - **SM value** - Even though it is considered a good practice, due to having individual sequencing, I have not added an SM value.
    - **BWA mem** - Both BWA and STAR were considered. Eventually, BWA was chosen because it is far faster than STAR. Efficiency was also considered when choosing the aligner algorithm within BWA. Therefore, ‘mem’ algorithm was chosen.
        - -t 8 → 8 threads.
    
    ```bash
    bwa mem -t 8 references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta reads/SRR5439568_1.fastq reads/SRR5439568_2.fastq -o SRR5439568.sam
    ```
    
    - **SAM → BAM** - Converting the SAM file into a BAM file.
        - -b → BAM output.
        - -h → include header in the output.
    
    ```bash
    samtools view -b -h SRR5439568.sam > SRR5439568.bam
    ```
    
    - **Sorting the BAM file**
    
    ```bash
    samtools sort SRR5439568.bam > SRR5439568_sorted.bam
    ```
    
    - **Indexing the BAM file** - This is mainly crucial for visualization of the variants later on within IGV.
    
    ```bash
    samtools index SRR5439568_sorted.bam
    ```
    
    - **One-liner** - In order to save memory, it is also possible to combine all the commands into a one-liner.
        - -t 8 → 8 threads.
    
    ```bash
    # one liner - was not in use.
    bwa mem -t 8 references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta reads/SRR5439568_1.fastq reads/SRR5439568_2.fastq -o SRR5439568.sam | samtools sort -o SRR5439568_sorted.bam -
    ```
    
    - **Alignment quality control** - Using samtools
    
    ```bash
    samtools flagstat SRR5439568_sorted.bam
    # output ->
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
    
    - **Filtering reads** - By using ‘samtools view’ command I have filtered reads that were unmapped and mapped not-properly. The differences can be seen in the new statistics output.
        - -b → BAM output.
        - -F 0x04 → Filters all un-mapped reads.
        - -f 0x2 → Includes properly paired only mapped reads.
    
    ```bash
    samtools view -b -F 0x04 -f 0x2 -o SRR5439568_sorted_filtered.bam SRR5439568_sorted.bam
    
    samtools flagstat SRR5439568_sorted_filtered.bam
    # output ->
    9378607 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    110203 + 0 supplementary
    0 + 0 duplicates
    9378607 + 0 mapped (100.00% : N/A)
    9268404 + 0 paired in sequencing
    4634202 + 0 read1
    4634202 + 0 read2
    9268404 + 0 properly paired (100.00% : N/A)
    9268404 + 0 with itself and mate mapped
    0 + 0 singletons (0.00% : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)
    ```
    
    ### Variant Calling
    
    - **bcftools mpileup** - provides a summary of the coverage of the mapped reads, to the reference genome.
        - -O b → BAM file output.
        - -f → reference genome (indexed).
    - **bcftools call**
        - -m → multiallelic caller. (the user should choose between multiallelic caller, and consensus caller). Unfortunately, there is no comprehensive data regarding the options, and no pros and cons details related to each. Nevertheless, a multiallelic caller is described as the one more suitable for rare-variant calling, and it is more modern and more common.
        - -v → output only variants.
    
    ```bash
    bcftools mpileup --threads 8 -f references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta -O b -o SRR5439568_raw.bcf SRR5439568_sorted.bam
    bcftools call -m -v SRR5439568_raw.bcf -o SRR5439568_raw_variants.vcf
    ```
    
    ### Variant filtering
    
    I decided to conduct hard filtering, rather than more sophisticated technics, such as VQSR.
    
    - -i → include (include when the expression is true).
    - QUAL → read quality. GATK recommends filtering every variant with QUAL less than 30. I followed this recommendation.
    - DP → depth. I decided to omit variants with depths lower than 10.
        - Later, I noticed this (the added quote from a relatively old GATK [documentation](https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/tutorials/2806-how-to-apply-hard-filters-to-a-call-set)), and wonder if I should remove the DP filtration. Because of the time limitation, I did not, but this issue should be reconsidered.
        
        > “The maximum DP (depth) filter only applies to whole genome data, where the probability of a site having exactly N reads given an average coverage of M is a well-behaved function. First principles suggest this should be a binomial sampling but in practice it is more a Gaussian distribution. Regardless, the DP threshold should be set a 5 or 6 sigma from the mean coverage across all samples, so that the DP > X threshold eliminates sites with excessive coverage caused by alignment artifacts. Note that for exomes, a straight DP filter shouldn’t be used because the relationship between misalignments and depth isn’t clear for capture data.”
        > 
    - MQ → the root mean square over all the reads within this sight. GATK considered variants above 40 MQ to be fine. But, surveying the web informed me that above 60-filtration is also very common, and scrutinizing my data led me to understand that the vast majority is above 60 anyway.
    
    ```bash
    bcftools filter -i 'QUAL>=30 & DP>=10 & MQ>=60' SRR5439568_raw_variants.vcf > SRR5439568_raw_variants_filtered.vcf
    ```
    
    ### Variant annotation
    
    To conduct a variant annotation I have used [wANNOAVR](https://wannovar.wglab.org/) to annotate the variants present in the filtered VCF file. I have submitted the VCF with the appropriate parameters and received as an input a CSV file, which I analyzed locally on my personal machine. 
    
    **I have filtered** **using Python**: 
    
    - Synonymous variants.
    - Benign annotations.
    - Uncertain or unknown annotations.
    
    ```python
    import pandas as pd 
    wnvr_df = pd.read_table(r'C:\Users\Yatir\OneDrive\code\Danny_Zeevi\wanvr_query_output_exome_summary.csv', sep=',')
    wnvr_df = wnvr_df[~(wnvr_df.loc[:,'ExonicFunc.refGene']=='synonymous SNV')]
    wnvr_df = wnvr_df[~(wnvr_df.loc[:,'ExonicFunc.refGene']=='unknown')]
    wnvr_df = wnvr_df[~(wnvr_df.loc[:,'ClinVar_SIG'].isnull())]
    wnvr_df = wnvr_df[~(wnvr_df.loc[:,'ClinVar_SIG'].str.contains('enign'))]
    wnvr_df = wnvr_df[~(wnvr_df.loc[:,'ClinVar_SIG'].isin(['Benign', 'benign', '.', 'not provided', 'Uncertain significance']))]
    wnvr_df.sort_values(by='Gene.refGene', inplace=True)
    wnvr_df.to_csv(r'C:\Users\Yatir\OneDrive\code\Danny_Zeevi\wanvr_query_output_exome_summary_filtered.csv', index=False, sep=',')
    
    # analyzing the data frame into informative (keeping the most interesting columns) :
    desired_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'dbSNP', 'ClinVar_SIG', 'ClinVar_DIS', 'Otherinfo']
    wnvr_df[desired_cols]
    ```
    
    | Chr | Start | End | Ref | Alt | Func.refGene | Gene.refGene | ExonicFunc.refGene | dbSNP | ClinVar_SIG | ClinVar_DIS | Otherinfo |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | chr1 | 100206504 | 100206504 | T | C | exonic | DBT | nonsynonymous SNV | rs12021720 | Pathogenic\x2cother | Intermediate_maple_syrup_urine_disease_type_2\... | hom |
    | chr9 | 34647227 | 34647227 | T | C | exonic | GALT | nonsynonymous SNV | rs111033663 | Pathogenic | Deficiency_of_UDPglucose-hexose-1-phosphate_ur... | het |
    | chr9 | 34648170 | 34648170 | A | G | exonic | GALT | nonsynonymous SNV | rs75391579 | Pathogenic|Pathogenic | Deficiency_of_UDPglucose-hexose-1-phosphate_ur... | het |
    | chr11 | 36574050 | 36574050 | A | G | exonic | RAG1 | nonsynonymous SNV | rs3740955 | Pathogenic | not_provided | het |
    | chr19 | 15879621 | 15879621 | C | T | exonic | CYP4F2 | nonsynonymous SNV | rs2108622 | drug response|drug response|drug response | acenocoumarol_response_-_Dosage|warfarin_respo... | het |
    
    ### Diagnosis
    
    - **RAG1** - is an heterozygous variant. The ClinVar accession RCV000536425.8 describes it as an autosomal recessive variant. Checking it on IGV affirms its heterozygosity.
    
    ![Untitled](Home%20Assignment%20-%20submission%20c4e42387a57c42a9a9fd9ecdc4ac63d0/Untitled%204.png)
    
    - **CYP4F2** - The data is not very clear. From the ClinVar accession RCV000211318.4, It seems that the gene is a Pharmacogenes. Specifically, having relation to drug response regarding Warfarin. As can be seen in IGV, its sequencing depth is not pleasing, but as far as can be seen from it - it is heterozygous.
    
    ![Untitled](Home%20Assignment%20-%20submission%20c4e42387a57c42a9a9fd9ecdc4ac63d0/Untitled%205.png)
    
    - **GALT** - My analysis reveals two different variants appearing at the same gene. It is a metabolic gene, that specifically holds a major role within glycolysis - **galactose-1-phosphate uridylyltransferase**. Its main job - converting galactose to glucose. The absence of this gene is causative for classic Galactosemia. In newborns who suffer from this disease, the consumption of lactose is life-threatening. **Galactosemia** is an autosomal recessive inherited rare disorder. In the IGV screenshot, both variants can be seen, both are heterozygous. For a disorder to take place, an individual should possess deleterious mutations in both chromosomes. Because the described symptoms of the newborn coincide with some of the classic galactosemia symptoms (vomiting, diarrhea, and jaundice), there is a great possibility that each of his chromosomes consists of a different mutation. That case is reasonable even with healthy parents. The other option, in it, the two variants are linked together on the same chromosome is also reasonable. In the latter option, the health situation of the infant is not related to the GALT variants.
    
    ![Untitled](Home%20Assignment%20-%20submission%20c4e42387a57c42a9a9fd9ecdc4ac63d0/Untitled%206.png)
    
    - **DBT** - That is an homozygous variant, which can be seen clearly in the IGV screenshot. The literature informs us, that the unfunctional DBT gene, is one of the causes of **Maple syrup urine disease (MSUD).** MSUD is an autosomal recessive, metabolic disorder, that is affecting branched-chain amino acids. The disease gets its name due to the sweet smell of urine and earwax of afflicted infants. The ClinVar accession ([RCV000012727.23](https://www.ncbi.nlm.nih.gov/clinvar/RCV000012727.23/)) related to the variant relates it to **Intermediate MSUD** called type 2. This is a disease variant that is less severe than the classic MSUD. When comparing the symptoms of the examined newborn, with common symptoms of those afflicted with Intermediate MSUD, it should be noted that the latter Intermediate MSUD is diagnosed later in life. In addition, the common symptoms of classic MSUD are rather different than the symptoms described in this case.
    
    ![Untitled](Home%20Assignment%20-%20submission%20c4e42387a57c42a9a9fd9ecdc4ac63d0/Untitled%207.png)
    
    ### Conclusions
    
    - Further examination of the CYP4F2 variant. First, to get a better sequencing for being way more sure regarding its actual presence. Later understand its nature - whether it is autosomal dominant or not. More, I would have checked if Warfarin was part of the newborn healthcare or its mother during pregnancy.
    - **Galactosemia** - Concerning the symptoms and the aggregated data, there is a great possibility that the child suffers from Galactosemia. First, urgently alter the child feeding to a lactose-free diet. Final acknowledgment of the disease should be obtained using NBS screening.
    - **Maple syrup urine disease** - There is a great possibility that the child suffers from it. Nevertheless, it is not quite sure if the disease is present now or not. That as well should be obtained precisely by screening. As an initial step, the child should move to a low-protein diet. This diet should be aimed to reduce amino-acids consumption, specifically leucine, valine, and isoleucinestart.