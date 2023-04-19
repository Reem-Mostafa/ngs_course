#!/bin/bash
#This bash script is about creating a functioning pipeline to analyze raw sequencing NGS data.


#1) First we need to create a directory tree to organize data:
#cd ~/ngs_course
#mkdir dnaseq && cd dnaseq

#mkdir data meta results logs


#Then we create a subdirectory to store raw files.
#cd data

#mkdir trimmed_fastq untrimmed_fastq

#Then We need to download the raw sequence files in the right directory (untrimmed_fastq)
cd data/untrimmed_fastq

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz

#Download bed file in the data directory
cd ..
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

#Then we download the Reference Genome file in the reference directory
#mkdir reference && cd reference

#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#cd ..
#This has been previously done to save time



#2)FASTQC
cd untrimmed_fastq

#Run FASTQC on each file

#This command line uses fastqc tool (v0.12.1) that runs QC on fastq files. The -t specifies the number of threads or processors running the analysis. The wildcard *.fastq.qz refers to the file extension of the raw data files.
fastqc -t 4 *.fastq.qz

#The html output results are moved to fastqc_untrimmed_reads and can be viewed via FileZilla (showing the sequence quality per base. per sequence, GC content, overrepresented sequences, among others)
#mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads



#To unzip the 2 fastqc output files

for zip in *.zip; do unzip $zip; done

#To save, concatenate the fastqc summary.txt into one folder
cat */summary.txt > ~/ngs_course/dnaseq/logs/fastqc_summaries.txt



#3)TRIMMOMATIC:v0.39-1
#First we need to decompress the fastq.qz files
cd ../..
cd data/untrimmed_fastq
zcat NGS0001.R1.fastq.qz > NGS0001.R1.fastq
zcat NGS0001.R2.fastq.qz > NGS0001.R2.fastq

#Now we can run trimmomatic version 0.39
#Trimmomatic is a tool for trimming and filtering reads from paired-end high-throughput sequencing.This command shows that 4 threads are being used in this process,only reads above phred quality score 33 are included, -baseout specifies the base name for output files and ILLUMINACLIP specifies the adapter sequences used for trimming.TRAILING:25: specifies that any bases with a quality score below 25 should be trimmed from the end of the reads.
#MINLEN:50: This option specifies that any reads that are shorter than 50 bases after trimming should be discarded.

bin/trimmomatic PE -threads 4 -phred33 /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq -baseout /home/ubuntu/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R.fastq ILLUMINACLIP:share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50

#The trimmomatic output is directed to trimmed_fastq
#Now we run fastqc on the trimmed files and direct the output to a new directory called fastqc_trimmed_reads
#mkdir ~/ngs_course/dnaseq/results/fastqc_trimmed_reads
cd ..
cd trimmed_fastq
fastqc -t 4 *.fastq
mv *.zip ~/ngs_course/dnaseq/results/fastqc_trimmed_reads/
mv *.html ~/ngs_course/dnaseq/results/fastqc_trimmed_reads/


#The fastqc output can be viewed using FileZilla like the previous one



#4)ALIGNMENT:
#Here we are going to index the Reference genome file to make alignment easier
#bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz (done previously)

#First we need to remove redundant and unneeded files to free up space:
cd ..
cd untrimmed_fastq
rm -r *.fastq.qz
rm -r *.fastq

#Align trimmed_fastq files to the indexed reference genome
cd ..
cd trimmed_fastq

bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1:111:D1375ACXX:1:1101\tSM:1101\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:11V6WR1' -I 250,50  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P.fastq ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P.fastq > ~/ngs_course/dnaseq/data/aligned_data/NGS0001.sam

#The output is a .sam file that has to be converted to .bam first since BAM files are compressed and binary-encoded, allowing for faster and more efficient storage and processing of large sequence alignment files.
cd ..
cd aligned_data
samtools view -b -h NGS0001.sam > NGS0001.bam

rm -r NGS0001.sam
#bam file is then sorted by genomic coordinates
samtools sort NGS0001.bam > NGS0001_sorted.bam

#Generate an index using samtools to allow for fast retrieval of data for specific regions of interest, such as a particular chromosome.
samtools index NGS0001_sorted.bam

#Mark duplicates. Duplicate reads arising from PCR amplification artefacts or sequencing biases can result in misleading analysis results. Picard MarkDuplicates identifies and marks these duplicate reads based on their alignment coordinates and other information such as the mapping quality and orientation.
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

#Index the sorted_marked_bam
samtools index NGS0001_sorted_marked.bam

#Filter BAM based on mapping quality and bitwise flags using samtools. "-F 1796": This option specifies the bitwise flag filter, indicating that reads with any of the specified flags indicating that reads that are either not primary alignments, fail quality checks, or are PCR or optical duplicates.
#"-q 20": This option specifies the mapping quality filter, indicating that only reads with a mapping quality score greater than or equal to 20 will be included in the output. Mapping quality scores range from 0 to 60, with higher scores indicating more confident mappings.
#"-o:" refers to the output file name, "NGS0001_sorted_filtered.bam", where the filtered BAM file will be written to.
#"NGS0001_sorted_marked.bam": This is the input BAM file that will be filtered based on the specified options.
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam

#Generate standard alignment statistics
#This command counts the number of reads in a BAM file that match different bitwise flag criteria and provides a summary report of the counts.
samtools flagstat NGS0001_sorted_filtered.bam

#This command is a utility in the samtools package that generates summary statistics for an indexed BAM file. It reports the number of reads that align to each reference sequence in the BAM file.The output contains four columns of information: the reference sequence name, the total length of the reference sequence, the number of reads in the BAM file that align to the reference sequence, and the number of reads in the BAM file that align to the reverse complement of the reference sequence.
samtools idxstats NGS0001_sorted_filtered.bam

#This command calculates the distribution of insert sizes for each library, where the insert size is defined as the distance between the 5' ends of the paired-end reads. This can be used to assess the quality of the sequencing library preparation and to identify any biases or anomalies in the insert size distribution that may affect downstream analysis.
picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5

#The output of this command (genomecov) is a text file that contains coverage information for each position in the reference genome. The file has columns:Chromosome, Position, Depth, Number of Bases, Percent Coverage (The percentage of bases at the position that are covered by at least one read),
#Depth (Forward Strand) and Depth (Reverse Strand).
bedtools genomecov -ibam NGS0001_sorted_filtered.bam > coverage.txt

#Remove some unwanted files:
rm -r NGS0001.bam
rm -r NGS001_sorted.bam
rm -r NGS0001_sorted_marked.bam

#5)Variant calling with FreeBayes

cd ..
cd reference
#First we have to uncompress the reference file hg19
zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa

#Then create an index for the Reference Genome in Fasta using Samtools (.fai)
samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa

#Call variants with Freebayes while restricting analysis to the target regions in the bed file using the --targets command. PS:It did not work on my version of freebayes (v0.9.21) so I updated it with 'conda install -c "bioconda/label/cf201901" freebayes'
freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --targets ~/ngs_course/dnaseq/data/annotation.bed --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/results/NGS0001.vcf

#Compress the resulting vcf and index it with tabix.
bgzip ~/ngs_course/dnaseq/results/NGS0001.vcf

tabix -f vcf ~/ngs_course/dnaseq/results/NGS0001.vcf.gz

#The VCF output is then filtered with 'vcffilter' to remove low quality calls, low-depth reads, alleles seen on one strand and unbalanced alleles.The -f option is used to specify the filtering criteria.
#QUAL > 1: This specifies that the variant must have a quality score greater than 1.
#QUAL / AO > 10: This specifies that the ratio of the quality score to the allele depth must be greater than 10.
#SAF > 0 & SAR > 0: This specifies that both the number of alternate alleles on the forward strand (SAF) and the number of alternate alleles on the reverse strand (SAR) must be greater than 0.
#RPR > 1 & RPL > 1: This specifies that both the number of reference reads on the forward strand (RPR) and the number of reference reads on the reverse strand (RPL) must be greater than 1.
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/ngs_course/dnaseq/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf


#PS:I had to install vcflib first using 'conda install vcflib'

#Compress and index the filtered vcf
bgzip ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf

tabix -p vcf NGS0001_filtered.vcf

#6)ANNOTATION WITH ANNOVAR
#Downloading annovar and uploading it via filezilla followed by uncompressing the tar.gz file using this command was done first: 'tar -zxvf annovar.latest.tar.gz'
#Annovar databases were then downloaded
#The vcf output from freebayes is converted into an annovar compatible input using 'convert2annovar' .
cd
cd annovar
./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.avinput

# Run Annovar table function to produce a csv file 'NGS0001_annovar'.
./table_annovar.pl ~/ngs_course/dnaseq/results/NGS0001_filtered.avinput humandb/ -buildver hg19 -out ~/ngs_course/dnaseq/results/NGS0001_annovar -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

#Download the resulting csv file via Filezilla and view it in excel.



#7)ANNOTATION WITH snpEFF:
#First download SnpEff using 'wget http://sourceforge.net/projects/snpeff/files/snpEff_v4_3_core.zip'and the hg19 database using 'java -jar snpEff.jar download -v hg19'
#Then run snpEff using the following command. - Xmx8g: sets the maximum heap size for the Java virtual machine to 8 gigabytes. This means that the Java program will be allowed to use up to 8 GB of RAM.:
cd
cd snpEff
java -Xmx8g -jar snpEff.jar hg19 ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf > snpeff_output.vcf

#Summary stats can be viewed via FileZilla.


