#!/bin/sh


############################################################
# download S. cerevisiae genome sequence and annotation
############################################################
# Ensembl R64-1-1
# ftp://ftp.ensembl.org/pub/release-84/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
# ftp://ftp.ensembl.org/pub/release-84/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.84.gtf.gz
# files were uncompressed manually into "./genome" folder



############################################################
# check quality of sequencing reads
############################################################
# FastQC version 0.11.5 (www.bioinformatics.babraham.ac.uk/projects/fastqc/)
# FASTQ files were downloaded manually from ArrayExpress/ENA to "./FASTQ" folder

FastQC_outdir="./FastQC_output/"
mkdir "${FastQC_outdir}"
FASTQ_dir="./FASTQ/"
ls -1 "${FASTQ_dir}" > FASTQ_files
FASTQ_files=`cat FASTQ_files`
for i in ${FASTQ_files};
do
	echo "############################################"
	date
	echo "FastQC processing file: ${FASTQ_dir}${i}"
	echo "############################################"
	fastqc --outdir "${FastQC_outdir}" --threads 4 "${FASTQ_dir}${i}"
	echo
done

# Inspect FastQC output & if needed, use tools such as Trimmomatic to remove adapters and poor-quality reads.



############################################################
# map reads to the genome
############################################################
# HISAT2 version 2.0.4 (https://ccb.jhu.edu/software/hisat2/index.shtml)
# samtools version 1.3.1 (http://www.htslib.org/)
#
# HISAT2 genome index dowloaded from ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/r64_tran.tar.gz
# Index archive was uncompressed manually to "./genome" folder
# Index contains genome sequence + transcript structures (S. cerevisiae, Ensembl R64-1-1).
#
# Single-end data are needed for downstream analyses, so only read_1 FASTQ files were used.
# Only keeping reads with mapping quality >= 10 (i.e., reads mapping well to a single location).
# Adjust the HISAT2 '--rna-strandness' parameter according to your sequencing library preparation protocol if needed.

min_intron=20
max_intron=10000
HISAT2_index="./genome/genome_tran"
BAM_dir="./BAM/"
mkdir "${BAM_dir}"
FASTQ_dir="./FASTQ/"
FASTQ_files=`cat FASTQ_files`
for i in ${FASTQ_files};
do
	infile="${FASTQ_dir}${i}"
	outfile="${BAM_dir}${i}.bam"
	echo "############################################"
	date
	echo "HISAT2 processing file: ${infile}"
	echo "############################################"
	hisat2 -x "${HISAT2_index}" -U "${infile}" --min-intronlen ${min_intron} --max-intronlen ${max_intron} --threads 4 --rna-strandness R | samtools view -b -q 10 --threads 4 - | samtools sort -o "${outfile}" - 
	samtools index "${outfile}"
	samtools flagstat "${outfile}"
	echo
done



############################################################
# extract transreads (spliced junction reads)
############################################################
# regtools 0.2.0 (https://regtools.readthedocs.io/en/latest/)

genome_fasta="./genome/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa"
genome_gtf="./genome/Saccharomyces_cerevisiae.R64-1-1.84.gtf"
regtools_outdir="./transreads/"
mkdir "${regtools_outdir}"
BAM_dir="./BAM/"
BAM_files=`ls "${BAM_dir}" | grep .bam$`
annotated_suffix=".annotated"
for i in $BAM_files;
do
	infile="${BAM_dir}${i}"
	outfile="${regtools_outdir}${i}.trans"
	echo "############################################"
	date
	echo "regtools processing file: ${infile}"
	echo "############################################"
	regtools junctions extract -i ${min_intron} -I ${max_intron} -o "${outfile}" "${infile}"
	regtools junctions annotate -o "${outfile}${annotated_suffix}" "${outfile}" "${genome_fasta}" "${genome_gtf}"
	echo
done



############################################################
# compile transread counts and extract splice site coordinates
############################################################
# R version 3.3.1 (https://www.r-project.org/)

Rscript junctions.R



############################################################
# extract splice site (intron) coverage
############################################################
# bedtools version 2.25.0 (http://bedtools.readthedocs.io/en/latest/)
# Adjust the -S/-s "strandness" parameter according to your sequencing library preparation protocol.

ss_dir="./genome/"
ss5_file="introns_known_5ss.bed"
ss3_file="introns_known_3ss.bed"
counts_suffix=".counts"
bedtools_outdir="./introns/"
mkdir "${bedtools_outdir}"
BAM_dir="./BAM/"
BAM_files=`find ./ | grep ${BAM_dir} | grep .bam$`
bedtools multicov -S -split -bed "${ss_dir}${ss5_file}" -bams ${BAM_files} > "${bedtools_outdir}${ss5_file}${counts_suffix}"
bedtools multicov -S -split -bed "${ss_dir}${ss3_file}" -bams ${BAM_files} > "${bedtools_outdir}${ss3_file}${counts_suffix}"
ls "${BAM_dir}" | grep .bam$ > "BAM_files"



############################################################
# calculate splicing efficiency and plot data
############################################################
# R version 3.3.1 (https://www.r-project.org/)

mkdir "./images/"
mkdir "./efficiency/"
Rscript efficiency.R
