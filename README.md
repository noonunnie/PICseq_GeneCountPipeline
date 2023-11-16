# PICseq dry lab pipeline
Raw fastq to cnv pipeline for PICseq data.

R1 holds the spatial barcode.
R2 holds the cDNA data.

# Required packages
fastqc - https://github.com/s-andrews/FastQC
  cutadapt

trimgalore - https://github.com/FelixKrueger/TrimGalore)https://github.com/FelixKrueger/TrimGalore

hisat2 - http://daehwankimlab.github.io/hisat2/download/

samtools - https://github.com/samtools/samtools

sed

umi-tools - https://github.com/CGATOxford/UMI-tools


# Must have downloaded
.sra sample file(s)

whitelist file of desired barcodes (one row per barcode)

gencode.v34
