#!/bin/bash

# Variables
####################################################################################################

min_read_size=50
min_phred=30
strand="reverse"

# Directory stuff
####################################################################################################
cd ..
source rnaseq_preproc/bin/activate

data_dir="data/"
ref_dir="${data_dir}reference/"
ercc_dir="${data_dir}ERCC92/"
idx_dir="${ref_dir}human_idx/"
ercc_idx_dir="${ercc_dir}ercc_idx/"
results_dir="results/preprocessing/"

qc_pretrim="${results_dir}qc_pretrim/"
qc_postrim="${results_dir}qc_postrim/"
trim_dir="${results_dir}trimmed/"
align_dir="${results_dir}aligned/"
align_ref_dir="${align_dir}ref/"
align_ercc_dir="${align_dir}ercc/"
counts_dir="${results_dir}counts/"
counts_ref_dir="${counts_dir}ref/"
counts_ercc_dir="${counts_dir}ercc/"

mkdir -p $ref_dir
mkdir -p $ercc_dir
mkdir -p $results
mkdir -p $qc_pretrim
mkdir -p $qc_postrim
mkdir -p $trim_dir
mkdir -p $idx_dir
mkdir -p $ercc_idx_dir
mkdir -p $align_dir
mkdir -p $align_ercc_dir
mkdir -p $align_ref_dir
mkdir -p $counts_dir
mkdir -p $counts_ref_dir
mkdir -p $counts_ercc_dir

####################################################################################################
# Pipeline                                                                                         #
####################################################################################################

# Download references
####################################################################################################

cd $ref_dir
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz

gzip GCF_000001405.26_GRCh38_genomic.fna.gz
gzip GCF_000001405.26_GRCh38_genomic.gff.gz
cd ../..

cd $ercc_dir
wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip

unzip ERCC92.zip
cd ../..

# Download FASTQs
####################################################################################################
cd $data_dir

while read -r accession; do
  echo "Downloading accession: $accession"
  prefetch $accession
  fasterq-dump $accession | gzip
done < accessions.txt

cd ..

# Preprocessing
####################################################################################################

## Transform annotation gff to gtf (cufflinks)
gffread "${ref_dir}GCF_000001405.26_GRCh38_genomic.gff" -T \
        -o "${ref_dir}GCF_000001405.26_GRCh38_genomic.gtf"

# Create an index with the genome file
if [ ! "$(ls -A $idx_dir)" ]; then
  STAR --runThreadN 8 \
       --runMode genomeGenerate \
       --genomeDir $idx_dir \
       --genomeFastaFiles "${ref_dir}GCF_000001405.26_GRCh38_genomic.fna" \
       --sjdbGTFfile "${ref_dir}GCF_000001405.26_GRCh38_genomic.gtf" \
       --sjdbOverhang 72
fi

# Create an index for the ERCC file
if [ ! "$(ls -A $ercc_idx_dir)" ]; then
  STAR --runThreadN 8 \
       --runMode genomeGenerate \
       --genomeDir $ercc_idx_dir \
       --genomeFastaFiles "${ercc_dir}ERCC92.fa" \
       --sjdbGTFfile "${ercc_dir}ERCC92.gtf" \
       --sjdbOverhang 72 \
       --genomeSAindexNbases 7
fi

# Read accessions.txt and for each accession run fastqc, trim reads and align with STAR
while read -r acc; do
  # Run FASTQC
  echo "Processing: $acc"
  fastqc "${data_dir}${acc}_1.fastq" "${data_dir}${acc}_2.fastq"  -o $qc_pretrim
  # Trim
  trim_galore --paired --quality $min_phred --length $min_read_size --phred33 --fastqc \
              --fastqc_args "-o ${qc_postrim}" -o $trim_dir \
              "${data_dir}${acc}_1.fastq" "${data_dir}${acc}_2.fastq"
  # Align trimmed reads to the index
  echo "Aligning ${acc} to reference genome..."
  align_ref_samp_dir="${align_ref_dir}${acc}/"
  mkdir -p $align_ref_samp_dir
  STAR --genomeDir $idx_dir \
       --readFilesIn "${trim_dir}${acc}_1_val_1.fq" "${trim_dir}${acc}_2_val_2.fq"\
       --runThreadN 10 \
       --outFileNamePrefix "${align_ref_samp_dir}${acc}_" \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMunmapped Within KeepPairs
done < "${data_dir}accessions.txt"

# Obtain alignment stats with picard
while read -r acc; do
  echo "Obtaining alignment stats of ${acc} with Picard..."
  align_ref_samp_dir="${align_ref_dir}${acc}/"
  # Perform some QCs
  picard CollectAlignmentSummaryMetrics \
         R="${ref_dir}GCF_000001405.26_GRCh38_genomic.fna" \
         I="${align_ref_samp_dir}${acc}_Aligned.sortedByCoord.out.bam" \
         O="${align_ref_samp_dir}${acc}_picard_alignment_summary_metrics.txt"

  picard CollectInsertSizeMetrics \
         I="${align_ref_samp_dir}${acc}_Aligned.sortedByCoord.out.bam" \
         O="${align_ref_samp_dir}${acc}_picard_insert_size_metrics.txt" \
         H="${align_ref_samp_dir}${acc}_insert_size_histogram.pdf" \
         M=0.5
done < "${data_dir}accessions.txt"

# Create FASTQs of unmapped reads and align them to the ERCC index
while read -r acc; do
  echo "Aligning ${acc} to ERCC..."
  # Get unmapped reads into a fastqs (paired ends)
  align_ref_samp_dir="${align_ref_dir}${acc}/"
  samtools fastq -f 4 \
           -1 "${align_ref_samp_dir}unmapped_1.fq" \
           -2 "${align_ref_samp_dir}unmapped_2.fq" \
           "${align_ref_samp_dir}Aligned.sortedByCoord.out.bam"

  # Align unmapped reads to ERCC index
  align_ercc_samp_dir="${align_ercc_dir}${acc}/"
  mkdir -p $align_ercc_samp_dir
  STAR --genomeDir $ercc_idx_dir \
       --readFilesIn "${align_ref_samp_dir}unmapped_1.fq" "${align_ref_samp_dir}unmapped_2.fq"\
       --runThreadN 10 \
       --outFileNamePrefix "${align_ercc_samp_dir}${acc}_" \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMunmapped Within KeepPairs
done < "${data_dir}accessions.txt"


# Obtain counts from the reference alignments and parse them in a CSV matrix
Rscript scripts/get_counts.R $align_ref_dir \
        --annotFile "${ref_dir}GCF_000001405.26_GRCh38_genomic.gtf" \
        --att "gene_name" \
        --strand $strand \
        --outDir $counts_ref_dir 

# Obtain counts from the ERCC alignments and parse them in a CSV matrix
Rscript scripts/get_counts.R $align_ercc_dir \
        --annotFile "${ercc_dir}ERCC92.gtf" \
        --att "gene_id" \
        --strand $strand \
        --outDir $counts_ercc_dir 

# Obtain MultiQC reports
multiqc $(find $qc_pretrim -type f -name "*.zip" -not -name "*_val_*") \
        --outdir "$qc_pretrim" -f

multiqc $(find $qc_postrim -type f -name "*.zip" -and -name "*_val_*") \
        $(find $align_ref_dir -type f -name "*_Log.final.out") \
        $(find $align_ref_dir -type f -name "*picard*") \
        --outdir "$qc_postrim" -f

# Differential expression
####################################################################################################
DE_dir="results/DE/"
mkdir -p $DE_dir

Rscript scripts/DE.R --gene_counts "${counts_ref_dir}counts.csv" \
        --ercc_counts "${counts_ercc_dir}counts.csv" \
        --sample_info "${data_dir}sample_info.csv" \
        --useERCCs \
        --outDir "${DE_dir}w_ercc/"

Rscript scripts/DE.R --gene_counts "${counts_ref_dir}counts.csv" \
        --ercc_counts "${counts_ercc_dir}counts.csv" \
        --sample_info "${data_dir}sample_info.csv" \
        --outDir "${DE_dir}wo_ercc/"

# Functional enrichment
####################################################################################################
enrich_dir="results/enrichment/"
mkdir -p $enrich_dir

Rscript scripts/enrichment.R "${DE_dir}w_ercc/DESeqRes_LFCShrunk.csv" \
        --alpha 0.05 \
        --outDir "${enrich_dir}w_ercc/"

Rscript scripts/enrichment.R "${DE_dir}wo_ercc/DESeqRes_LFCShrunk.csv" \
        --alpha 0.05 \
        --outDir "${enrich_dir}wo_ercc/"

deactivate