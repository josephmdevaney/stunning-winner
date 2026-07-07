#!/usr/bin/env bash
# Fetch real HG002 PacBio HiFi reads at the 80 STRchive loci from NIST/GIAB's public
# S3 bucket (remote-indexed, no full download), then genotype them with the real
# trgt v5.1.0 binary. Produces fixtures/hg002_real/HG002.strchive_loci.trgt.vcf.gz.
#
# Requires: samtools, curl, a trgt v5.1.0 binary on PATH (or set TRGT_BIN), and
# strchive_trgt.bed (see build_strchive_trgt_catalog.py) in the working directory.
set -euo pipefail

BAM_URL="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002_GRCh38/HG002_GRCh38.haplotag.10x.bam"
TRGT_BIN="${TRGT_BIN:-trgt}"
WORKDIR="${1:-./hg002_real_workdir}"
BED="strchive_trgt.bed"

mkdir -p "$WORKDIR/region_bams"
cd "$WORKDIR"

# 1. Padded per-locus regions, then slice each from the remote BAM (indexed range
#    requests -- never downloads the full 69GB file).
awk -F'\t' '{print $1":"($2-2000<0?0:$2-2000)"-"($3+2000)}' "../$BED" > regions.txt
i=0
while read -r region; do
  i=$((i+1))
  safe=$(echo "$region" | tr ':-' '__')
  samtools view -b "$BAM_URL" "$region" > "region_bams/${i}_${safe}.bam"
done < regions.txt

samtools merge -f merged_unsorted.bam region_bams/*.bam
samtools sort -o hg002_strchive_loci.bam merged_unsorted.bam
samtools index hg002_strchive_loci.bam

# 2. N-padded reference FASTA with real GRCh38 sequence patched in only at the loci
#    windows (see build_padded_reference.py for why full chromosome lengths matter).
curl -s "https://giab.s3.amazonaws.com/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz.fai" -o ref.fasta.gz.fai
cp "../$BED" strchive_trgt.bed
python3 "$(dirname "$0")/build_padded_reference.py"
samtools faidx padded_reference.fasta

# 3. Real trgt v5.1.0 genotyping.
"$TRGT_BIN" genotype \
  --genome padded_reference.fasta \
  --reads hg002_strchive_loci.bam \
  --repeats strchive_trgt.bed \
  --output-prefix hg002_real \
  --sample-name HG002 \
  --karyotype XX \
  --threads 4

cp hg002_real.vcf.gz ../fixtures/hg002_real/HG002.strchive_loci.trgt.vcf.gz
echo "done: fixtures/hg002_real/HG002.strchive_loci.trgt.vcf.gz"
