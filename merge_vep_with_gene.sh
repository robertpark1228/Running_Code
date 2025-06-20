#!/bin/bash

# Usage:
# ./merge_vep_with_genelevel.sh vep_output.tsv dbNSFP5.1_gene.gz merged_output.tsv

# Inputs
VEP_FILE="$1"
GENE_FILE="$2"
OUTPUT_FILE="$3"

# Check input
if [[ ! -f "$VEP_FILE" ]] || [[ ! -f "$GENE_FILE" ]]; then
  echo "❌ Missing input file. Usage: ./merge_vep_with_genelevel.sh <vep_output.tsv> <dbNSFP5.1_gene.gz> <output.tsv>"
  exit 1
fi

# Extract gene column name from VEP (Gene or Gene_Name)
GENE_COL=$(head -n1 "$VEP_FILE" | tr '\t' '\n' | grep -E '^Gene(_Name)?$')
if [[ -z "$GENE_COL" ]]; then
  echo "❌ Could not find a Gene or Gene_Name column in $VEP_FILE"
  exit 1
fi

# Decompress gene-level data
zcat "$GENE_FILE" > dbNSFP5.1_gene.tsv

# Merge using csvtk (recommended), or fallback to join if installed
if command -v csvtk &> /dev/null; then
  echo "✅ Using csvtk to merge on column $GENE_COL"
  csvtk join -t -f "${GENE_COL},Gene_name" "$VEP_FILE" dbNSFP5.1_gene.tsv > "$OUTPUT_FILE"
else
  echo "❌ csvtk not installed. Please install it with: conda install -c bioconda csvtk"
  exit 1
fi

echo "✅ Merged file created: $OUTPUT_FILE"
