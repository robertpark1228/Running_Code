#!/bin/bash

# Input VCF
VCF_IN=${1}
CSV_OUT=${VCF_IN}_annotate.csv

# Extract the CSQ field format from the VCF header
CSQ_FORMAT=$(bcftools view -h "$VCF_IN" | grep "ID=CSQ" | sed -E 's/.*Format: //; s/">//')
IFS='|' read -r -a CSQ_FIELDS <<< "$CSQ_FORMAT"

# Print CSV header
echo -e "CHROM,POS,REF,ALT,DP,${CSQ_FORMAT//|/,}" > "$CSV_OUT"

# Extract fields with bcftools and process with awk
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/CSQ\n' "$VCF_IN" | \
awk -F'\t' -v OFS=',' -v csq_format="$CSQ_FORMAT" '
BEGIN {
    n = split(csq_format, fields, "|");
}
{
    chrom = $1; pos = $2; ref = $3; alt = $4; dp = $5; csq = $6;
    split(csq, csqs, ",");
    for (i in csqs) {
        split(csqs[i], values, "|");
        printf "%s,%s,%s,%s,%s", chrom, pos, ref, alt, dp;
        for (j = 1; j <= n; j++) {
            printf ",%s", (j in values ? values[j] : "");
        }
        print "";
    }
}' >> "$CSV_OUT"
