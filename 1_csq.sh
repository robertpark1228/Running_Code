#!/bin/bash

# Usage: ./extract_all_csq_entries.sh input.vep.vcf output.tsv

vcf="$1"
out="$2"

# Extract CSQ format string from VCF
csq_header=$(grep -oP '^##INFO=<ID=CSQ.*?Format: \K[^"]+' "$vcf")
csq_header_tab=$(echo "$csq_header" | tr "|" "\t")
csq_field_count=$(echo "$csq_header" | awk -F"|" '{print NF}')

# Output header
echo -e "CHROM\tPOS\tREF\tALT\tDP\tAD\tGT\tMQ\tQUAL\tFILTER\tZygosity\t$csq_header_tab" > "$out"

# Extract CSQ with multiple transcript annotations
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD\t%GT\t%MQ\t%QUAL\t%FILTER\t%CSQ\n]' "$vcf" \
| awk -F'\t' -v OFS='\t' -v csq_count="$csq_field_count" '
{
    z = ($7=="0/1"||$7=="1/0"||$7=="0|1"||$7=="1|0") ? "Heterozygous" :
        ($7=="1/1"||$7=="1|1") ? "Homozygous" : "Unknown";

    # CSQ field may contain multiple comma-separated annotations
    split($11, all_csq, ",");
    for (j in all_csq) {
        split(all_csq[j], csq, "|");
        for (i = 1; i <= csq_count; i++) {
            if (!(i in csq)) csq[i] = ".";
        }
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,z;
        for (i = 1; i <= csq_count; i++) {
            printf "\t%s", csq[i];
        }
        printf "\n";
    }
}' >> "$out"
