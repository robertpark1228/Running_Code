#!/bin/bash
input=${1}

vep -i ${input} -o ${input}.vep.dbnsfp.vcf --vcf --cache --offline --assembly GRCh38 --plugin dbNSFP,/disk1/oneomics_analysis/references/dbNSFP5.1a/dbNSFP5.1a_all.tsv.gz,ALL --dir_plugins /disk1/oneomics_analysis/references/vep/VEP_plugins --fork 80 --offline --dir_cache /disk1/oneomics_analysis/references/vep/ --everything --force_overwrite --fasta /disk1/oneomics_analysis/references/vep/Homo_sapiens.GRCh38.dna.toplevel.fa --hgvs --variant_class --sift b --polyphen b --humdiv --gene_phenotype --hgvsg


vcf="${input}.vep.dbnsfp.vcf"
out="${input}.vep.dbnsfp.vcf.tsv"

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



csvtk join -f Gene,Gene ${out} dbNSFP5.1_gene.tsv -t -H -k > ${out}_with_gene.tsv
#grep -v "^##" ${input}.vep.tsv > ${input}.vep.tsv.temp

#csvtk join -f Gene,Gene ${input}.vep.tsv.temp dbNSFP5.1_gene.tsv -t -H -k > ${input}.vep.tsv.txt
