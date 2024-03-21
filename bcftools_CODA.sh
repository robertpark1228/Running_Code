bcftools view -I -s '' chr3.vcf.gz --force-samples --threads 60 -O z -o chr3.only.vcf.gz

bcftools filter --threads 40 -i'INFO/AF>0.01 & INFO/AF<0.05' KoGES5000_chrY_vqsr_99.9_PASS_VEP.vcf.gz | grep -v "^#"
