vep \
  -i input.vcf \
  -o output_vep_dbnsfp.txt \
  --tab \
  --cache \
  --offline \
  --assembly GRCh38 \
  --dir_plugins ~/.vep/Plugins \
  --plugin dbNSFP,/path/to/dbNSFP5.1a_all.tsv.gz,SIFT_score,Polyphen2_HDIV_pred,REVEL_score,CADD_phred \
  --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,HGVSc,HGVSp,SIFT_score,Polyphen2_HDIV_pred,REVEL_score,CADD_phred"
