#!/bin/bash

# Usage: ./vep_full_dbnsfp_annotation.sh input.vcf output.tsv

INPUT_VCF=$1
OUTPUT_TSV=$2

# Paths (update these!)
DBNSFP_PATH="/disk1/oneomics_analysis/references/dbNSFP5.1a/dbNSFP5.1a_all.tsv.gz"
FASTA_PATH="/disk1/oneomics_analysis/references/vep/Homo_sapiens.GRCh38.dna.toplevel.fa"
PLUGIN_DIR="/disk1/oneomics_analysis/references/vep/VEP_plugins"

# Run VEP with full dbNSFP and expanded field set
vep \
  -i "$INPUT_VCF" \
  -o "$OUTPUT_TSV" \
  --cache \
  --offline \
  --fasta "$FASTA_PATH" \
  --assembly GRCh38 \
  --tab \
  --fork 80 \
  --offline \
  --dir_cache /disk1/oneomics_analysis/references/vep/ \
  --everything
  --force_overwrite \
  --dir_plugins "$PLUGIN_DIR" \
  --plugin dbNSFP,"$DBNSFP_PATH",\
  SIFT_score,SIFT_converted_rankscore,SIFT_pred,\
  LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,\
  MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,\
  MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,\
  FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,\
  PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,\
  MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,\
  MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,\
  M-CAP_score,M-CAP_rankscore,M-CAP_pred,\
  MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,\
  fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,\
  Eigen-raw_coding,Eigen-phred_coding,Eigen-PC-raw_coding,Eigen-PC-phred_coding,Eigen-PC-raw_coding_rankscore,\
  integrated_fitCons_score,integrated_fitCons_rankscore,integrated_confidence_value,\
  GERP++_NR,GERP++_RS,GERP++_RS_rankscore,\
  gnomAD_exomes_AF,gnomAD_exomes_AFR_AF,gnomAD_exomes_AMR_AF,gnomAD_exomes_EAS_AF,gnomAD_exomes_FIN_AF,gnomAD_exomes_NFE_AF,gnomAD_exomes_SAS_AF,\
  gnomAD_genomes_AF,gnomAD_genomes_AFR_AF,gnomAD_genomes_AMR_AF,gnomAD_genomes_EAS_AF,gnomAD_genomes_FIN_AF,gnomAD_genomes_NFE_AF,\
  GDI,GDI-Phred,\
  Gene_damage_prediction(all_disease-causing_genes),Gene_damage_prediction(all_Mendelian_disease-causing_genes),\
  Gene_damage_prediction(Mendelian_AD_disease-causing_genes),Gene_damage_prediction(Mendelian_AR_disease-causing_genes),\
  Gene_damage_prediction(all_PID_disease-causing_genes),Gene_damage_prediction(PID_AD_disease-causing_genes),\
  Gene_damage_prediction(PID_AR_disease-causing_genes),Gene_damage_prediction(all_cancer_disease-causing_genes),\
  Gene_damage_prediction(cancer_recessive_disease-causing_genes),Gene_damage_prediction(cancer_dominant_disease-causing_genes"\
  
  --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,IMPACT,HGVSc,HGVSp,Existing_variation,CLINVAR_CLNSIG,CLINVAR_CLNDN,CLINVAR_CLNREVSTAT,REF_AA,ALT_AA,SIFT_score,SIFT_pred,PolyPhen_score,PolyPhen_pred,REVEL_score,CADD_phred,MetaSVM_pred,MetaLR_pred,GDI,GDI-Phred,Gene_damage_prediction(all_disease-causing_genes),Gene_damage_prediction(Mendelian_AR_disease-causing_genes),gnomAD_exomes_AF,gnomAD_genomes_AF,gnomAD_exomes_AFR_AF,gnomAD_exomes_EAS_AF,gnomAD_genomes_EAS_AF,GO_biological_process,GO_cellular_component,GO_molecular_function,GTEx_V8_tissue,GTEx_V8_gene,Pathway(KEGG)_id,Pathway(KEGG)_full,Pathway(Uniprot),Trait_association(GWAS),MIM_id,MIM_disease,Function_description,Disease_description,Expression(egenetics),Expression(GNF/Atlas),Interactions(IntAct),Interactions(BioGRID),Interactions(ConsensusPathDB)"

echo -e "\n✅ Annotation complete.\nOutput written to: $OUTPUT_TSV"
