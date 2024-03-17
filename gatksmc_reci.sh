#납품용_분석용_논문용
1) Trim Galore
${trim_galore} --trim1 --illumina --paired fastq_R_1.fq Fastq_R_2.fq -O [파일명]

2) bwa-mem 
bwa mem -M -t 16 -R ${readgroup} hg38.fa fastq_R_1.fq Fastq_R_2.fq | samtools view -@ 16 -Sb | samtools sort -@ 16 -O bam -o ${output}_sorted.bam
READGROUP : ID:[파일명] PL:ILLUMINA     LB:Truseq       SM:WGS0XXX PI:350

3) gatk - MarkDuplicate
java -jar -Xmx8g /disk1/oneomics_analysis/beta_pipeline/modules/gatk.jar MarkDuplicates REMOVE_DUPLICATES=true I=[파일명].bam O=[파일명]_sorted.bam_sorted_MarkDuplicate.bam M= [파일명]  _sorted.bam_markduplicate.txt COMPRESSION_LEVEL=1 CREATE_INDEX=true

4) gatk -  BaseRecalibrator
필요 래퍼런스 파일명
#wgs_interval_contig_grch38.list
#hs38.fa
#Homo_sapiens_assembly38.dbsnp138.vcf.gz
#Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#Homo_sapiens_assembly38.known_indels.vcf.gz
#1000G_phase1.snps.high_confidence.hg38.vcf.gz

java -jar -Xmx8g /disk1/oneomics_analysis/beta_pipeline/modules/gatk.jar BaseRecalibrator 
-L /disk1/oneomics_analysis/beta_pipeline/references/wgs_interval_contig_grch38.list 
-R /disk1/oneomics_analysis/beta_pipeline/references/hg38/hs38.fa -I  [파일명]_sorted.bam_sorted_MarkDuplicate.bam   
--known-sites /disk1/oneomics_analysis/beta_pipeline/references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz 
--known-sites /disk1/oneomics_analysis/beta_pipeline/references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
--known-sites /disk1/oneomics_analysis/beta_pipeline/references/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz 
--known-sites /disk1/oneomics_analysis/beta_pipeline/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz 
-O  [파일명]_sorted.bam_sorted_MarkDuplicate.bam.table

5)  gatk -  ApplyBQSR
  필요 래퍼런스 파일명  
#wgs_interval_contig_grch38.list
#hs38.fa

java -jar -Xmx8g /disk1/oneomics_analysis/beta_pipeline/modules/gatk.jar ApplyBQSR 
-L /disk1/oneomics_analysis/beta_pipeline/references/wgs_interval_contig_grch38.list 
-R /disk1/oneomics_analysis/beta_pipeline/references/hg38/hs38.fa 
-I  [파일명]_sorted.bam_sorted_MarkDuplicate.bam   
-bqsr [파일명]_sorted.bam_sorted_MarkDuplicate.bam.table 
 -O  [파일명]  _sorted_MarkDuplicate_bqsr.bam

6) gatk -  HaplotypeCaller
  필요 래퍼런스 파일명  
#wgs_interval_contig_grch38.list
#hs38.fa

java -jar -Xmx8g /disk1/oneomics_analysis/beta_pipeline/modules/gatk.jar HaplotypeCaller 
-L /disk1/oneomics_analysis/beta_pipeline/references/wgs_interval_contig_grch38.list 
-R /disk1/oneomics_analysis/beta_pipeline/references/hg38/hs38.fa 
-I [파일명]  _sorted_MarkDuplicate_bqsr.bam  
-O [파일명]  _sorted_MarkDuplicate_bqsr.bam.g.vcf --native-pair-hmm-threads 32 --ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation

7)  gatk -  GenotypeGVCFs 
  #hs38.fa  
java -jar -Xmx8g /disk1/oneomics_analysis/beta_pipeline/modules/gatk.jar GenotypeGVCFs 
-R /disk1/oneomics_analysis/beta_pipeline/references/hg38/hs38.fa 
-V [파일명]  _sorted_MarkDuplicate_bqsr.bam.vcf
