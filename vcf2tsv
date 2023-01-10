#!/bin/bash
#filelocation
input=${1}
#name=$(echo ${input} | sed '{s/\.vcf\.gz//g}')
#uniqID=$(date "+%Y%m%d_%H%M")
#echo ${input}

help()
{
   echo "VCF2tsv/LowerFreq."
   echo
   echo "Syntax: [-h|v]"
   echo "options:"
   echo "-h     Print this Help."
   echo "-v     use with [bgzipped VCF]"
   echo



}

inputvcf()
{
echo "---------------------------------------------------------------------------------------------------------------------"
echo "Convertor runs with file : ${input}"
uniqID=$(date "+%Y%m%d_%H%M%S")
name=$(echo ${input} | sed '{s/\.vcf\.gz//g}')
echo "Result file name will be: ${name}.${uniqID}.patho.lower0.01.tsv"
echo "---------------------------------------------------------------------------------------------------------------------"
echo "----------------------------------------Ignore Warning Message-------------------------------------------------------"
#tabix
#if tbi exist : skip with warning message
echo "TABIX GENERATE"
tabix -p vcf ${input}

echo "GRCh38 KRGDB Database Annotation..."
#KRGDB Annotation
bcftools annotate -a /disk1/references/KRGDB_GRCh38.tabixed.vcf.gz -c INFO/KRGDB_ALF ${input} --threads 40 -O z -o ${name}.krgdb.vcf.gz


#Header Part
echo "Writing header..."
printf "Gene\tCoordinate&Variant\tGenotype\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tAllele_Freq\tgnomAD_EAS_AF\tClinVar Sifnificance KRGDB\n" > ${name}.${uniqID}.patho.lower0.01.tsv

#Report Part with uniqID
echo "Filling file with varaitns information..."
bcftools +split-vep ${name}.krgdb.vcf.gz -f '%SYMBOL %CHROM:%POS"_"%REF/%ALT [%GT] %Consequence %SIFT %PolyPhen %AF %gnomAD_EAS_AF %CLIN_SIG %KRGDB_ALF %ANN\n' -d -A space | awk -F "|" '{ print $1,$10,$11 }' | awk -F " " '{print $1,$2,$3,$4,$5,$6,$12,$13,$7,$8,$9,$10}' | sed '{s/ /\t/g}' | grep "patho" | awk -F " " '{if($12<0.01)print $0}' | sed '{s/1\/1/hom/g}' | sed '{s/0\/1/het/g}' | sed '{s/ /\t/g}' | sed '{s/"//g}' >> ${name}.${uniqID}.patho.lower0.01.tsv


echo "----------------------------------------------------------------------------------------------------------------------"

echo "DONE!"
echo "Result:${name}.${uniqID}.patho.lower0.01.tsv"
}


while getopts "v:h" option; do
   case $option in
      h) # display Help
         help
         exit;;
      v) # incorrect option
         input=${OPTARG}
         inputvcf
         exit;;
   esac
done

echo "help:-h"
echo "example bash convert.sh -v [file:vcf.gz]"

