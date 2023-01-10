input=${1}
date_folder=$(date "+%Y_%m_%d")
mkdir -p ./${date_folder}
while read line;do name=$(echo $line | cut -d " " -f1 | sed '{s/\.\///g}');f1=$(echo $line | cut -d " " -f2);f2=$(echo $line | cut -d " " -f3);mkdir ./${date_folder}/${name}; hlahd.sh -t 60 -m 150 -c 0.7 -f /disk1/oneomics_analysis/HLA_Genotype/freq_data/ ${f1} ${f2} /disk1/oneomics_analysis/HLA_Genotype/HLA_gene.split.manipulate.txt /disk1/oneomics_analysis/HLA_Genotype/dictionary $name ./${date_folder}/${name};done < ${input}
