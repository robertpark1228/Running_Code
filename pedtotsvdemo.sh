#!/bin/bash

cat merged_imputed_exome_plink.ped | cut -d " " -f1-6 > ID_COL
cat ID_COL | cut -d " " -f2 > ID
(echo ID;cat ID) > ID_HEADER
cat merged_imputed_exome_plink.ped | cut -d " " -f7- > GENO
#분기
cat GENO | awk -F " " '{for(i=1; i<=(NF / 2); i++) printf ("%s ", $(2i-1)"/"$(2i));print " " }' > GENO_OUT
cat merged_imputed_exome_plink.map | cut -f2 | tr "\n" "\t" | sed '{s/\t$/\n/g}' > HEADER
cat HEADER GENO_OUT > HEADER_GENO_OUT
paste ID_HEADER HEADER_GENO_OUT >  merged_imputed_exome_plink.tsv
