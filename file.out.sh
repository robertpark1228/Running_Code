for i in $(find -name "as5*.tsv");do (head -n 1 ${i} && grep -F -f samples ${i}) > ${i}_filtered.tsv;done
