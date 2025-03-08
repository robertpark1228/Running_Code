 awk '{ for (i=1; i<=NF; i+=2) printf "%s\t", $i; print "" }' EPIC.txt > out.txt
