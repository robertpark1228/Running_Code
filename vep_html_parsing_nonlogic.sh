#!/bin/bash
file=${1}
awk '/General statistics/,/Variant classes/' ${file} > 1
awk '/Variant classes/,/Consequences \(most severe\)/' ${file} > 2
awk '/Consequences \(most severe\)/,/Consequences \(all\)/' ${file} > 3
awk '/Consequences \(all\)/,/Coding consequences/' ${file} > 4
awk '/Coding consequences/,/Variants by chromosome/' ${file} > 5

while read line;do cat 1 | grep -P "${line}\t";done < cata1 | sed '{s/\t/ /g}' | sed '{s/  */ /g}' | uniq > ${file}_1.txt
while read line;do cat 2 | grep -P "${line}\t";done < cata2 | sed '{s/\t/ /g}' | sed '{s/  */ /g}' | uniq > ${file}_2.txt
while read line;do cat 3 | grep -P "${line}\t";done < cata3 | sed '{s/\t/ /g}' | sed '{s/  */ /g}' | uniq > ${file}_3.txt
while read line;do cat 4 | grep -P "${line}\t";done < cata4 | sed '{s/\t/ /g}' | sed '{s/  */ /g}' | uniq > ${file}_4.txt
while read line;do cat 5 | grep -P "${line}\t";done < cata5 | sed '{s/\t/ /g}' | sed '{s/  */ /g}' | uniq > ${file}_5.txt


while read line;do cat ${file}_1.txt | grep -q "${line}" && grep "${line}" ${file}_1.txt  || echo ;done < cata1
while read line;do cat ${file}_2.txt | grep -q "${line}" && grep "${line}" ${file}_2.txt  || echo ;done < cata2
while read line;do cat ${file}_3.txt | grep -q "${line}" && grep "${line}" ${file}_3.txt  || echo ;done < cata3
while read line;do cat ${file}_4.txt | grep -q "${line}" && grep "${line}" ${file}_4.txt  || echo ;done < cata4
while read line;do cat ${file}_5.txt | grep -q "${line}" && grep "${line}" ${file}_5.txt  || echo ;done < cata5
