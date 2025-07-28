while read line;do a=$(echo $line | cut -d " " -f1);b=$(echo $line | cut -d " " -f2);bash id.sh ${a} ${b};done < CCDC63 > CCDC63.json
