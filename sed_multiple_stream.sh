while read line;do n=$(echo ${line} | cut -d "," -f1);c=$(echo ${line} | cut -d "," -f2);echo sed \"{s/${n}$/${c}/g}\";done < zout_genus.csv | sed -z 's/\n/ | /g'
