for i in $(cat id);do printf "${i}\n%.0s" {1..2};done

scp -rvP 3030 ./WGS_*_output
