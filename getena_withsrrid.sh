#!/bin/bash
input=${1}
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${input}&result=read_run&fields=fastq_ftp" \
	| tail -n +2 \
	| cut -f2 \
	| tr ';' '\n' \
	| sed 's|^|https://|'

