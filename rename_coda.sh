#!/bin/bash

# Usage function with example
usage() {
    echo "Usage: $0 <input_vcf> <name_map> <output_vcf>"
    echo "Example: $0 input.vcf.gz name_map.txt output.vcf.gz"
    echo
    echo "Expected format of name_map.txt (tab-separated):"
    echo "-------------------------------------"
    echo "old_sample1    new_sample1"
    echo "old_sample2    new_sample2"
    echo "old_sample3    new_sample3"
    echo "-------------------------------------"
    exit 1
}

# Check for three arguments
if [ "$#" -ne 3 ]; then
    usage
fi

# Assign arguments to variables
INPUT_VCF="$1"
NAME_MAP="$2"
OUTPUT_VCF="$3"

# Check if input VCF exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF file '$INPUT_VCF' not found!"
    exit 1
fi

# Check if name mapping file exists
if [ ! -f "$NAME_MAP" ]; then
    echo "Error: Name mapping file '$NAME_MAP' not found!"
    exit 1
fi


# Check if bcftools is installed
if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools is not installed. Install it with: sudo apt install bcftools"
    exit 1
fi

# Check if name mapping_2 file exists
if [ ! -f ${PWD}/${NAME_MAP%.txt}/work/$NAME_MAP.vcf.txt ]
	then
	echo "MAPPING_VCF FILE IS NOT EXIST"
	echo "sleep 5sec"
	sleep 5
	mkdir -p ${PWD}/${NAME_MAP%.txt}
	mkdir -p ${PWD}/${NAME_MAP%.txt}/work
        cat $NAME_MAP | cut -f2,4 > ${PWD}/${NAME_MAP%.txt}/work/$NAME_MAP.vcf.txt
	
	        echo "check ${PWD}/${NAME_MAP%.txt}/work/$NAME_MAP.vcf.txt is existed"
        sleep 3
	echo "echo ${INPUT_VCF}" 
        echo "Renaming samples in VCF..."
	echo

##SAMPLENAMECHANGE
        bcftools reheader -s ${PWD}/${NAME_MAP%.txt}/work/$NAME_MAP.vcf.txt -o ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF} ${INPUT_VCF} --threads 4

##HEADER_CUT
        bcftools view -h ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF} > ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}.header.txt
        
##HEADER_CHANGE
	cat ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}.header.txt | grep -v "##bcftools" | sed "4i ##REVISION=CODA250228" > ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}_changed_header.txt
        
##REHEADER_WITH_NEWID_CODA
	bcftools reheader -h ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}_changed_header.txt -o ${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF} ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}


        if [[ "${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF}" == *.gz ]]
                then
                        echo "Indexing output VCF..."
                        bcftools index ${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF}
        fi

                echo "Sample renaming complete. Output saved to: $OUTPUT_VCF"
	
	

	else

	echo "check ${PWD}/${NAME_MAP%.txt}/work/$NAME_MAP.vcf.txt is existed"
	sleep 3
	echo WORK_VCF: "${INPUT_VCF}"
	echo "Renaming samples in VCF..."
##SAMPLENAMECHANGE
echo "MAP_NAME_ORI2CODA"
        bcftools reheader -s ${PWD}/${NAME_MAP%.txt}/work/$NAME_MAP.vcf.txt -o ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF} ${INPUT_VCF} --threads 4

##HEADER_CUT
echo "TAKING_HEADER"
        bcftools view -h ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF} > ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}.header.txt

##HEADER_CHANGE
echo "DELETE_INFO_HEADER_ADD_CODA"
        cat ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}.header.txt | grep -v "##bcftools" | sed "4i ##REVISION=CODA250228" > ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}_changed_header.txt

##REHEADER_WITH_NEWID_CODA
echo "REHEADER"
        bcftools reheader -h ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}_changed_header.txt -o ${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF} ${PWD}/${NAME_MAP%.txt}/work/${OUTPUT_VCF}








	if [[ "${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF}" == *.gz ]]
		then
    			echo "Indexing output VCF..."
    			bcftools index ${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF}
	fi

		echo "Sample renaming complete. Output saved to: $OUTPUT_VCF"


fi



# Rename samples in the VCF file
#echo "Renaming samples in VCF..."
#cat $NAME_MAP | cut -f2,4 > $NAME_MAP.vcf.txt
#bcftools reheader -s $NAME_MAP.vcf.txt -o ${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF} ${INPUT_VCF}

# Index the new VCF file if it's compressed
#if [[ "${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF}" == *.gz ]]; then
#    echo "Indexing output VCF..."
#    bcftools index ${PWD}/${NAME_MAP%.txt}/${OUTPUT_VCF}
#fi

#echo "Sample renaming complete. Output saved to: $OUTPUT_VCF"
