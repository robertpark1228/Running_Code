#!/bin/bash
# merged_rename_vcf_menu.sh
# Revised version: This script renames VCF sample names and updates headers using bcftools.
# In VCF list mode, you can now select multiple list files.
# For each selected list file, its VCF entries are processed concurrently.
# When all VCFs in a list file complete, the script then proceeds to the next list file.
#
# The sample mapping file (name_map) must be a tab-delimited file with columns:
#   old_fam1   old_sample1   new_fam1   new_sample1
#
# Dependencies: bcftools must be installed.
#
# Usage:
#   ./merged_rename_vcf_menu.sh
#

# Ensure bcftools is installed
if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools is not installed. Please install it (e.g., sudo apt install bcftools)."
    exit 1
fi

today_date=$(date +"%Y%m%d")

# Function to perform renaming for a single VCF file
rename_vcf() {
    local input_vcf="$1"
    local sample_file="$2"
    local output_vcf="$3"

    # Check if input VCF exists
    if [ ! -f "$input_vcf" ]; then
        echo "Error: Input VCF file '$input_vcf' not found!"
        return 1
    fi

    # Use the basename of the mapping file (without .txt) for base folder creation.
    local mapping_filename
    mapping_filename=$(basename "$sample_file")
    local mapping_base="${mapping_filename%.txt}"  # e.g. name_map_D24023

    # Base directory is fixed by the mapping file:
    local base_dir="${PWD}/${mapping_base}"

    # Extract CODA_D[XXXXX] from the list file name if available,
    # otherwise fallback to extracting from the output VCF file name.
    local extracted=""
    if [ -n "$LIST_FILE_NAME" ]; then
        extracted=$(echo "$LIST_FILE_NAME" | grep -oE 'CODA_D[0-9]{5}' | head -n1)
    else
        extracted=$(echo "$output_vcf" | grep -oE 'CODA_D[0-9]{5}')
        if [ -z "$extracted" ]; then
            local dnum
            dnum=$(echo "$output_vcf" | grep -oE 'D[0-9]{5}')
            if [ -n "$dnum" ]; then
                extracted="CODA_${dnum}"
            else
                extracted="CODA_${today_date}"
            fi
        fi
    fi

    # Create the working and final directories.
    # Final output is placed under:
    # $(pwd)/name_map_D24023/CODA_D[XXXXX]/CODA_Korea_Biobank_Array_VCF/
    local work_dir="${base_dir}/work"
    local final_dir="${base_dir}/${extracted}/CODA_Korea_Biobank_Array_VCF"
    mkdir -p "$work_dir"
    mkdir -p "$final_dir"

    # Create mapping file: extract columns 2 and 4 from sample_file; save it in work directory.
    local mapping_file="${work_dir}/${mapping_filename}.vcf.txt"
    if [ ! -f "$mapping_file" ]; then
        echo "Mapping file not found. Creating mapping file from $sample_file ..."
        cut -f2,4 "$sample_file" > "$mapping_file"
        sleep 3
    else
        echo "Using existing mapping file: $mapping_file"
        sleep 2
    fi

    echo "Renaming samples in VCF: $input_vcf"
    # Step 1: Rename samples in the VCF using bcftools reheader (output goes to work directory)
    bcftools reheader -s "$mapping_file" -o "${work_dir}/${output_vcf}" "$input_vcf" --threads 4

    # Step 2: Extract header and modify it: remove any bcftools-specific header lines,
    # then insert a revision line at line 4 using the extracted pattern and current date.
    bcftools view -h "${work_dir}/${output_vcf}" > "${work_dir}/${output_vcf}.header.txt"
    cat "${work_dir}/${output_vcf}.header.txt" | grep -v "##bcftools" | sed "4i ##REVISION=${extracted}_${today_date}" > "${work_dir}/${output_vcf}_changed_header.txt"

    # Step 3: Reheader the VCF with the updated header, writing the final file to final_dir.
    bcftools reheader -h "${work_dir}/${output_vcf}_changed_header.txt" -o "${final_dir}/${output_vcf}" "${work_dir}/${output_vcf}"

    # Index the output if it is compressed (.gz)
    if [[ "${final_dir}/${output_vcf}" == *.gz ]]; then
        echo "Indexing output VCF: ${final_dir}/${output_vcf}"
        bcftools index "${final_dir}/${output_vcf}"
    fi

    echo "Sample renaming complete. Output saved to: ${final_dir}/${output_vcf}"
}

# Function for processing a single VCF file
process_single() {
    echo "---- Single VCF Mode ----"
    read -p "Enter sample mapping file (name_map) path: " sample_file
    if [ ! -f "$sample_file" ]; then
        echo "Error: Sample file '$sample_file' not found!"
        exit 1
    fi

    while true; do
        read -p "Enter input VCF file path: " input_vcf
        if [ ! -f "$input_vcf" ]; then
            echo "Error: Input VCF file '$input_vcf' not found!"
            continue
        fi

        read -p "Enter output VCF file name (e.g., output.vcf.gz): " output_vcf
        rename_vcf "$input_vcf" "$sample_file" "$output_vcf"

        read -p "Do you want to process another single VCF file? (y/n): " answer
        if [[ "$answer" != "y" && "$answer" != "Y" ]]; then
            break
        fi
    done
}

# Function for processing one or more VCF list files concurrently per file
process_list() {
    echo "---- VCF List Mode ----"
    read -p "Enter sample mapping file (name_map) path: " sample_file
    if [ ! -f "$sample_file" ]; then
        echo "Error: Sample file '$sample_file' not found!"
        exit 1
    fi

    # Display default list file options
    echo "Select VCF list file(s) (you can choose multiple, e.g., 1,2,3):"
    echo "  1) CODA_D24023_KBA_KA_5493_VCF.txt"
    echo "  2) CODA_D24033_KBA_KN_8105_VCF.txt"
    echo "  3) CODA_D24038_KBA_KC_58693_VCF.txt"
    echo "  4) Enter custom file path"
    read -p "Enter your choice(s): " list_choices_raw

    # Convert comma-separated values to space separated tokens
    list_choices=$(echo "$list_choices_raw" | tr ',' ' ')

    for choice in $list_choices; do
        case "$choice" in
            1)
                list_file="./metadata/CODA_D24023_KBA_KA_5493_VCF.txt"
                ;;
            2)
                list_file="./metadata/CODA_D24033_KBA_KN_8105_VCF.txt"
                ;;
            3)
                list_file="./metadata/CODA_D24038_KBA_KC_58693_VCF.txt"
                ;;
            4)
                read -p "Enter VCF list file path: " list_file
                ;;
            *)
                echo "Invalid choice: $choice. Skipping."
                continue
                ;;
        esac

        # If the list file does not exist, auto-generate a default template.
        if [ ! -f "$list_file" ]; then
            echo "VCF list file '$list_file' not found. Auto-generating default list file."
            {
                echo "# Format: <input_vcf_path> <output_vcf_file_name>"
                echo "/path/to/input.vcf output.vcf.gz"
            } > "$list_file"
            echo "Default template created at '$list_file'. Please edit this file with your VCF paths."
        fi

        export LIST_FILE_NAME="$list_file"
        TOTAL_LINES=$(grep -cvE '^\s*$|^#' "$list_file")
        if [ "$TOTAL_LINES" -eq 0 ]; then
            echo "Error: No valid entries found in $list_file!"
            continue
        fi

        CURRENT_LINE=0
        echo "Processing VCF list file: $list_file"
        while read -r line || [ -n "$line" ]; do
            # Skip empty lines and comment lines starting with #
            if [[ -z "$line" || "$line" =~ ^# ]]; then
                continue
            fi

            # Extract input VCF (first column) and output VCF (second column)
            input_vcf_line=$(echo "$line" | awk '{print $1}')
            output_vcf_line=$(echo "$line" | awk '{print $2}')

            if [[ -n "$input_vcf_line" && -n "$output_vcf_line" ]]; then
                # Process each VCF concurrently (background execution)
                rename_vcf "$input_vcf_line" "$sample_file" "$output_vcf_line" &
            else
                echo "Warning: Skipping invalid line: $line"
            fi

            ((CURRENT_LINE++))
            PERCENT=$((CURRENT_LINE * 100 / TOTAL_LINES))
            printf "\rProgress for this list: [%-50s] %d%%" "$(printf '#%.0s' $(seq 1 $((PERCENT / 2))))" "$PERCENT"
        done < "$list_file"

        # Wait for all background processes (for this list file) to finish
        wait
        echo -e "\nCompleted processing VCF list file: $list_file"
    done
}

# Main interactive menu
echo "=============================="
echo " VCF Sample Renaming Tool Menu"
echo "=============================="
echo "Select processing mode:"
echo "  1) Single VCF file"
echo "  2) VCF list file (multiple list files can be selected)"
read -p "Enter choice (1 or 2): " mode_choice

case "$mode_choice" in
    1)
        process_single
        ;;
    2)
        process_list
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo "Processing complete."
