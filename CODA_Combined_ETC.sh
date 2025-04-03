#!/bin/bash
# combined_processing.sh
# Combined processing script for CODA datasets.
#
# Usage: ./combined_processing.sh <type> <metadata_file>
#
# If you run the script with the --help or -h option, it will display this help message.
#
# Parameters:
#   <type>: Specifies which dataset to process. The following full file names are supported:
#
#     CODA_D24023_CGH_Array_chip_CNV_TXT_1000
#       - Explanation: Project "CODA_D24023", data type "CGH_Array_chip_CNV_TXT", 1000 samples.
#
#     CODA_D24023_Illumina_Hiseq_Exome_Fastq_100
#       - Explanation: Project "CODA_D24023", data type "Illumina_Hiseq_Exome_Fastq", 100 samples.
#
#     CODA_D24023_Infinium_Methylation_450K_CSV_50
#       - Explanation: Project "CODA_D24023", data type "Infinium_Methylation_450K_CSV", 50 samples.
#
#     CODA_D24023_Infinium_Methylation_450K_IDAT_4
#       - Explanation: Project "CODA_D24023", data type "Infinium_Methylation_450K_IDAT", 4 samples.
#
#     CODA_D24027_Infinium_Methylation_450K_CSV_50
#       - Explanation: Project "CODA_D24027", data type "Infinium_Methylation_450K_CSV", 50 samples.
#
#     CODA_D24027_Infinium_Methylation_450K_IDAT_400
#       - Explanation: Project "CODA_D24027", data type "Infinium_Methylation_450K_IDAT", 400 samples.
#
#     CODA_D24027_Infinium_Methylation_850K_IDAT_1528
#       - Explanation: Project "CODA_D24027", data type "Infinium_Methylation_850K_IDAT", 1528 samples.
#
#     CODA_D24038_HEXA_methylation850K_IDAT_822
#       - Explanation: Project "CODA_D24038", data type "HEXA_methylation850K_IDAT", 822 samples.
#
#     CODA_D24038_Infinium_Methylation_850K_IDAT_822
#       - Explanation: Project "CODA_D24038", data type "Infinium_Methylation_850K_IDAT", 822 samples.
#
#   <metadata_file>: Path to the metadata text (TXT) file containing sample mapping information.
#                    This file must be tab-separated with the columns:
#                    ORIGINAL_FAM_ID, ORIGINAL_SAMPLE_ID, CONVERTED_FAM_ID, CONVERTED_SAMPLE_ID.
#
# Interactive Mode:
#   If you do not supply both arguments, you will be prompted with a numbered menu.
#   You may enter a comma-separated list (e.g., "1,3,5") to select multiple dataset types.
#   Then you will be asked for one metadata file which is used for all selections.
#   If the metadata file does not exist, an error is printed and the script exits.
#
# Example:
#   ./combined_processing.sh CODA_D24023_CGH_Array_chip_CNV_TXT_1000 sample_metadata.txt
#

# Display usage information.
show_usage() {
    echo "Usage: $0 <type> <metadata_file>"
    echo ""
    echo "Parameters:"
    echo "  <type>         : The dataset type to process. Choose one of:"
    echo "                   CODA_D24023_CGH_Array_chip_CNV_TXT_1000"
    echo "                     (Project: CODA_D24023, Data: CGH_Array_chip_CNV_TXT, 1000 samples)"
    echo "                   CODA_D24023_Illumina_Hiseq_Exome_Fastq_100"
    echo "                     (Project: CODA_D24023, Data: Illumina_Hiseq_Exome_Fastq, 100 samples)"
    echo "                   CODA_D24023_Infinium_Methylation_450K_CSV_50"
    echo "                     (Project: CODA_D24023, Data: Infinium_Methylation_450K_CSV, 50 samples)"
    echo "                   CODA_D24023_Infinium_Methylation_450K_IDAT_4"
    echo "                     (Project: CODA_D24023, Data: Infinium_Methylation_450K_IDAT, 4 samples)"
    echo "                   CODA_D24027_Infinium_Methylation_450K_CSV_50"
    echo "                     (Project: CODA_D24027, Data: Infinium_Methylation_450K_CSV, 50 samples)"
    echo "                   CODA_D24027_Infinium_Methylation_450K_IDAT_400"
    echo "                     (Project: CODA_D24027, Data: Infinium_Methylation_450K_IDAT, 400 samples)"
    echo "                   CODA_D24027_Infinium_Methylation_850K_IDAT_1528"
    echo "                     (Project: CODA_D24027, Data: Infinium_Methylation_850K_IDAT, 1528 samples)"
    echo "                   CODA_D24038_HEXA_methylation850K_IDAT_822"
    echo "                     (Project: CODA_D24038, Data: HEXA_methylation850K_IDAT, 822 samples)"
    echo "                   CODA_D24038_Infinium_Methylation_850K_IDAT_822"
    echo "                     (Project: CODA_D24038, Data: Infinium_Methylation_850K_IDAT, 822 samples)"
    echo ""
    echo "  <metadata_file>: Path to the metadata TXT file containing sample mapping information."
    echo "                   (Must be tab-separated with header columns: ORIGINAL_FAM_ID, ORIGINAL_SAMPLE_ID,"
    echo "                    CONVERTED_FAM_ID, CONVERTED_SAMPLE_ID)"
    echo ""
    echo "Interactive Mode:"
    echo "  If you do not supply both arguments, you will be prompted with a numbered menu."
    echo "  You may enter a comma-separated list (e.g., '1,3,5') to select multiple dataset types."
    echo "  Then you will be asked for a single metadata file which is used for all selections."
    exit 1
}

# Check for help option.
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_usage
fi

# Declare associative array for metadata file(s) and selected types.
declare -A type_input_files
selected_types=()

if [[ $# -lt 2 ]]; then
    echo "Not enough arguments provided. Running in interactive mode."
    echo ""
    echo "Please select one or more dataset types (enter numbers separated by commas):"
    options=(
      "CODA_D24023_CGH_Array_chip_CNV_TXT_1000"
      "CODA_D24023_Illumina_Hiseq_Exome_Fastq_100"
      "CODA_D24023_Infinium_Methylation_450K_CSV_50"
      "CODA_D24023_Infinium_Methylation_450K_IDAT_4"
      "CODA_D24027_Infinium_Methylation_450K_CSV_50"
      "CODA_D24027_Infinium_Methylation_450K_IDAT_400"
      "CODA_D24027_Infinium_Methylation_850K_IDAT_1528"
      "CODA_D24038_HEXA_methylation850K_IDAT_822"
      "CODA_D24038_Infinium_Methylation_850K_IDAT_822"
      "Quit"
    )
    for i in "${!options[@]}"; do
        printf "  %d) %s\n" "$((i+1))" "${options[$i]}"
    done

    read -p "Enter choices (e.g., 1,3,5): " choice_input
    IFS=',' read -ra choices <<< "$(echo "$choice_input" | tr -d ' ')"
    for choice in "${choices[@]}"; do
        case "$choice" in
            1) selected_types+=("D24023_CGH") ;;
            2) selected_types+=("D24023_Exome") ;;
            3) selected_types+=("D24023_450K_CSV") ;;
            4) selected_types+=("D24023_450K_IDAT") ;;
            5) selected_types+=("D24027_450K_CSV") ;;
            6) selected_types+=("D24027_450K_IDAT") ;;
            7) selected_types+=("D24027_850K_IDAT") ;;
            8) selected_types+=("D24038_HEXA") ;;
            9) selected_types+=("D24038_Infinium") ;;
            10) echo "Exiting."; exit 0 ;;
            *) echo "Invalid choice: $choice. Skipping." ;;
        esac
    done

    if [[ ${#selected_types[@]} -eq 0 ]]; then
        echo "No valid dataset type selected. Exiting."
        exit 1
    fi

    read -p "Enter metadata file path (must be a TXT file): " metadata_path
    if [[ ! -f "$metadata_path" ]]; then
        echo "Error: Metadata file '$metadata_path' does not exist. Exiting."
        exit 1
    fi
    for dtype in "${selected_types[@]}"; do
        type_input_files["$dtype"]="$metadata_path"
    done
else
    selected_types+=("$1")
    if [[ ! -f "$2" ]]; then
        echo "Error: Metadata file '$2' does not exist. Exiting."
        exit 1
    fi
    type_input_files["$1"]="$2"
fi

##################################
# Processing Functions
##################################
process_D24023_CGH() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24023/CGH_Array_chip_CNV_txt_1000/"
    local ref_metadata="/ANALYSIS/ycpark/masking_script/metadata/CGH_ARRAY_SAMPLE_LIST.txt"
    if [[ ! -f "$metadata" ]]; then
        echo "Error: Metadata file '$metadata' not found."
        return 1
    fi
    if [[ ! -f "$ref_metadata" ]]; then
        echo "Error: Reference metadata file '$ref_metadata' not found."
        return 1
    fi
    mkdir -p "$output_dir"
    grep -F -f "$ref_metadata" "$metadata" | cut -f1,4 | while read -r line; do
        local a b
        a=$(echo "$line" | cut -f1)
        b=$(echo "$line" | cut -f2)
        local found_file
        found_file=$(find /ARCHIVE/CODA/CODA_D24023/GENETICS/CODA_CGH_ARRAY_Chip_CNV_TXT_1000/ -name "*${a}_*.txt" 2>/dev/null | head -n 1)
        if [[ -z "$found_file" ]]; then
            echo "Warning: No file found matching '*${a}_*.txt'. Skipping..."
            continue
        fi
        local new_filename="${b}_segMNT.txt"
        cp -v "$found_file" "${output_dir}/${new_filename}"
    done
}

process_D24023_Exome() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24023/CODA_Illumina_Hiseq_Exome_seq_FASTQ_100/"
    local ref_metadata="${PWD}/metadata/CODA_Illumina_Hiseq_Exome_seq_FASTQ_100.list"
    local fastq_dir="/ARCHIVE/CODA/CODA_D24023/GENETICS/CODA_Illumina_Hiseq_Exome_seq_FASTQ_100/"
    if [[ ! -f "$metadata" ]]; then
        echo "Error: Metadata file '$metadata' not found."
        return 1
    fi
    if [[ ! -f "$ref_metadata" ]]; then
        echo "Error: Reference metadata file '$ref_metadata' not found."
        return 1
    fi
    mkdir -p "$output_dir"
    grep -F -f "$ref_metadata" "$metadata" | cut -f1,4 | while IFS=$'\t' read -r a b; do
        a=$(echo "$a" | xargs)
        b=$(echo "$b" | xargs)
        local r1_file="${fastq_dir}${a}_R1.fastq.gz"
        local r2_file="${fastq_dir}${a}_R2.fastq.gz"
        echo "Checking files:"
        echo " - Expected: $r1_file"
        echo " - Expected: $r2_file"
        if [[ -f "$r1_file" && -f "$r2_file" ]]; then
            echo "Copying: $r1_file -> ${output_dir}${b}_R1.fastq.gz"
            echo "Copying: $r2_file -> ${output_dir}${b}_R2.fastq.gz"
            cp -v "$r1_file" "${output_dir}${b}_R1.fastq.gz"
            cp -v "$r2_file" "${output_dir}${b}_R2.fastq.gz"
        else
            echo "Warning: Missing FASTQ files for '${a}', skipping..." >&2
        fi
    done
}

process_D24023_450K_CSV() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24023/OGTT1_5th_methylation_50_1gI_work/"
    local ref_metadata="${PWD}/metadata/OGTT1_5th_methylation_50_1gI.list"
    mkdir -p "$output_dir/CODA_D24023/"
    grep -Ff "$ref_metadata" "$metadata" | cut -f4 | tr '\n' ',' | sed 's/,$/\n/' > "$output_dir/header"
    echo "KoBB,$(cat "$output_dir/header")" > "$output_dir/header.csv"
    tail -n +3 /ANALYSIS/SCRIPT/DATA/CODA_D24023/GENETICS/CODA_Infinium_HumanMethylation_450K_CSV_50/OGTT1_5th_methylation_50_1gI.csv > "$output_dir/${metadata}.body.tmp"
    local output_file="${PWD}/${metadata%.txt}/CODA_D24023/OGTT1_5th_methylation_50_1gI.${metadata%.txt}.masked.csv"
    cat "$output_dir/header.csv" "$output_dir/${metadata}.body.tmp" > "$output_file"
    echo "Masked file created: $output_file"
}

process_D24023_450K_IDAT() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24023/Infinium_Methylation_450K_IDAT_4/"
    local ref_metadata="./metadata/Infinium_Methylation_450K_idat.list"
    local idat_dir="/ARCHIVE/CODA/CODA_D24023/GENETICS/CODA_Infinium_HumanMethylation_450K_IDAT_4/"
    if [[ ! -f "$metadata" ]]; then
        echo "Error: Metadata file '$metadata' not found."
        return 1
    fi
    if [[ ! -f "$ref_metadata" ]]; then
        echo "Error: Reference metadata file '$ref_metadata' not found."
        return 1
    fi
    mkdir -p "$output_dir"
    grep -F -f "$ref_metadata" "$metadata" | cut -f1,4 | while IFS=$'\t' read -r a b; do
        a=$(echo "$a" | xargs)
        b=$(echo "$b" | xargs)
        local grn_file="${idat_dir}${a}_Grn.idat"
        local red_file="${idat_dir}${a}_Red.idat"
        echo "Checking files:"
        echo " - Expected: $grn_file"
        echo " - Expected: $red_file"
        if [[ -f "$grn_file" && -f "$red_file" ]]; then
            echo "Copying: $grn_file -> ${output_dir}${b}_Grn.idat"
            echo "Copying: $red_file -> ${output_dir}${b}_Red.idat"
            cp -v "$grn_file" "${output_dir}${b}_Grn.idat"
            cp -v "$red_file" "${output_dir}${b}_Red.idat"
        else
            echo "Warning: Missing IDAT files for '${a}', skipping..." >&2
        fi
    done
}

process_D24027_450K_CSV() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24027/OGTT1_5th_methylation_50_5gI_work/"
    local ref_metadata="${PWD}/metadata/OGTT1_5th_methylation_50_5gI.list"
    mkdir -p "$output_dir/CODA_D24027/"
    grep -Ff "$ref_metadata" "$metadata" | cut -f4 | tr '\n' ',' | sed 's/,$/\n/' > "$output_dir/header"
    echo "KoBB,$(cat "$output_dir/header")" > "$output_dir/header.csv"
    tail -n +3 /ANALYSIS/SCRIPT/DATA/CODA_D24027/GENETICS/CODA_Infinium_HumanMethylation_450K_csv_50/OGTT1_5th_methylation_50_5gI.csv > "$output_dir/${metadata}.body.tmp"
    local output_file="${PWD}/${metadata%.txt}/CODA_D24027/OGTT1_5th_methylation_50_5gI.${metadata%.txt}.masked.csv"
    cat "$output_dir/header.csv" "$output_dir/${metadata}.body.tmp" > "$output_file"
    echo "Masked file created: $output_file"
}

process_D24027_450K_IDAT() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24027/CODA_Infinium_HumanMethylation_450K_IDAT_400/"
    local ref_metadata="./metadata/CODA_Infinium_HumanMethylation_450K_IDAT_400.list"
    local idat_dir="/ARCHIVE/CODA/CODA_D24027/GENETICS/CODA_Infinium_HumanMethylation_450K_IDAT_400/"
    if [[ ! -f "$metadata" ]]; then
        echo "Error: Metadata file '$metadata' not found."
        return 1
    fi
    if [[ ! -f "$ref_metadata" ]]; then
        echo "Error: Reference metadata file '$ref_metadata' not found."
        return 1
    fi
    mkdir -p "$output_dir"
    grep -F -f "$ref_metadata" "$metadata" | cut -f1,4 | while IFS=$'\t' read -r a b; do
        a=$(echo "$a" | xargs)
        b=$(echo "$b" | xargs)
        local grn_file="${idat_dir}${a}_Grn.idat"
        local red_file="${idat_dir}${a}_Red.idat"
        echo "Checking files:"
        echo " - Expected: $grn_file"
        echo " - Expected: $red_file"
        if [[ -f "$grn_file" && -f "$red_file" ]]; then
            echo "Copying: $grn_file -> ${output_dir}${b}_Grn.idat"
            echo "Copying: $red_file -> ${output_dir}${b}_Red.idat"
            cp -v "$grn_file" "${output_dir}${b}_Grn.idat"
            cp -v "$red_file" "${output_dir}${b}_Red.idat"
        else
            echo "Warning: Missing IDAT files for '${a}', skipping..." >&2
        fi
    done
}

process_D24027_850K_IDAT() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24027/CODA_Infinium_Methylation_EPIC_Array_850K_DI_OB_IDAT_1528/"
    local ref_metadata="./metadata/CODA_Infinium_Methylation_EPIC_Array_850K_DI_OB_IDAT_1528.list"
    local idat_dir="/ARCHIVE/CODA/CODA_D24027/GENETICS/CODA_Infinium_Methylation_EPIC_Array_850K_DI_OB_IDAT_1528/"
    if [[ ! -f "$metadata" ]]; then
        echo "Error: Metadata file '$metadata' not found."
        return 1
    fi
    if [[ ! -f "$ref_metadata" ]]; then
        echo "Error: Reference metadata file '$ref_metadata' not found."
        return 1
    fi
    mkdir -p "$output_dir"
    grep -F -f "$ref_metadata" "$metadata" | cut -f1,4 | while IFS=$'\t' read -r a b; do
        a=$(echo "$a" | xargs)
        b=$(echo "$b" | xargs)
        local grn_file="${idat_dir}${a}_Grn.idat"
        local red_file="${idat_dir}${a}_Red.idat"
        echo "Checking files:"
        echo " - Expected: $grn_file"
        echo " - Expected: $red_file"
        if [[ -f "$grn_file" && -f "$red_file" ]]; then
            echo "Copying: $grn_file -> ${output_dir}${b}_Grn.idat"
            echo "Copying: $red_file -> ${output_dir}${b}_Red.idat"
            cp -v "$grn_file" "${output_dir}${b}_Grn.idat"
            cp -v "$red_file" "${output_dir}${b}_Red.idat"
        else
            echo "Warning: Missing IDAT files for '${a}', skipping..." >&2
        fi
    done
}

process_D24038_HEXA() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24038/CODA_D24038_HEXA_methylation850K_IDAT_822/"
    local ref_metadata="./metadata/CODA_HEXA_methylation850K_IDAT_822_idat.list"
    local idat_dir="/ARCHIVE/CODA/CODA_D24038/GENETICS/CODA_HEXA_methylation850K_IDAT_822/"
    if [[ ! -f "$metadata" ]]; then
        echo "Error: Metadata file '$metadata' not found."
        return 1
    fi
    if [[ ! -f "$ref_metadata" ]]; then
        echo "Error: Reference metadata file '$ref_metadata' not found."
        return 1
    fi
    mkdir -p "$output_dir"
    grep -F -f "$ref_metadata" "$metadata" | cut -f1,4 | while IFS=$'\t' read -r a b; do
        a=$(echo "$a" | xargs)
        b=$(echo "$b" | xargs)
        local grn_file="${idat_dir}${a}_Grn.idat"
        local red_file="${idat_dir}${a}_Red.idat"
        echo "Checking files:"
        echo " - Expected: $grn_file"
        echo " - Expected: $red_file"
        if [[ -f "$grn_file" && -f "$red_file" ]]; then
            echo "Copying: $grn_file -> ${output_dir}${b}_Grn.idat"
            echo "Copying: $red_file -> ${output_dir}${b}_Red.idat"
            cp -v "$grn_file" "${output_dir}${b}_Grn.idat"
            cp -v "$red_file" "${output_dir}${b}_Red.idat"
        else
            echo "Warning: Missing IDAT files for '${a}', skipping..." >&2
        fi
    done
}

process_D24038_Infinium() {
    local metadata="$1"
    local output_dir="${PWD}/${metadata%.txt}/CODA_D24038/CODA_HEXA_methylation850K_IDAT_822/"
    local ref_metadata="./metadata/CODA_HEXA_methylation850K_IDAT_822.list"
    local idat_dir="/ARCHIVE/CODA/CODA_D24038/GENETICS/CODA_HEXA_methylation850K_IDAT_822/"
    if [[ ! -f "$metadata" ]]; then
        echo "Error: Metadata file '$metadata' not found."
        return 1
    fi
    if [[ ! -f "$ref_metadata" ]]; then
        echo "Error: Reference metadata file '$ref_metadata' not found."
        return 1
    fi
    mkdir -p "$output_dir"
    grep -F -f "$ref_metadata" "$metadata" | cut -f1,4 | while IFS=$'\t' read -r a b; do
        a=$(echo "$a" | xargs)
        b=$(echo "$b" | xargs)
        local grn_file="${idat_dir}${a}_Grn.idat"
        local red_file="${idat_dir}${a}_Red.idat"
        echo "Checking files:"
        echo " - Expected: $grn_file"
        echo " - Expected: $red_file"
        if [[ -f "$grn_file" && -f "$red_file" ]]; then
            echo "Copying: $grn_file -> ${output_dir}${b}_Grn.idat"
            echo "Copying: $red_file -> ${output_dir}${b}_Red.idat"
            cp -v "$grn_file" "${output_dir}${b}_Grn.idat"
            cp -v "$red_file" "${output_dir}${b}_Red.idat"
        else
            echo "Warning: Missing IDAT files for '${a}', skipping..." >&2
        fi
    done
}

##################################
# Main: Process each selected dataset type in parallel
##################################
pids=()
for dtype in "${selected_types[@]}"; do
    metadata_file="${type_input_files[$dtype]}"
    echo "Processing dataset type: $dtype with metadata file: $metadata_file"
    case "$dtype" in
        D24023_CGH)
            process_D24023_CGH "$metadata_file" &
            pids+=($!)
            ;;
        D24023_Exome)
            process_D24023_Exome "$metadata_file" &
            pids+=($!)
            ;;
        D24023_450K_CSV)
            process_D24023_450K_CSV "$metadata_file" &
            pids+=($!)
            ;;
        D24023_450K_IDAT)
            process_D24023_450K_IDAT "$metadata_file" &
            pids+=($!)
            ;;
        D24027_450K_CSV)
            process_D24027_450K_CSV "$metadata_file" &
            pids+=($!)
            ;;
        D24027_450K_IDAT)
            process_D24027_450K_IDAT "$metadata_file" &
            pids+=($!)
            ;;
        D24027_850K_IDAT)
            process_D24027_850K_IDAT "$metadata_file" &
            pids+=($!)
            ;;
        D24038_HEXA)
            process_D24038_HEXA "$metadata_file" &
            pids+=($!)
            ;;
        D24038_Infinium)
            process_D24038_Infinium "$metadata_file" &
            pids+=($!)
            ;;
        *)
            echo "Error: Unknown dataset type '$dtype'"
            ;;
    esac
done

# Wait for all background processes to complete.
for pid in "${pids[@]}"; do
    wait "$pid"
done
