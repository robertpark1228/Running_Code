#!/bin/bash
# bed_ped_multiple_updates.sh
#
# Usage: $0 -m id_mapping.txt [-b input.bed] [-p input.ped] [-e exception_list.txt] [-s coda_id] [-f bed_file_list]
#
# This script updates a PLINK binary fileset (.bed, .bim, .fam) or a .ped file with new Family and Individual IDs.
#
# Options:
#   -m  (Required) ID mapping table (format: old_family_id old_individual_id new_family_id new_individual_id)
#   -b  (Optional) Input PLINK .bed file (requires corresponding .bim and .fam)
#   -f  (Optional) File containing a list of .bed files to process
#   -p  (Optional) Input PLINK .ped file (if renaming individuals in PED format)
#   -e  (Optional) Exception list (IDs to exclude from renaming)
#   -s  (Optional) Specify a CODA_D[XXXXX] dataset manually
#   -h  Show this help message
#
# If neither a bed file (-b) nor a bed file list (-f) is provided,
# the script uses the predefined list of bed file paths below and presents a numbered menu for selection.
#
# The mapping file (-m) must be a tab-separated TXT file with the columns:
# ORIGINAL_FAM_ID, ORIGINAL_SAMPLE_ID, CONVERTED_FAM_ID, CONVERTED_SAMPLE_ID.
#
# All selected bed file updates run in parallel.
#
# Example:
#   ./bed_ped_multiple_updates.sh -m id_mapping.txt
#

usage() {
    echo "Usage: $0 -m id_mapping.txt [-b input.bed] [-p input.ped] [-e exception_list.txt] [-s coda_id] [-f bed_file_list]"
    echo ""
    echo "Options:"
    echo "  -m  (Required) ID mapping table (format: old_family_id old_individual_id new_family_id new_individual_id)"
    echo "  -b  (Optional) Input PLINK .bed file (requires corresponding .bim and .fam)"
    echo "  -f  (Optional) File containing a list of .bed files to process"
    echo "  -p  (Optional) Input PLINK .ped file (if renaming individuals in PED format)"
    echo "  -e  (Optional) Exception list (IDs to exclude from renaming)"
    echo "  -s  (Optional) Specify a CODA_D[XXXXX] dataset manually"
    echo "  -h  Show this help message"
    exit 1
}

# Predefined list of bed files from the provided list.
PREDEFINED_BED_FILES=(
"/ARCHIVE/CODA/CODA_D24027/GENETICS/CODA_Exomechip_BED_5098/CODA_Exomechip_BED_5098_KA_giban5.bed"
"/ARCHIVE/CODA/CODA_D24033/GENETICS/CODA_Affy6_BED_1816/CODA_CAVAS_Affy6_1816_maskingID.bed"
"/ARCHIVE/CODA/CODA_D24033/GENETICS/CODA_BED_Exomechip_KN_3068/CODA_Exomechip_14025_KN.bed"
"/ARCHIVE/CODA/CODA_D24033/GENETICS/CODA_Illumina_Omni1_BED_3666/CODA_CAVAS_IlluminaOmni1_3666_maskingID.bed"
"/ARCHIVE/CODA/CODA_D24033/GENETICS/CODA_Korea_Biobank_Array_GNT_BED_8105/CODA_Korea_Biobank_Array_GNT_BED_8105.bed"
"/ARCHIVE/CODA/CODA_D24038/GENETICS/CODA_HEXA_Affy6_BED_3693/CODA_HEXA_Affy6_BED_3693_maskingID.bed"
"/ARCHIVE/CODA/CODA_D24038/GENETICS/CODA_Exomechip_BED_3433/CODA_Exomechip_3433_KC.bed"
"/ARCHIVE/CODA/CODA_D24038/GENETICS/CODA_Korea_Biobank_Array_GNT_BED_58693/CODA_Korea_Biobank_Array_GNT_BED_58693.bed"
"/ARCHIVE/CODA/CODA_D24023/GENETICS/CODA_Korea_Biobank_Array_GNT_BED_5493/CODA_Korea_Biobank_Array_GNT_BED_5493.bed"
"/ARCHIVE/CODA/CODA_D24023/GENETICS/CODA_Affy5.0_BED_8840/CODA_KARE_affy5_8840_maskingID.bed"
"/ARCHIVE/CODA/CODA_D24023/GENETICS/CODA_Affy5.0_HapMap_BED_8840/CODA_KARE_HapMap_8840_maskingID.bed"
"/ARCHIVE/CODA/CODA_D24023/GENETICS/CODA_Exomechip_BED_2426/CODA_Exomechip_2426.bed"
"/ARCHIVE/CODA/CODA_D24023/GENETICS/CODA_Affy5.0_1KG_BED_8840/CODA_KARE_1KG_8840_maskingID.bed"
)

# Parse arguments.
parse_arguments() {
    while getopts "b:m:e:p:s:f:h" opt; do
        case ${opt} in
            b ) BED_FILE="$OPTARG" ;;
            m ) MAPPING_FILE="$OPTARG" ;;
            e ) EXCEPTION_FILE="$OPTARG" ;;
            p ) PED_FILE="$OPTARG" ;;
            s ) CODA_DIR_OPTION="$OPTARG" ;;
            f ) BED_FILE_LIST="$OPTARG" ;;
            h ) usage ;;
            * ) usage ;;
        esac
    done

    # If mapping file is not provided, prompt for it.
    if [[ -z "$MAPPING_FILE" ]]; then
        read -p "Enter the ID mapping file path: " MAPPING_FILE
    fi
    if [[ ! -f "$MAPPING_FILE" ]]; then
        echo "Error: Mapping file '$MAPPING_FILE' does not exist. Exiting."
        exit 1
    fi
    PROJECT_NAME="${MAPPING_FILE%.txt}"

    # Determine bed files.
    if [[ -n "$BED_FILE_LIST" ]]; then
        if [[ ! -f "$BED_FILE_LIST" ]]; then
            echo "Error: BED file list '$BED_FILE_LIST' does not exist. Exiting."
            exit 1
        fi
        mapfile -t BED_FILES < "$BED_FILE_LIST"
    elif [[ -n "$BED_FILE" ]]; then
        BED_FILES=("$BED_FILE")
    else
        # Use the predefined list.
        BED_FILES=("${PREDEFINED_BED_FILES[@]}")
        echo "Using predefined bed file list:"
        for i in "${!BED_FILES[@]}"; do
            printf "  %d) %s\n" "$((i+1))" "${BED_FILES[$i]}"
        done
        read -p "Enter your choices (e.g., 1,3,5): " choice_input
        IFS=',' read -ra choices <<< "$(echo "$choice_input" | tr -d ' ')"
        selected=()
        for choice in "${choices[@]}"; do
            index=$((choice-1))
            if [[ $index -ge 0 && $index -lt ${#BED_FILES[@]} ]]; then
                selected+=("${BED_FILES[$index]}")
            else
                echo "Invalid selection: $choice"
            fi
        done
        if [[ ${#selected[@]} -eq 0 ]]; then
            echo "Error: No valid .bed file selected. Exiting."
            exit 1
        fi
        BED_FILES=("${selected[@]}")
    fi

    echo "Mapping file: $MAPPING_FILE"
    echo "Selected .bed files: ${BED_FILES[@]}"
}

# Process a single BED file.
process_bed_file() {
    local BED_FILE="$1"
    echo "Processing: $BED_FILE"
    if [[ ! -f "$BED_FILE" ]]; then
        echo "Error: File not found - $BED_FILE"
        return
    fi

    local BASE_NAME
    BASE_NAME=$(basename "$BED_FILE" .bed)

    # Determine CODA_D[XXXXX] dataset from the file path or use the specified option.
    local CODA_DIR
    CODA_DIR=$(echo "$BED_FILE" | grep -oE 'CODA_D[0-9]+' | head -n 1)
    if [[ -z "$CODA_DIR" && -n "$CODA_DIR_OPTION" ]]; then
        CODA_DIR="$CODA_DIR_OPTION"
    fi
    if [[ -z "$CODA_DIR" ]]; then
        echo "Error: Could not determine CODA_D[XXXXX] dataset for $BED_FILE"
        return
    fi

    local ORIGINAL_DIR="./${PROJECT_NAME}/${CODA_DIR}/original"
    local UPDATED_DIR="./${PROJECT_NAME}/${CODA_DIR}/updated"
    mkdir -p "$ORIGINAL_DIR" "$UPDATED_DIR"

    local FAM_FILE="${BED_FILE%.bed}.fam"
    local BIM_FILE="${BED_FILE%.bed}.bim"

    if [[ ! -f "$FAM_FILE" || ! -f "$BIM_FILE" ]]; then
        echo "Error: Missing PLINK files for $BED_FILE (.bim, .fam). Skipping."
        return
    fi

    # Copy original files if not already present.
    if [[ ! -f "${ORIGINAL_DIR}/${BASE_NAME}.bed" ]]; then
        cp "$BED_FILE" "${ORIGINAL_DIR}/${BASE_NAME}.bed"
        cp "$FAM_FILE" "${ORIGINAL_DIR}/${BASE_NAME}.fam"
        cp "$BIM_FILE" "${ORIGINAL_DIR}/${BASE_NAME}.bim"
    fi

    # Generate update file for PLINK using the mapping file.
    awk 'NR==FNR {map[$1" "$2] = $3" "$4; next} 
         ($1" "$2) in map {print $1, $2, map[$1" "$2]}' "$MAPPING_FILE" "${ORIGINAL_DIR}/${BASE_NAME}.fam" > "${UPDATED_DIR}/${BASE_NAME}_update_ids.txt"

    # Run PLINK to update the files.
    plink --bfile "${ORIGINAL_DIR}/${BASE_NAME}" --update-ids "${UPDATED_DIR}/${BASE_NAME}_update_ids.txt" --make-bed --out "${UPDATED_DIR}/${BASE_NAME}"
    
    if [[ $? -eq 0 ]]; then
        echo "PLINK files updated successfully for $BED_FILE."
        echo "Original files stored in: ${ORIGINAL_DIR}"
        echo "Updated files stored in: ${UPDATED_DIR}"
    else
        echo "Error: PLINK update failed for $BED_FILE."
    fi
}

# Main: parse arguments and process each selected bed file in parallel.
main() {
    parse_arguments "$@"
    pids=()
    for file in "${BED_FILES[@]}"; do
        process_bed_file "$file" &
        pids+=($!)
    done
    for pid in "${pids[@]}"; do
        wait "$pid"
    done
}

main "$@"
