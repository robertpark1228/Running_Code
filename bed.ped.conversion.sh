#!/bin/bash

# Help message
usage() {
    echo "Usage: $0 -m id_mapping.txt [-b input.bed] [-p input.ped] [-e exception_list.txt] [-s coda_id] [-f bed_file_list]"
    echo ""
    echo "This script updates a PLINK binary fileset (.bed, .bim, .fam) or a .ped file with new Family and Individual IDs."
    echo "Options:"
    echo "  -b  (Optional) Input PLINK .bed file (requires corresponding .bim and .fam)"
    echo "  -m  ID mapping table (format: old_family_id old_individual_id new_family_id new_individual_id) (Required)"
    echo "  -e  (Optional) Exception list (IDs to exclude from renaming)"
    echo "  -p  (Optional) Input PLINK .ped file (if renaming individuals in a PED format)"
    echo "  -s  (Optional) Specify a CODA_D[XXXXX] dataset manually"
    echo "  -f  (Optional) File containing a list of .bed files to process"
    echo "  -h  Show this help message"
    exit 1
}

# Process arguments
parse_arguments() {
    while getopts "b:m:e:p:s:f:h" opt; do
        case ${opt} in
            b ) BED_FILE="$OPTARG" ;;
            m ) MAPPING_FILE="$OPTARG" ;;
            e ) EXCEPTION_FILE="$OPTARG" ;;
            p ) PED_FILE="$OPTARG" ;;
            s ) CODA_DIR="$OPTARG" ;;
            f ) BED_FILE_LIST="$OPTARG" ;;
            h ) usage ;;
            * ) usage ;;
        esac
    done

    if [[ -z "$MAPPING_FILE" ]]; then
        echo "Error: Missing required argument -m (ID mapping file)."
        usage
    fi

    # Extract project folder name from the mapping file (removing .txt)
    PROJECT_NAME="${MAPPING_FILE%.txt}"

    if [[ -n "$BED_FILE_LIST" ]]; then
        mapfile -t BED_FILES < "$BED_FILE_LIST"
    elif [[ -n "$BED_FILE" ]]; then
        BED_FILES=("$BED_FILE")
    else
        echo "Error: No .bed file provided. Use -b or -f."
        usage
    fi

    echo "Processing BED files: ${BED_FILES[@]}"
}

# Process each BED file
process_bed_file() {
    local BED_FILE="$1"
    echo "Processing: $BED_FILE"
    if [[ ! -f "$BED_FILE" ]]; then
        echo "Error: File not found - $BED_FILE"
        return
    fi
    
    local BASE_NAME=$(basename "$BED_FILE" .bed)

    # Extract CODA_D[XXXXX] from file path or use the specified one
    local CODA_DIR=$(echo "$BED_FILE" | grep -oE 'CODA_D[0-9]+' | head -n 1)
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

    if [[ ! -f "$BIM_FILE" || ! -f "$FAM_FILE" ]]; then
        echo "Error: Missing PLINK files for $BED_FILE (.bim, .fam). Skipping."
        return
    fi

    # Copy original files to the original folder if not already copied
    if [[ ! -f "${ORIGINAL_DIR}/${BASE_NAME}.bed" ]]; then
        cp "$BED_FILE" "${ORIGINAL_DIR}/${BASE_NAME}.bed"
        cp "$FAM_FILE" "${ORIGINAL_DIR}/${BASE_NAME}.fam"
        cp "$BIM_FILE" "${ORIGINAL_DIR}/${BASE_NAME}.bim"
    fi

    # Generate update file for PLINK
    awk 'NR==FNR {map[$1" "$2] = $3" "$4; next} 
        ($1" "$2) in map {print $1, $2, map[$1" "$2]} 
        ' "$MAPPING_FILE" "${ORIGINAL_DIR}/${BASE_NAME}.fam" > "${UPDATED_DIR}/${BASE_NAME}_update_ids.txt"

    # Run PLINK to generate updated files
    plink --bfile "${ORIGINAL_DIR}/${BASE_NAME}" --update-ids "${UPDATED_DIR}/${BASE_NAME}_update_ids.txt" --make-bed --out "${UPDATED_DIR}/${BASE_NAME}"

    if [[ $? -eq 0 ]]; then
        echo "PLINK files updated successfully for $BED_FILE."
        echo "Original files stored in: ${ORIGINAL_DIR}"
        echo "Updated files stored in: ${UPDATED_DIR}"
    else
        echo "Error: PLINK update failed for $BED_FILE."
    fi
}

# Main script execution
main() {
    parse_arguments "$@"
    for BED_FILE in "${BED_FILES[@]}"; do
        process_bed_file "$BED_FILE" &
    done
}

main "$@"
