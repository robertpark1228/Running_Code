import pandas as pd
import pyreadstat
import argparse

def convert_sas_to_tsv(input_file, output_file):
    # Read the SAS dataset
    df, meta = pyreadstat.read_sas7bdat(input_file)

    # Save as TSV
    df.to_csv(output_file, sep="\t", index=False)

    print(f"Conversion completed: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert SAS7BDAT file to TSV format.")
    parser.add_argument("input_file", help="Path to the input SAS7BDAT file")
    parser.add_argument("output_file", help="Path to the output TSV file")

    args = parser.parse_args()

    convert_sas_to_tsv(args.input_file, args.output_file)
