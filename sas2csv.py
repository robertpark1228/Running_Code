import pandas as pd
import pyreadstat
import argparse
import os

def convert_csv_to_sas(csv_file, sas_file):
    try:
        # Check if the input file exists
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"Error: The file '{csv_file}' does not exist.")

        # Load CSV into DataFrame
        df = pd.read_csv(csv_file)

        # Convert to SAS7BDAT
        pyreadstat.write_sas(df, sas_file, format="sas7bdat")
        
        print(f"✅ Successfully converted '{csv_file}' to '{sas_file}'")

    except FileNotFoundError as e:
        print(e)
    except pd.errors.EmptyDataError:
        print("Error: The CSV file is empty.")
    except pd.errors.ParserError:
        print("Error: The CSV file contains malformed data.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a CSV file to SAS7BDAT format.",
        epilog="Example usage: python csv_to_sas.py input.csv output.sas7bdat"
    )
    
    parser.add_argument("csv_file", help="Path to the input CSV file")
    parser.add_argument("sas_file", help="Path to the output SAS7BDAT file")

    args = parser.parse_args()

    convert_csv_to_sas(args.csv_file, args.sas_file)
