import pandas as pd
import os

directory = "/home/georgina/Documents/ATTPC/scalers"
parquet_file = os.path.join(directory, "run_0033_scaler.parquet")  # Change to any file you want to check

# Load the Parquet file
df = pd.read_parquet(parquet_file)

# Show column information
print(df.info())  # Shows column types and number of nulls
print(df.head())  # Shows the first few rows

