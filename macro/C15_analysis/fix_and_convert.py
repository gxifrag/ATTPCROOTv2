import os
import pandas as pd
from parquet_to_root import parquet_to_root

# Paths (use the directories you've specified)
input_directory = "/home/georgina/Documents/ATTPC/scalers"
output_directory = "/home/georgina/Documents/ATTPC/rootfilesScalers"

# Create output directory if it doesn't exist
if not os.path.isdir(output_directory):
    os.makedirs(output_directory)
    print(f"{output_directory} created successfully!")

# Loop through all .parquet files
for file in os.listdir(input_directory):
    filename, extension = os.path.splitext(file)
    
    if extension == ".parquet":
        parquet_path = os.path.join(input_directory, file)
        root_path = os.path.join(output_directory, filename + ".root")

        print(f"Processing: {file}")

        # Load Parquet file
        df = pd.read_parquet(parquet_path)

        # Ensure all columns are correctly typed
        df = df.astype({col: "int64" for col in df.select_dtypes(include=["int", "float"]).columns})

        # Save fixed Parquet file
        fixed_parquet_path = os.path.join(input_directory, f"fixed_{file}")
        df.to_parquet(fixed_parquet_path, index=False)
        
        # Convert to ROOT
        parquet_to_root(fixed_parquet_path, root_path)
        print(f"Converted {file} -> {filename}.root ✅")
        
    else:
        print(f"{file} is not a .parquet file.")

