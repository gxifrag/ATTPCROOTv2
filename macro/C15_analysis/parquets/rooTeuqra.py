from parquet_to_root import parquet_to_root
import os

directory = "/home/georgina/Documents/ATTPC/scalers"
output_directory = "/home/georgina/Documents/ATTPC/rootfilesScalers"

if not os.path.isdir(output_directory):
    os.makedirs(output_directory)
    print(f"{output_directory} created successfully!")

for file in os.listdir(directory):
    print(file)
    filename, extension = os.path.splitext(file)
    if extension == ".parquet":
        parquet_to_root(os.path.join(directory, file), os.path.join(output_directory, filename + ".root"))
    else:
        print(f"{file} is not a .parquet file.")
