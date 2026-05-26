import os
import numpy as np
import uproot
import polars as pl
import pyarrow as pa

def parquet_to_root_uproot(parquet_path: str, root_path: str, tree_name: str = "kinematics") -> None:
    df = pl.read_parquet(parquet_path)

    arrow_table = df.to_arrow()
    for i, field in enumerate(arrow_table.schema):
        if field.type == pa.large_utf8():
            arrow_table = arrow_table.set_column(i, field.name, arrow_table.column(field.name).cast(pa.string()))

    data = {
        name: arrow_table.column(name).to_pylist() if arrow_table.schema.field(name).type == pa.string()
        else arrow_table.column(name).to_numpy()
        for name in arrow_table.schema.names
    }

    branch_types = {
        name: arr.dtype if isinstance(arr, np.ndarray) else str
        for name, arr in data.items()
    }

    with uproot.recreate(root_path) as root_file:
        root_file.mktree(tree_name, branch_types)
        root_file[tree_name].extend(data)

directory = "/home1/georgina.xifra/attpc_engine/"
output_directory = "/home1/georgina.xifra/attpc_engine/"

if not os.path.isdir(output_directory):
    os.makedirs(output_directory)
    print(f"{output_directory} created successfully!")

parquet_files = ["output_16Cpd.parquet", "output_16Cpt.parquet"]

for file in parquet_files:
    print(f"Converting {file}...")
    filename = os.path.splitext(file)[0]
    try:
        parquet_to_root_uproot(os.path.join(directory, file), os.path.join(output_directory, filename + ".root"))
        print(f"  Done -> {filename}.root")
    except Exception as e:
        print(f"  Failed: {e}")
