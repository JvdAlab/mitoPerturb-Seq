##### pySCENIC 
##### add metadata columns to loom file produced by aucell step of pyscenic 

# load pyscenic environment 
module load miniconda/2024-02-20
conda deactivate
conda activate /../../.conda/envs/pyscenic/

# load python
python

# load packages
import h5py
import pandas as pd
import numpy as np

auc_loom_file = "/../../project_rna_assay_pyscenic_output_filtered.loom"

metadata_file = "/../../project_metadata.csv"
metadata = pd.read_csv(metadata_file, index_col=0)
metadata.index = metadata.index.astype(str)

with h5py.File(auc_loom_file, "r+") as f:
   loom_cells = [c.decode("utf-8") for c in f["col_attrs/CellID"][:]]
   metadata = metadata.loc[loom_cells]
   for col in metadata.columns:
           values = metadata[col].astype(str).values.astype("S")
           if f"col_attrs/{col}" in f:
                   print(f"Overwriting existing column: {col}")
                   del f[f"col_attrs/{col}"]
           print(f"Adding metadata column: {col}")
           f.create_dataset(f"col_attrs/{col}", data=values)

print("✅ Metadata successfully added to loom file!")

with h5py.File(auc_loom_file, "r") as f:
   print(list(f["col_attrs"].keys()))

