# Imports

import pandas as pd

# Remove irrelevant data types

df = pd.read_csv("/content/patient28_adata_ccs.csv")
celltypes = []

for row in df.itertuples(index=False):
  if row.leiden_r06 == "CAF" or row.leiden_r06 == "CD8 T cell" or row.leiden_r06 == "Tumour epithelial" or row.leiden_r06 == "Tumour epithelial (proliferative)":
    celltypes.append({
        "leiden_r06": row.leiden_r06,
        "x_aligned": row.x_aligned,
        "y_aligned": row.y_aligned,
        "CD274": row.CD274
    })

pd.DataFrame(celltypes).to_csv("patient28_celltypes.csv", index=False)

# Filter out some cells, normalize, reorder

df = pd.read_csv("patient28_celltypes.csv")
filtered = []
proportion = 4

counter = 0
for row in df.itertuples(index=False):
    if row.leiden_r06 == "CAF":
      if counter % proportion == 0:
        filtered.append({
          "leiden_r06": row.leiden_r06,
          "x_aligned": (row.x_aligned),
          "y_aligned": (row.y_aligned),
          "CD274": row.CD274
        })
      counter+=1

counter = 0
for row in df.itertuples(index=False):
    if row.leiden_r06 == "Tumour epithelial" or row.leiden_r06 == "Tumour epithelial (proliferative)":
      if counter % proportion == 0:
        filtered.append({
          "leiden_r06": row.leiden_r06,
          "x_aligned": (row.x_aligned),
          "y_aligned": (row.y_aligned),
          "CD274": row.CD274
        })
      counter+=1

counter = 0
for row in df.itertuples(index=False):
    if row.leiden_r06 == "CD8 T cell":
      if counter % proportion == 0:
        filtered.append({
          "leiden_r06": row.leiden_r06,
          "x_aligned": (row.x_aligned),
          "y_aligned": (row.y_aligned),
          "CD274": row.CD274
        })
      counter+=1

pd.DataFrame(filtered).to_csv("patient28_filtered.csv", index=False)

df = pd.read_csv("patient28_filtered.csv")
normalized_reordered = []

for row in df.itertuples(index=False):
    if row.leiden_r06 == "CAF":
      normalized_reordered.append({
        "leiden_r06": row.leiden_r06,
        "x_aligned": (row.x_aligned)/4,
        "y_aligned": (row.y_aligned)/4,
        "CD274": row.CD274
      })

for row in df.itertuples(index=False):
    if row.leiden_r06 == "Tumour epithelial" or row.leiden_r06 == "Tumour epithelial (proliferative)":
      normalized_reordered.append({
        "leiden_r06": row.leiden_r06,
        "x_aligned": (row.x_aligned)/4,
        "y_aligned": (row.y_aligned)/4,
        "CD274": row.CD274
      })

for row in df.itertuples(index=False):
    if row.leiden_r06 == "CD8 T cell":
      normalized_reordered.append({
        "leiden_r06": row.leiden_r06,
        "x_aligned": (row.x_aligned)/4,
        "y_aligned": (row.y_aligned)/4,
        "CD274": row.CD274
      })

pd.DataFrame(normalized_reordered).to_csv("patient28_filtered_edited.csv", index=False)
