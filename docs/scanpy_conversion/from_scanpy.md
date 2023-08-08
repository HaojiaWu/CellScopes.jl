## Converstion of scanpy anndata to CellScopes object
If existing anndata generated from the scanpy standard pipeline was saved in the hard drive (hdf5 format), we provide a function to directly read this file and convert it into a corresponding CellScope object. The current version support conversion of scRNA-seq, Xenium and Visium data from scanpy. Note that scanpy needs to be pre-installed into your current environment following the scanpy tutorial:
https://scanpy.readthedocs.io/en/stable/installation.html
```julia
import CellScopes as cs
```
#### 1. From scRNA anndata
We created a function called ```from_scanpy``` to read the h5ad file and convert it to a ```scRNAObject``` in CellScopes. The h5ad file ```pbmc3k.h5ad``` below was created following a tutorial for processing a PBMC dataset.
https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
```julia
scanpy_sc = cs.from_scanpy("pbmc3k.h5ad"; data_type="scRNA", anno="leiden")
```
This should create a scRNAObject that is suitable for downstream analysis or visiualization using CellScopes.
```
scRNAObject in CellScopes.jl
Genes x Cells = 1838 x 2638
Available data:
- Raw count
- Normalized count
- Scaled count
- Clustering data
- UMAP data
All fields:
- rawCount
- normCount
- scaleCount
- metaData
- varGene
- dimReduction
- clustData
- undefinedData
```
##### 1.1 Visualize cells


