## Conversion of Seurat object to CellScopes object
If you are a Seurat user and interested in exploring the functionalities within ```CellScopes```, the ```from_seurat``` function is at your disposal. It enables the direct reading of Seurat RDS files, converting them into a CellScope object. As of now, the ```from_seurat``` function can read Seurat objects that have been prepared from scRNA-seq, Xenium, and Visium analyses. Note that the Seurat package must be pre-installed in your current R environment for this functionality to work properly. For installation of Seurat, please visit Satija's lab for more detail. https://satijalab.org/seurat/
```julia
import CellScopes as cs
```
### 1. From a scRNA-seq Seurat object
Here is how to read the rds file from a scRNA-seq Seurat object. The Seurat object ```pbmc_tutorial.rds``` was prepared using the tutorial from a PBMC analysis.
https://satijalab.org/seurat/articles/pbmc3k_tutorial
```julia
seurat_sc = cs.from_seurat("pbmc_tutorial.rds";
        data_type="scRNA", 
        assay="RNA")
```
This should create a scRNAObject.
```
scRNAObject in CellScopes.jl
Genes x Cells = 13714 x 2700
Available data:
- Raw count
- Normalized count
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
#### 1.1 Visualize cell annotaion in umap
```julia
cs.dim_plot(pbmc; marker_size=5)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seurat_sc_umap.png" width="400"> <br>

#### 1.2 Visualize gene expression
```julia
cs.feature_plot(seurat_sc, ["CST3", "NKG7", "PPBP"]; marker_size=6, dim_type="umap")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seurat_sc_genes.png" width="1200"> <br>

### 2. From a Xenium Seurat object




