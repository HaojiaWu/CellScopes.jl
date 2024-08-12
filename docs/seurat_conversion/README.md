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
The xenium rds file is prepared following this guide: https://satijalab.org/seurat/articles/spatial_vignette_2.html#mouse-brain-10x-genomics-xenium-in-situ. 
After the rds file is ready, you can convert it to a CellScopes XeniumObject.
```julia
seurat_xenium=cs.from_seurat("xenium_seurat.rds"; 
    data_type="spatial",
    tech="xenium", 
    assay="Xenium")
```
This is how it looks after the object is created.
```
XeniumObject in CellScopes.jl
Genes x Cells = 248 x 36553
Available data:
- rawCount
- normCount
- metaData
- spmetaData
- dimReduction
- clustData
- polynormCount
- coordData
- polygonData
All fields:
- rawCount
- normCount
- scaleCount
- metaData
- spmetaData
- varGene
- dimReduction
- clustData
- polyCount
- polynormCount
- coordData
- imputeData
- imageData
- polygonData
```
#### 2.1 plot cell annoation on umap
```julia
cs.dim_plot(seurat_xenium; legend_ncol=2)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seurat_xenium_umap.png" width="400"> <br>

#### 2.2 plot cell annoation on tissue
```julia
cs.sp_dim_plot(seurat_xenium, "cluster"; 
    width=600, height = 500, 
    marker_size=3,legend_ncol=2)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seurat_xenium_tissue.png" width="400"> <br>

#### 2.3 plot gene expression on tissue
```julia
cs.sp_feature_plot(seurat_xenium, ["Rorb", "Bcl11b"]; 
    color_keys=["gray94", "lemonchiffon", "red"])
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seurat_xenium_genes.png" width="800"> <br>

### 3. From a Visium Seurat object
Here is a short tuturial on converting a Visium Seurat object to a CellScopes Visium Object. The Visium seurat object is prepared by following this tutorial: https://satijalab.org/seurat/articles/spatial_vignette.html. Then you can read the rds file and convert it into a CellScopes object:
```julia
seurat_visium = cs.from_seurat("brain_visium.rds"; 
    data_type="spatial",tech="visium", assay="SCT")
```
#### 3.1 plot cell annoation on tissue
```julia
cs.sp_dim_plot(seurat_visium, "cluster"; 
    marker_size = 8,  width=600, height = 500,  
    do_label=false, alpha=0.5, img_res="high")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seurat_visium_tissu.png" width="400"> <br>

#### 3.2 plot gene expression on tissue
```julia
cs.sp_feature_plot(seurat_visium, ["Gng4","Nrgn", "Ttr"]; 
    marker_size = 8, color_keys=["gray90", "lemonchiffon" ,"red"])
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seurat_visium_genes.png" width="1200"> <br>

