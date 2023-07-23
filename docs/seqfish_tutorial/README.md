# Analysis and visualization of seqFISH data with CellScopes.jl
[seqFISH (Sequential Fluorescence In Situ Hybridization)](https://spatial.caltech.edu/seqfish/) is an advanced microscopy technique used to visualize the spatial arrangement of individual mRNA molecules within cells. This tutorial walks you through how to use ```CellScopes``` to analyze [a embryo seqFISH dataset](https://www.nature.com/articles/s41587-021-01006-2).

## 1. Read input files and create a seqFISH object
The seqFISH data were provided as R object format (https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/) so we need a proprecessing step in R to create the transcript, cell, and gene-by-cell count files. Here is the example R script to process the original files and save as csv formats.

```R
##The following codes were run in R
library(Matrix)
## Read the transcripts, cells and count files
mRNA <- readRDS("mRNA.Rds")
counts <- readRDS("counts.Rds")
cells <- readRDS("metadata.Rds")
### Save as CSV files
write.csv(mRNA, "transcripts.csv")
cells<- cells[,1:14]
write.csv(cells, "cells.csv")
counts <- data.frame(as.matrix(counts))
write.csv(counts, "counts.csv")
```

## 2. Read input files and create a seqFISH object
At this point, we can transition back to Julia to continue processing the input files in the appropriate format that CellScopes can further process. The seqFISH data contain three embryo dataset. We only use one of them (embryo1) in this tutorial.

### a. Read all files and reformat them. Note that only "embryo1" was used for downstream analysis.
```julia
import CellScopes as cs
using DataFrames, CSV, SparseArrays
transcripts = CSV.read("transcripts.csv", DataFrame)
cells = CSV.read("cells.csv", DataFrame)
em1_cell = filter(:embryo => ==("embryo1"), cells)
rename!(em1_cell, :x_global=>:x, :y_global => :y);
em1_transcript = filter(:uniqueID => ∈(Set(em1_cell.uniqueID)), transcripts)
rename!(em1_transcript, :geneID => :gene,:x => :x_local, :y => :y_local, :x_global=>:x, :y_global => :y)
em1_transcript.x = em1_transcript.x .- minimum(em1_transcript.x)
em1_transcript.y = em1_transcript.y .- minimum(em1_transcript.y)
em1_cell.x = em1_cell.x .- minimum(em1_cell.x)
em1_cell.y = em1_cell.y .- minimum(em1_cell.y)
counts = CSV.read("counts.csv", DataFrame)
em1_counts = counts[:, ["Column1"; unique(em1_cell.uniqueID)]]
em1_cell.cell = string.(collect(1:length(em1_cell.uniqueID)))
em1_transcript = cs.map_values(em1_transcript, :uniqueID, :cell, em1_cell.uniqueID, em1_cell.cell)
em1_cell=em1_cell[!, [["cell","x","y"]; setdiff(names(em1_cell), ["cell","x","y"])]]
em1_transcript=em1_transcript[!, [["cell","x","y"]; setdiff(names(em1_transcript), ["cell","x","y"])]];
```
### b. Create a seqFISH object

```julia
gene_name = em1_counts.Column1
cell_name = em1_cell.cell
count_df = em1_counts[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = cs.RawCountObject(count_df, cell_name, gene_name)
embryo1 = cs.CartanaObject(em1_transcript, em1_cell, raw_count;
    prefix = "embryo1", min_gene = 0, min_cell = 3)
```

This should create a ```seqFISHObject``` object.
```
Adding prefix embryo1 to all cells...
seqFISHObject in CellScopes.jl
Genes x Cells = 351 x 19451
Available data:
- rawCount
- metaData
- spmetaData
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

### 2. Normalization, dimension reduction and cell clustering
We then perform data normalization, dimension reduction and cell clustering on the Object using the general functions developed for all data types. Since MERSCOPE output does not include cell clutering and umap dimensional reduction, we generate these results using internal function in ```CellScopes```.
```julia
merfish = cs.normalize_object(merfish)
merfish = cs.scale_object(merfish)
merfish = cs.find_variable_genes(merfish; nFeatures=500)
merfish = cs.run_pca(merfish;  method=:svd, pratio = 0.9, maxoutdim = 20)
merfish = cs.run_umap(merfish; min_dist=0.4, n_neighbors=20)
merfish = cs.run_clustering(merfish; res=0.00005)
```

### 3. Visiualization
We provide a set of functions to help visualize the cell and gene expression on the dimension reduction space and the spatial space.

#### 3.1 Cell annotation

##### a. Visualize the cell annotation on umap.
Note that some outlier data points were usually seen in the MERFISH data. We set the axis limit to trim away these outlier points (usually only several points so not affect the overall results). We are working on including a similar clipping approach (e.g. clip.range parameter) from [Seurat](https://satijalab.org/seurat/articles/spatial_vignette_2.html) in the normalization step to solve the outlier datapoint issue in umap.
```julia
cs.dim_plot(merfish; dim_type = "umap", marker_size = 4, 
    x_lims = (-15, 14), y_lims = (-14, 14))
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/merfish_umap.png" width="600"> 

##### b. Visualize the cell annotation on tissue.
```julia
colors =["blue", "red","magenta", "seagreen", "deepskyblue3", 
    "green2", "purple","darkgoldenrod1", "cyan", "dodgerblue"]
celltypes = string.(unique(merfish.metaData.cluster))
anno_color=Dict(celltypes .=> colors)
cs.sp_dim_plot(merfish, "cluster"; anno_color = anno_color,
    do_label = false,marker_size = 5, 
    canvas_size = (2500, 2000), do_legend=false
    )
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/merfish_tissue.png" width="600"> 

##### c. Select a field of view to visualize the cell annotation
```julia
cs.plot_fov(merfish, 20, 20)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/merfish_fov.png" width="600"> 

The grid plot above allows us to choose the specific region of interest for visualization. For instance, we can crop the view to focus on the hippocampus, as shown in the example below.

```julia
crop_fov = cs.subset_fov(merfish, [188, 189, 208, 209],20,20)
x1 = minimum(crop_fov.x)
x2 = maximum(crop_fov.x)
y1 = minimum(crop_fov.y)
y2 = maximum(crop_fov.y)
cs.sp_dim_plot(merfish, "cluster"; anno_color = anno_color,
    do_label = false,marker_size = 4, x_lims = (x1, x2),
    y_lims = (y1, y2),
    canvas_size = (500, 400), do_legend=false
    )
```

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/merfish_crop2.png" width="600"> 

The function below is for visualizing the cell type annotation directly on cell polygons.
<br>
```julia
cs.plot_cell_polygons(merfish, "cluster"; cell_colors=colors, stroke_color="gray40",
    x_lims=(x1, x2), y_lims = (y1, y2),stroke_width=0.5,
    width = 500, height = 400)
```

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/merfish_crop1.png" width="600"> 

#### 3.2 Gene expression visiualization
Plot genes on whole tissue:

```julia
cs.sp_feature_plot(merfish, ["COL1A1", "CDH1"]; 
    color_keys=["gray94", "lemonchiffon", "red"], 
    height=800, width=1000, marker_size = 4)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/merfish_gene_whole.png" width="600"> 

Plot genes on the selected region:

```julia
cs.plot_gene_polygons(merfish, ["COL1A1", "CDH1"]; x_lims = (x1, x2),y_lims = (y1, y2),
    width = 500, height = 400, color_keys=["gray94", "lemonchiffon", "red","darkred"])
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/merfish_gene_polygon.png" width="600"> 

#### 3.3 Save the object

To save all data in the ```MerfishObject```, simply run:

```julia
cs.save(merfish; filename="merfish.jld2")
```

