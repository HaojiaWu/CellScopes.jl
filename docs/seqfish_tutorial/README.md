# Analysis and visualization of seqFISH data with CellScopes.jl
[seqFISH (Sequential Fluorescence In Situ Hybridization)](https://spatial.caltech.edu/seqfish/) is an advanced microscopy technique used to visualize the spatial arrangement of individual mRNA molecules within cells. This tutorial walks you through how to use ```CellScopes``` to analyze [a embryo seqFISH dataset](https://www.nature.com/articles/s41587-021-01006-2).

### 1. Read input files and create a seqFISH object
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

### 2. Read input files and create a seqFISH object
At this point, we can transition back to Julia to continue processing the input files in the appropriate format that CellScopes can further process. The seqFISH data contain three embryo dataset. We only use one of them (embryo1) in this tutorial.

#### a. Read all files and reformat them. Note that only "embryo1" was used for downstream analysis.
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
#### b. Create a seqFISH object

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

### 3. Normalization, dimension reduction and cell clustering
The seqFISH object can be further processed using the general functions in ```CellScopes``` for normalization, PCA, cell clustering and low dimentional cell projection.
```julia
embryo1 = cs.normalize_object(embryo1)
embryo1 = cs.scale_object(embryo1)
embryo1 = cs.find_variable_genes(embryo1; nFeatures=351)
embryo1 = cs.run_pca(embryo1;  method=:svd, pratio = 0.9, maxoutdim = 20)
embryo1 = cs.run_umap(embryo1; min_dist=0.4, n_neighbors=20)
embryo1 = cs.run_clustering(embryo1; res=0.005);
```

### 4. Create cell polygons for visualization
We use Baysor(https://github.com/kharchenkolab/Baysor) to draw cell boundary polygons to mimic cell segmentation. We then map the polygons to the closest cells based on the Euclidean distance. Finally, we create a gene-by-polygon count matrix for gene visualization purpose.
```julia
import Baysor as B
em1_transcript.cell = parse.(Int64, em1_transcript.cell)
polygons = B.boundary_polygons(
    em1_transcript, em1_transcript.cell, grid_step=0.1, bandwidth=0.2
)
embryo1.polygonData = polygons;
embryo1 = cs.polygons_cell_mapping(embryo1)
embryo1 = cs.generate_polygon_counts(embryo1)
```

### 5. Visiualization

#### 5.1 Cell annotation

##### a. Visualize the cell annotation on umap.

```julia
colors =[ "mediumorchid4","darkseagreen1", "palegreen1","blue", "lightsalmon", "lightyellow1", 
     "slategray4", "yellow3", "darkseagreen1", "olivedrab2", "green", "slateblue1", 
    "#978bf0", "black", "#a7d64e", "#788c3b", "#57e24c", "dodgerblue", 
    "purple", "#db3c18", "cyan", "#f87197", "#ff0087", "#42ebd1", "#28997e", "pink", 
     "red", "fuchsia", "yellow","deepskyblue3", "#aa2afa"];

cs.dim_plot(embryo1; dim_type = "umap", marker_size = 4, anno_color= anno_color)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seqfish_umap2.png" width="600"> 

##### b. Visualize the cell annotation on tissue.
```julia
celltypes = string.(unique(embryo1.metaData.cluster))
anno_color=Dict(celltypes .=> colors)

cs.sp_dim_plot(embryo1, "cluster"; anno_color = anno_color,
    do_label = false,marker_size = 7, do_legend=false
    )
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seqfish_sp.png" width="600"> 

##### c. Visualize cell boundary in tissue

```julia
cs.plot_cell_polygons(embryo1, "cluster";  cell_colors= colors,
    stroke_color="gray70", stroke_width=0.5, 
    width=800, height=900)
```

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seqfish_cellpoly.png" width="600"> 

#### 5.2 Gene expression visiualization
Plot genes on whole tissue:

```julia
cs.sp_feature_plot(embryo1, ["Pou3f1","Popdc2", "Six3"]; 
    color_keys=["gray94", "lemonchiffon", "red","darkred"], marker_size=4
)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seqfish_geneplot.png" width="600"> 

Plot genes on cell polygons:

```julia
cs.plot_gene_polygons(embryo1, ["Popdc2"];
    width = 800, height = 1000, 
    color_keys=["midnightblue", "green", "orange","red"])
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/seqfish_gene.png" width="600"> 

#### 5.3 Save the object

To save all data in the ```seqFISHObject```, simply run:

```julia
cs.save(embryo1;filename="embryo1.jld2")
```

