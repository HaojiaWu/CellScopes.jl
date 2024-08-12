# Analysis and visualization of STARmap data with CellScopes.jl
[STARmap](https://www.science.org/doi/10.1126/science.aat5691) is an imaging-based spatial technique that enables high-resolution, 3D spatial profiling of gene expression within individual cells in their native tissue environment. Below is a simple tutorial to apply ```CellScopes``` STARmap.


## 1. Download STARmap data and perform cell segmentation with Baysor
We followed the Baysor tutorial to segment the cells for STARmap.
```julia
## download data
;wget -P data http://pklab.med.harvard.edu/viktor/baysor/starmap/molecules.csv
;wget -P data http://pklab.med.harvard.edu/viktor/baysor/starmap/segmentation.tiff
## cell segmentation
;export JULIA_NUM_THREADS=18
;baysor  run -x x -y y -z z -g gene -s 50 -m 50 -o ./starmap -p --save-polygons=geojson ./data/molecules.csv ./data/segmentation.tiff
```
## 2. Read all files and prepare STARmap object

```julia
import CellScopes as cs
using CSV, DataFrames
using SparseArrays
molecules =  DataFrame(CSV.File("/mnt/sdd/starmap/segmentation.csv"))
count_df =  DataFrame(CSV.File("/mnt/sdd/starmap/segmentation_counts.tsv"))
cells =  DataFrame(CSV.File("/mnt/sdd/starmap/segmentation_cell_stats.csv"))
gene_name = count_df.gene
cell_name = string.(names(count_df)[2:end])
count_df = count_df[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = cs.RawCountObject(count_df, cell_name, gene_name)
molecules.cell = string.(molecules.cell)
cells.cell = string.(cells.cell)
brain = cs.starMapObject(molecules, cells, raw_count;
    prefix = "brain", min_gene = 0, min_cell = 3)
```
This will create a STARmap object
```
Adding prefix brain to all cells...
STARmapObject in CellScopes.jl
Genes x Cells = 1020 x 5385
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
- polygonData
```

### 2. Normalization, dimension reduction and cell clustering
Then we can perform normalization, dimension reduction and cell clustering using the internal methods in ```CellScopes```.
```julia
brain = cs.normalize_object(brain)
brain = cs.scale_object(brain)
brain = cs.find_variable_genes(brain; nFeatures=1020)
brain = cs.run_pca(brain;  method=:svd, pratio = 1, maxoutdim = 20)
brain = cs.run_umap(brain; min_dist=0.2, n_neighbors=30)
brain = cs.run_clustering(brain; res=0.008, n_neighbors=30)
```

### 3. Visiualization
We provide a set of functions to help visualize the cell and gene expression on the dimension reduction space and the spatial space.

#### 3.1 Cell annotation

##### a. Visualize the cell annotation on umap.

```julia
cs.dim_plot(brain; dim_type = "umap", marker_size = 4)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/starmap_clustering.png" width="600"> 

##### b. Visualize the cell annotation on tissue.
```julia
cs.sp_dim_plot(brain, "cluster"; #anno_color = anno_color,
    do_label = false,marker_size = 8, 
    width=800, height=300, do_legend=false
    )
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/starmap_cell_tissue.png" width="600"> 


The function below is for visualizing the cell type annotation directly on cell polygons.
<br>
```julia
cs.plot_cell_polygons(brain, "cluster"; stroke_color="black",
    width = 800, height = 300)
```

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/starmap_cell_polygons.png" width="600"> 

#### 3.2 Gene expression visiualization
Plot genes on whole tissue:

```julia
cs.sp_feature_plot(brain, ["Pcp4","Slc17a7"]; 
    color_keys=["gray90", "dodgerblue1", "blue"], 
    height=300, width=800, marker_size = 8)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/starmap_gene_plot.png" width="600"> 

Plot genes on cell polygons:

```julia
cs.plot_gene_polygons(brain, ["Pcp4"]; color_keys=["gray96","lemonchiffon","red","darkred"],width = 800, height = 300,stroke_width=0.2)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/starmap_gene_polygon.png" width="600"> 

#### 3.3 Save the object

```julia
cs.save(brain, filename="brain_starmap.jld2")
```
