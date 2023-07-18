# Analysis and visualization of SlideSeq data with CellScopes.jl
Slide-seq is a sequencing-based spatial transcriptomics technique jointly developed by Macosko lab and Chen lab.  In this tutorial, we'll be exploring how to utilize```CellScopes.jl``` for the analysis of Slide-seq data derived from a [mouse brain](https://singlecell.broadinstitute.org/single_cell/study/SCP815/highly-sensitive-spatial-transcriptomics-at-near-cellular-resolution-with-slide-seqv2#study-download).

## 1. Read the SlideSeq data

To read the slide-seq data in Julia, just need to feed CellScopes with the gene expression matrix and the beads coordinates files.  

```julia
import CellScopes as cs
slide = cs.read_slideseq("Puck_200115_08_bead_locations.csv", "Puck_200115_08.digital_expression.txt.gz"; min_gene=100, min_cell = 3)```
```
```
This should create an data type called ```SlideObject```.
SlideseqObject in CellScopes.jl
Genes x Cells = 19979 x 41674
Available data:
- rawCount
- metaData
- spmetaData
All fields:
- rawCount
- normCount
- scaleCount
- metaData
- spmetaData
- varGene
- dimReduction
- clustData
```

### 2. Normalization, dimension reduction and cell clustering
We then perform data normalization, dimension reduction and cell clustering on the SlideSeq Object.
```julia
slide = cs.normalize_object(slide; scale_factor = 10000)
slide = cs.scale_object(slide)
slide = cs.find_variable_genes(slide; nFeatures = 1000)
slide = cs.run_pca(slide;  method=:svd, pratio = 1, maxoutdim = 10)
slide = cs.run_umap(slide; min_dist=0.1, n_neighbors=20)
slide = cs.run_clustering(slide; res=0.0005, n_neighbors=20)
```
### 3. Visiualization

#### 3.1 Cell type visualization
##### a. Visualize the cell clustering on UMAP.
```julia
cs.dim_plot(slide; dim_type ="umap", marker_size=2)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/slide_umap.png" width="600"> 

##### b. Visualize the cell type distribution on tissue
```julia
colors =["plum", "red","magenta", "darkred", "deepskyblue3", 
    "green2", "purple","darkgoldenrod1", "cyan", "blue", 
    "yellow", "pink","dodgerblue"]
celltypes = string.(unique(slide.metaData.cluster))
anno_color=Dict(celltypes .=> colors)

cs.sp_dim_plot(slide, "cluster"; anno_color=anno_color,
    marker_size = 2, canvas_size = (600,500), 
    do_label=false)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/slide_tissue.png" width="600"> 

##### c. Highlight the cell type of interest

#### 3.2 Gene expression visiualization

```julia
cs.sp_feature_plot(slide, ["Kctd8", "Cplx3", "Cpne6", "Neurod6", "Prox1", "Necab1"]; 
    marker_size = 3)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/slide_genes.png" width="600"> 
