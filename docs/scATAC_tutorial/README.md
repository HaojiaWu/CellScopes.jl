# Tutorial for scATAC-seq analysis with CellScopes.jl

The following tutorial illustrates a standard analysis for scATAC-seq data. The current version only supports the data output produced by 10x cellranger-atac. As it continues to evolve, ```CellScopes.jl``` will support more single-cell chromatin modalities. The scATAC analysis below was inspired heavily by the [Signac package](https://stuartlab.org/signac/) from Stuart's lab.

## 1 Tutorial: Mouse kidney 12K cells
This tutorial uses a mouse kidney dataset for demo purposes. To test the scATAC-seq functionalies in ```CellScopes.jl``` using your own working environment, please download an example data from the [10x website](https://www.10xgenomics.com/resources/datasets?query=&page=1&configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000). 

### 1.1 Load data
We provided an easy-to-use function ```read_atac``` to directly read the cellrange-atac output into Julia and construct a ```scATACObject```. ```scATACObject``` is a novel data structure in ```CellScopes.jl``` designed for storing the original and processed data to facilitate the downstream analysis. The only parameter required to construct a ```scATACObject``` is the path to the cellranger-atac output. Similar to the single cell RNA-seq analysis, you can also set the min_peak and min_cell parameters to filter cell cells and peaks, respectively.

```julia
import CellScopes as cs
atac_path = "/mnt/sdd/multiomics/atac_raw/m12hr_run/outs/"
@time atac_obj = cs.read_atac(atac_path; min_peak=2500)
```
This should read the peak count, peak annotation, and fragment file into Julia and construct a complete ```scATACObject```. Note that it will take about a while to complete this step since some files are big in size. The messages below shows that a ```scATACObject``` has been successfully constracted.
```julia
This step reads all information directly from cellranger-atac output for downstream analysis. It may take 10 - 15 mins to complete as certain files (e.g. the fragment file) can be large in size.
1/3 Reading peak count data...
1/3 Peak count was loaded!
2/3 Reading the peak annotation file...
2/3 Peak annotation was loaded!
3/3 Reading the fragment file...
3/3 Fragments were loaded!
scATACObject was successfully constructed!
542.737890 seconds (1.27 G allocations: 78.840 GiB, 2.45% gc time, 10.05% compilation time: 19% of which was recompilation)
scATACObject in CellScopes.jl
Peaks x Cells = 128626 x 12620
Available data:
- Raw count
- Fragment data
- Peak data
All fields:
- rawCount
- normCount
- scaleCount
- metaData
- varPeak
- dimReduction
- clustData
- peakAnno
- fragmentData
- activityData
- undefinedData
```
### 1.2 Preprocess the scATACObject
This step performs data normalization on the raw peak count. Shadowed by the normalization methods developed by Signac, the raw read count was normalized using the term frequency inverse document frequency approach (TF-IDF). Here we computes the TF-IDF matrix as log(TF x IDF) following the **Method 1** in the RunTFIDF function of Signac.
```julia
atac_obj = cs.run_tf_idf(atac_obj)
```

### 1.3 Find top features
We use a similar approach implemented in the ```FindTopFeatures``` function in Signac to identify the top peaks. Users can set the min_cutoff parameters to choose the n% features to be selected for top features. For example, "q35" means 35% of the top peaks to be selected for dimensional reduction analysis. 

```julia
atac_obj = cs.find_top_features(atac_obj; min_cutoff="q35")
```

### 1.4 Dimensional reduction
Next, we perform linear demansional reduction using latent semantic indexing (LSI) approach in Signac (TF-IDF + SVD are knowns as LSI). We then embed the cells into low dimensional space using the UMAP methods as described in scRNA-seq and used the same graph-based approach to partition the cell distnace matrix into cell clusters.
```julia
atac_obj = cs.run_svd(atac_obj; method=:svd, pratio=0.99, maxoutdim=20)
atac_obj = cs.run_umap(atac_obj; dims_use=2:20, min_dist=0.2)
atac_obj = cs.run_clustering(atac_obj; res=0.0015)
```

## 2 Data visualization
Inspired by Seurat and Scanpy, we utilize various methods to visualize cell annotations and gene expression. 
### 2.1 Visualize cell annotaiton.
a. Dim plot on PCA
```julia
cs.dim_plot(pbmc; dim_type = "pca", marker_size = 4)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/pca.png" width="600"> <br>

b. Dim plot on tSNE
```julia
cs.dim_plot(pbmc; dim_type = "tsne", marker_size = 4)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/tsne.png" width="600"> <br>

c. Dim plot on UMAP
```julia
cs.dim_plot(pbmc; dim_type = "umap", marker_size = 4)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/umap.png" width="600"> <br>

d. Dim plot on selected cluster
```julia
cs.highlight_cells(pbmc, "6")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/highlight.png" width="500"> <br>

### 2.2 Visualize gene expression.
a. Feature plot
```julia
cs.feature_plot(pbmc, ["CST3","IL32","CD79A"]; 
    order=false, marker_size = 4, 
    count_type ="norm", color_keys=("black","indianred1","red"))
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/featureplot.png" width="800"> <br>

b. Feature plot (split by condition)
```julia
pbmc.metaData.fake_group = repeat(["group1","group2","group3"],900) # Create a fake condition
cs.feature_plot(pbmc, ["CST3","IL32","CD79A"]; 
    order=false, marker_size = 7, 
    count_type ="norm", color_keys=("black","indianred1","red"), split_by="fake_group")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/split_by1.png" width="800"> <br>

c. Dot plot 
```julia
cs.dot_plot(pbmc, ["GZMB","GZMA", "CD3D","CD68","CD79A"], "cluster";
               count_type="norm",height=300, width=150,  expr_cutoff = 1, 
                cell_order=["1","5","4","3","8","2","7","6"])
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/dotgraph.png" width="300"> <br>

d. Dot plot (split by condition)
```julia
cs.dot_plot(pbmc, ["GZMB","GZMA", "CD3D","CD68","CD79A"], 
                "cluster"; split_by="fake_group",
                count_type="norm",height=300, width=150,  expr_cutoff = 1, 
                cell_order=["1","5","4","3","8","2","7","6"])
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/split_by2.png" width="600"> <br>

e. Violin plot
```julia
from = ["1","5","4","3","8","2","7","6"]
to = ["c1","c2","c3","c4","c5","c6","c7","c8"]
pbmc.metaData = cs.mapvalues(pbmc.metaData, :cluster, :cluster2, from, to);
cs.violin_plot(pbmc, ["GZMB","GZMA", "CD3D","CD68","CD79A"]; group_by="cluster2",
height = 500,alpha=0.5, col_use = :tab10)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/violin.png" width="600"> <br>

## 3. Tutorial: MCA 400K cells
```CellScopes.jl``` can analyze atlas-scale single cell data as well. Below are some example codes to complete the analysis of the [MCA dataset](https://figshare.com/articles/MCA_DGE_Data/5435866) which contains ~400K cells. This takes about 50 minutes in a linux server with 128GB RAM and 16 cores.

```julia
import CellScopes as cs
using MatrixMarket, CSV, DataFrames
using SparseArrays
```

```julia
counts = MatrixMarket.mmread("mca_mtx/matrix.mtx");
cells = CSV.File("mca_mtx/barcodes.tsv", header = false) |> DataFrame
cells = string.(cells.Column1)
genes = CSV.File("mca_mtx/genes.tsv", header = false) |> DataFrame
genes = string.(genes.Column2);
@time gene_kept = (vec ∘ collect)(cs.rowSum(counts).> 0.0);
genes = genes[gene_kept];
```
*1.811749 seconds (2.46 M allocations: 131.292 MiB, 37.24% compilation time)*
```julia
@time cell_kept = (vec ∘ collect)(cs.colSum(counts) .> 0.0)
cells = cells[cell_kept];
```
*0.180891 seconds (18.55 k allocations: 4.606 MiB, 3.59% compilation time)*
```julia
@time counts = counts[gene_kept, cell_kept];
```
*4.641869 seconds (631.47 k allocations: 3.878 GiB, 28.75% gc time, 7.56% compilation time)*
```julia
@time rawcount = cs.RawCountObject(counts, cells, genes);
```
*0.011193 seconds (5.12 k allocations: 290.205 KiB, 99.64% compilation time)*
```julia
@time mca = cs.scRNAObject(rawcount)
```
*4.828527 seconds (1.61 M allocations: 3.954 GiB, 6.74% gc time, 16.00% compilation time)*
```julia
@time mca = cs.normalize_object(mca; scale_factor = 10000)
```
*15.791931 seconds (3.82 M allocations: 11.933 GiB, 20.05% gc time, 2.43% compilation time)*
```julia
@time mca = cs.find_variable_genes(mca)
```
*217.251548 seconds (21.15 M allocations: 126.109 GiB, 4.52% gc time, 2.32% compilation time)*
```julia
@time mca = cs.scale_object(mca; features = mca.varGene.var_gene)
```
*147.311905 seconds (4.10 M allocations: 97.858 GiB, 8.31% gc time, 1.23% compilation time)*
```julia
@time mca = cs.run_pca(mca; maxoutdim = 30)
```
*236.710203 seconds (4.22 M allocations: 24.603 GiB, 0.52% gc time, 0.74% compilation time)*
```julia
@time mca = cs.run_umap(mca; reduce_dims = 30, min_dist = 0.6, n_neighbors=30, n_epochs=100)
```
*1075.675636 seconds (63.08 M allocations: 23.239 GiB, 1.64% gc time, 0.37% compilation time)*
```julia
@time mca = cs.run_clustering(mca; res=0.0001,n_neighbors=30) # To-do list: runtime optimization
```
*590.371976 seconds (40.33 M allocations: 1.199 TiB, 0.43% gc time, 0.13% compilation time)*

```julia
cs.dim_plot(mca; marker_size =1, do_label=false, do_legend=false)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/umap2.png" width="600"> <br>
