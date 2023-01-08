# CellScopes.jl <img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/logo.png" width="60" height="60"> <br>
```CellScopes.jl``` is a Julia-based toolkit for analyzing single cell data. It accepts a gene by cell count matrix as input and applies data normalization, dimensional reduction, cell clustering, and visualization techniques similar to those used in Seurat and Scanpy. Currently, ```CellScopes.jl``` only supports scRNA-seq data, but support for spatial transcriptomics and scATAC-seq is planned for future releases. This is our proposal for the package's development.

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/CellScopes.png" width="600"> <br>

### 1. Installation
To install ```CellScopes.jl```, you will need to have Julia 1.6 or higher installed. It is recommended to use Julia 1.7.3 or higher to avoid issues with dependencies. To install all of the necessary dependencies, run the following command line in Julia. Note that this will not install the unregisterd package ```Leiden.jl```, which you may need to install manually from the GitHub repository first.

```julia
using Pkg;
Pkg.add("https://github.com/bicycle1885/Leiden.jl") # Install the unregistered dependency Leiden.jl
Pkg.add("https://github.com/HaojiaWu/CellScopes.jl") # Install CellScopes.jl
```
### 2. Tutorial: PBMC 3K
This tutorial uses the pbmc3k dataset from 10x Genomics, which has been previously used by [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) for demo purpose. This will read in the data and create a RawCountObject that can be used as input for ```CellScopes.jl```.
#### 2.1 Download the pbmc3k data (in Terminal)
```bash
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar xvf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

#### 2.2 Read the data (in Julia)
The cells and genes can be filtered by setting the parameters ```min_gene``` and ```min_cell```, respectively.
```julia
import CellScopes as cs
raw_counts = cs.read_10x("filtered_gene_bc_matrices/hg19"; min_gene = 3);
```
This should create an object type called ```RawCountObject```.
```julia
raw_counts
#=
CellScopes.RawCountObject
Genes x Cells = 13100 x 2700
Available fields:
- count_mtx
- cell_name
- gene_name
=#
```
#### 2.3 Create a scRNAObject
We then create a ```scRNAObject``` using the count object above. The ```scRNAObject``` serves as a container to store all the data needed for and generated from the downstream analysis. The cells and genes can be further filtered by setting the parameters ```min_gene``` and ```min_cell```, respectively.
```julia
pbmc = cs.scRNAObject(raw_counts)
#=
scRNAObject in CellScopes.jl
Genes x Cells = 13100 x 2700
Available data:
- Raw count
=#
```
#### 2.4 Normalize the scRNAObject
We use a normalization method called global-scaling, which is similar to Seurat's "LogNormalize" method. This normalization method scales the feature expression measurements for each cell by the total expression, multiplies the result by a default scale factor of 10,000, and log-transforms the final value. The normalized values are stored as a ```NormCountObject```.
```julia
pbmc = cs.NormalizeObject(pbmc; scale_factor = 10000)
#=
scRNAObject in CellScopes.jl
Genes x Cells = 13100 x 2700
Available data:
- Raw count
- Normalized count
=#
```
We then use the ```ScaleObject``` function to scale and center the data.

```julia
pbmc = cs.ScaleObject(pbmc)
#=
scRNAObject in CellScopes.jl
Genes x Cells = 13100 x 2700
Available data:
- Raw count
- Normalized count
- Scaled count
=#
```

#### 2.5 Find variable genes
We use the ```vst``` approach implemented in the ```FindVariableFeatures``` function in Seurat or the ```pp.highly_variable_genes``` function in Scanpy to identify the variable genes. To standardize the counts, we use the following formula:
```math
Z_{ij} = \frac{x_{ij} - \bar{x}_i}{σ_i}
```
x<sub>ij</sub> is observed UMI, x̄ is the gene mean (rowMean) and σ<sub>i</sub> is the expected variance from the loess fit.

```julia
pbmc = cs.FindVariableGenes(pbmc)
#=
scRNAObject in CellScopes.jl
Genes x Cells = 13100 x 2700
Available data:
- Raw count
- Normalized count
- Scaled count
- Variable genes
=#
```

#### 2.6 Run principal component analysis (PCA).
Next, we perform PCA on the scaled data using only the previously identified variable genes as input. This is completed by the [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl) package.
```julia
pbmc = cs.RunPCA(pbmc;  method=:svd, pratio = 1, maxoutdim = 10)
#=
scRNAObject in CellScopes.jl
Genes x Cells = 13100 x 2700
Available data:
- Raw count
- Normalized count
- Scaled count
- Variable genes
- PCA data
=#
```
#### 2.7 Cell clustering.
We use a graph-based approach to identify the clusters. We first construct a KNN graph based on the significant components using the [NearestNeighborDescent.jl](https://github.com/dillondaudert/NearestNeighborDescent.jl) package. We then extract the KNN matrix from the graph and convert it into an adjacency matrix. This adjacent matrix is used as input for the [Leiden.jl](https://github.com/bicycle1885/Leiden.jl) package, which performs community detection. The entire process is implemented in the RunClustering function.

```julia
pbmc = cs.RunClustering(pbmc; res=0.015)
#=
scRNAObject in CellScopes.jl
Genes x Cells = 13100 x 2700
Available data:
- Raw count
- Normalized count
- Scaled count
- Variable genes
- Clustering data
- PCA data
=#
```
#### 2.8 Run UMAP or tSNE.
CellScapes.jl provides various non-linear dimensionality reduction techniques, including tSNE and UMAP, to allow for visualization and exploration of datasets. In the current version, UMAP is much faster than tSNE for large datasets, so it is highly recommended to use UMAP. We use the [TSne.jl](https://github.com/lejon/TSne.jl) and [UMAP.jl](https://github.com/dillondaudert/UMAP.jl) packages for tSNE and UMAP analysis, respectively.
```julia
pbmc = pbmc = cs.RunTSNE(pbmc; max_iter = 3000, perplexit = 50)
pbmc = cs.RunUMAP(pbmc; min_dist=0.4)
#=
scRNAObject in CellScopes.jl
Genes x Cells = 13100 x 2700
Available data:
- Raw count
- Normalized count
- Scaled count
- Variable genes
- Clustering data
- PCA data
- tSNE data
- UMAP data
=#
```
#### 2.9 Find markers.
```CellScopes.jl``` can help you find markers that define clusters through differential expression. Same as Seurat and Scanpy, we perform wilcoxon rank sum test on each pair of cell types to identify the differential genes. This is implemented by the [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) and [MultipleTesting.jl](https://github.com/juliangehring/MultipleTesting.jl)
```julia
markers = cs.FindMarkers(pbmc, "7", "6")
```
### 3. Data visualization
Inspired by Seurat and Scanpy, we utilize various methods to visualize cell annotations and gene expression. 
#### 3.1 Visualize cell annotaiton.
a. DimGraph on PCA
```julia
cs.DimGraph(pbmc; dim_type = "pca")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/pca.png" width="600"> <br>

b. DimGraph on tSNE
```julia
cs.DimGraph(pbmc; dim_type = "tsne")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/tsne.png" width="600"> <br>

c. DimGraph on UMAP
```julia
cs.DimGraph(pbmc; dim_type = "umap")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/umap.png" width="600"> <br>

d. DimGraph on selected cluster
```julia
cs.HighlightCells(pbmc, "6")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/highlight.png" width="500"> <br>

#### 3.2 Visualize gene expression.
a. GeneDimGraph
```julia
cs.GeneDimGraph(pbmc, ["CST3","IL32","CD79A"]; 
    order=false, marker_size = 4, 
    count_type ="norm", color_keys=("black","indianred1","red"))
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/featureplot.png" width="800"> <br>

b. GeneDimGraph (split by condition)
```julia
pbmc.metaData.fake_group = repeat(["group1","group2","group3"],900) # Create a fake condition
cs.GeneDimGraph(pbmc, ["CST3","IL32","CD79A"]; 
    order=false, marker_size = 7, 
    count_type ="norm", color_keys=("black","indianred1","red"), split_by="fake_group")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/split_by1.png" width="800"> <br>

c. GeneDotGraph 
```julia
cs.GeneDotGraph(pbmc, ["GZMB","GZMA", "CD3D","CD68","CD79A"], "cluster";
               count_type="norm",height=300, width=150,  expr_cutoff = 1, 
                cell_order=["1","5","4","3","8","2","7","6"])
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/dotgraph.png" width="300"> <br>

d. GeneDotGraph (split by condition)
```julia
cs.GeneDotGraph(pbmc, ["GZMB","GZMA", "CD3D","CD68","CD79A"], 
                "cluster"; split_by="fake_group",
                count_type="norm",height=300, width=150,  expr_cutoff = 1, 
                cell_order=["1","5","4","3","8","2","7","6"])
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/split_by2.png" width="600"> <br>

e. GeneVlnGraph
```julia
from = ["1","5","4","3","8","2","7","6"]
to = ["c1","c2","c3","c4","c5","c6","c7","c8"]
pbmc.metaData = mapvalues(pbmc.metaData, :cluster, :cluster2, from, to);
cs.GeneVlnGraph(pbmc, ["GZMB","GZMA", "CD3D","CD68","CD79A"]; group_by="cluster2",
height = 500,alpha=0.5, col_use = :tab10)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/violin.png" width="600"> <br>

### 4. Tutorial: MCA 400K cells
```CellScopes.jl``` can analyze atlas-scale single cell data as well. Below are some example codes to complete the analysis of the MCA dataset which contained ~400K cells. This take about 2.5 hours in a linux server with 256GB RAM and 16 cores.

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
rowSum(mtx::AbstractMatrix{<:Real}) = sum(mtx, dims=2)
@time gene_kept = (vec ∘ collect)(rowSum(counts).> 0.0);
genes = genes[gene_kept];
```
1.811749 seconds (2.46 M allocations: 131.292 MiB, 37.24% compilation time)
```julia
colSum(mtx::AbstractMatrix{<:Real}) = sum(mtx, dims=1)
@time cell_kept = (vec ∘ collect)(colSum(counts) .> 0.0)
cells = cells[cell_kept];
```
0.180891 seconds (18.55 k allocations: 4.606 MiB, 3.59% compilation time)
```julia
@time counts = counts[gene_kept, cell_kept];
```
4.641869 seconds (631.47 k allocations: 3.878 GiB, 28.75% gc time, 7.56% compilation time)
```julia
@time rawcount = cs.RawCountObject(counts, cells, genes);
```
0.011193 seconds (5.12 k allocations: 290.205 KiB, 99.64% compilation time)
```julia
@time mca = cs.scRNAObject(rawcount)
```
4.828527 seconds (1.61 M allocations: 3.954 GiB, 6.74% gc time, 16.00% compilation time)
```julia
@time mca = cs.NormalizeObject(mca; scale_factor = 10000)
```
15.791931 seconds (3.82 M allocations: 11.933 GiB, 20.05% gc time, 2.43% compilation time)
```julia
@time mca = cs.FindVariableGenes(mca)
```
360.900305 seconds (22.33 M allocations: 239.714 GiB, 2.77% gc time, 1.68% compilation time)
```julia
@time mca = cs.ScaleObject(mca; features = mca.varGene.var_gene)
```
169.916280 seconds (2.00 M allocations: 73.590 GiB, 6.25% gc time, 0.72% compilation time)
```julia
@time mca = cs.RunPCA(mca; maxoutdim = 30)
```
616.548600 seconds (4.26 M allocations: 36.830 GiB, 1.16% gc time, 0.94% compilation time)
```julia
@time mca = cs.RunUMAP(mca; reduce_dims = 30, min_dist = 0.6, n_neighbors=30)
```
2554.080855 seconds (63.06 M allocations: 22.751 GiB, 0.57% gc time, 0.14% compilation time)
```julia
@time mca = cs.RunClustering(mca; res=0.0001,n_neighbors=30) # To-do list: runtime optimization
```
3466.776112 seconds (2.30 M allocations: 1.198 TiB, 0.04% gc time, 0.02% compilation time)

```julia
cs.DimGraph(mca; marker_size =1, do_label=false, do_legend=false)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/umap2.png" width="600"> <br>
