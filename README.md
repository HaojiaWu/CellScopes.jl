# CellScopes.jl

```CellScopes.jl``` is a Julia-based toolkit for analyzing single cell data. It accepts a gene by cell count matrix as input and applies data normalization, dimensional reduction, cell clustering, and visualization techniques similar to those used in Seurat and Scanpy. Currently, CellScopes.jl only supports scRNA-seq data, but support for spatial transcriptomics and scATAC-seq is planned for future releases. This is our proposal for the package's development.

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/CellScopes.png" width="600"> <br>

### 1. Installation
To install ```CellScopes.jl```, you will need to have Julia 1.6 or higher installed. It is recommended to use Julia 1.7.3 or higher to avoid issues with dependencies. To install all of the necessary dependencies, run the following command line in Julia. Note that this will not install the unregisterd package ```Leiden.jl```, which you may need to install manually from the GitHub repository first.

```julia
using Pkg;
Pkg.add("https://github.com/bicycle1885/Leiden.jl") # Install the unregistered dependency Leiden.jl
Pkg.add("https://github.com/HaojiaWu/CellScopes.jl") # Install CellScopes.jl
```
### 2. Tutorial for scRNA-seq analysis
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
This should create an object type called RawCountObject.
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
We then create a scRNAObject using the count object above. The scRNAObject serves as a container to store all the data needed for and generated from the downstream analysis. The cells and genes can be further filtered by setting the parameters ```min_gene``` and ```min_cell```, respectively.
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
We use a normalization method called global-scaling, which is similar to Seurat's "LogNormalize" method. This normalization method scales the feature expression measurements for each cell by the total expression, multiplies the result by a default scale factor of 10,000, and log-transforms the final value. The normalized values are stored as a NormCountObject.
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
#### 2.7 Find clusters.
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
pbmc = cs.RunTSNE(pbmc)
pbmc = cs.RunUMAP(pbmc)
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





