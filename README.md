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


