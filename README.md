# CellScopes.jl <img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/logo.png" width="60" height="60"> <br>
```CellScopes.jl``` is a Julia-based toolkit for analyzing single cell and spatial transcriptomic data. It performs data normalization, dimensional reduction, cell clustering, and visualization on various data types. Here is the structural design of ```CellScopes.jl``` in its current version.

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/CellScopes-version-1.png" width="600"> <br>

## 1. Installation
To install ```CellScopes.jl```, you will need to have Julia 1.6 or higher installed. It is recommended to use Julia 1.7.3 or higher to avoid issues with dependencies. To install all of the necessary dependencies, run the following command line in Julia. Note that this will not install the unregisterd package ```Leiden.jl```, which you may need to install manually from the GitHub repository first.

```julia
using Pkg;
Pkg.add("https://github.com/bicycle1885/Leiden.jl") # Install the unregistered dependency Leiden.jl
Pkg.add("https://github.com/HaojiaWu/CellScopes.jl") # Install CellScopes.jl
```
## 2. Tutorials

In this markdown page, We will show you how we use ```CellScopes.jl``` to analyze a dataset generated 10x Cartana (Xenium). CellScopes also supports analysis for scRNA-seq and Visium datasets. Please refer to the tutorials below: <br>
a. 10x Cartana: https://github.com/HaojiaWu/CellScopes.jl/tree/main/cartana_tutorial
<br>
b. scRNA-seq (PBMC 3K and MCA 400K cells): https://github.com/HaojiaWu/CellScopes.jl/tree/main/scRNA_tutorial
<br>
c. 10x Visium: https://github.com/HaojiaWu/CellScopes.jl/tree/main/visium_tutorial
<br>

