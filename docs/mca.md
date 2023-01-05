```julia
import CellScopes as cs
```


```julia
using MatrixMarket, CSV, DataFrames
using SparseArrays
```


```julia
counts = MatrixMarket.mmread("mca_mtx/matrix.mtx");
```


```julia
cells = CSV.File("mca_mtx/barcodes.tsv", header = false) |> DataFrame
cells = string.(cells.Column1)
genes = CSV.File("mca_mtx/genes.tsv", header = false) |> DataFrame
genes = string.(genes.Column2);
```


```julia
rowSum(mtx::AbstractMatrix{<:Real}) = sum(mtx, dims=2)
@time gene_kept = (vec ∘ collect)(rowSum(counts).> 0.0);
genes = genes[gene_kept];
```

      1.106898 seconds (2.42 M allocations: 127.514 MiB, 54.09% compilation time)



```julia
colSum(mtx::AbstractMatrix{<:Real}) = sum(mtx, dims=1)
@time cell_kept = (vec ∘ collect)(colSum(counts) .> 0.0)
cells = cells[cell_kept];
```

      0.179723 seconds (18.55 k allocations: 4.606 MiB, 3.34% compilation time)



```julia
@time counts = counts[gene_kept, cell_kept];
```

      3.100293 seconds (631.44 k allocations: 3.878 GiB, 2.58% gc time, 10.20% compilation time)



```julia
@time rawcount = cs.RawCountObject(counts, cells, genes);
```

      0.009917 seconds (5.12 k allocations: 290.771 KiB, 99.44% compilation time)



```julia
@time mca = cs.scRNAObject(rawcount)
```

      4.201256 seconds (1.42 M allocations: 3.940 GiB, 0.80% gc time, 11.55% compilation time)





    scRNAObject in CellScopes.jl




    Genes x Cells = 38884 x 405191
    Available data:
    - Raw count
    Available fields:
    - rawCount
    - normCount
    - scaleCount
    - metaData
    - varGene
    - dimReduction
    - clustData
    - undefinedData



```julia
@time mca = cs.NormalizeObject(mca; scale_factor = 10000) # To-do list: memory usage optimization
```

    208.063749 seconds (1.95 M allocations: 410.981 GiB, 3.19% gc time, 0.51% compilation time)





    scRNAObject in CellScopes.jl




    Genes x Cells = 38884 x 405191
    Available data:
    - Raw count
    - Normalized count
    Available fields:
    - rawCount
    - normCount
    - scaleCount
    - metaData
    - varGene
    - dimReduction
    - clustData
    - undefinedData



```julia
@time mca = cs.FindVariableGenes(mca) # To-do list: memory usage optimization
```

    360.900305 seconds (22.33 M allocations: 239.714 GiB, 2.77% gc time, 1.68% compilation time)





    scRNAObject in CellScopes.jl




    Genes x Cells = 38884 x 405191
    Available data:
    - Raw count
    - Normalized count
    - Variable genes
    Available fields:
    - rawCount
    - normCount
    - scaleCount
    - metaData
    - varGene
    - dimReduction
    - clustData
    - undefinedData



```julia
@time mca = cs.ScaleObject(mca; features = pbmc.varGene.var_gene) # To-do list: runtime optimization
```

    866.051116 seconds (1.63 G allocations: 63.882 GiB, 12.91% gc time, 0.33% compilation time)





    scRNAObject in CellScopes.jl




    Genes x Cells = 38884 x 405191
    Available data:
    - Raw count
    - Normalized count
    - Scaled count
    - Variable genes
    Available fields:
    - rawCount
    - normCount
    - scaleCount
    - metaData
    - varGene
    - dimReduction
    - clustData
    - undefinedData



```julia
@time mca = cs.RunPCA(mca; maxoutdim = 75) # To-do list: runtime optimization
```

    616.548600 seconds (4.26 M allocations: 36.830 GiB, 1.16% gc time, 0.94% compilation time)





    scRNAObject in CellScopes.jl




    Genes x Cells = 38884 x 405191
    Available data:
    - Raw count
    - Normalized count
    - Scaled count
    - Variable genes
    - PCA data
    Available fields:
    - rawCount
    - normCount
    - scaleCount
    - metaData
    - varGene
    - dimReduction
    - clustData
    - undefinedData



```julia
@time mca = cs.RunClustering(mca; reduce_dims = 75, res=0.015) # To-do list: runtime optimization
```

    3390.240210 seconds (23.96 M allocations: 1.204 TiB, 0.27% gc time, 0.02% compilation time)





    scRNAObject in CellScopes.jl




    Genes x Cells = 38884 x 405191
    Available data:
    - Raw count
    - Normalized count
    - Scaled count
    - Variable genes
    - Clustering data
    - PCA data
    Available fields:
    - rawCount
    - normCount
    - scaleCount
    - metaData
    - varGene
    - dimReduction
    - clustData
    - undefinedData



```julia
@time mca = cs.RunUMAP(mca; reduce_dims = 75)
```

    3463.918168 seconds (64.92 M allocations: 63.651 GiB, 0.61% gc time, 0.06% compilation time)





    scRNAObject in CellScopes.jl




    Genes x Cells = 38884 x 405191
    Available data:
    - Raw count
    - Normalized count
    - Scaled count
    - Variable genes
    - Clustering data
    - PCA data
    - UMAP data
    Available fields:
    - rawCount
    - normCount
    - scaleCount
    - metaData
    - varGene
    - dimReduction
    - clustData
    - undefinedData



```julia
using JLD2
```


```julia
cs.DimGraph(pbmc)
```


```julia

```
