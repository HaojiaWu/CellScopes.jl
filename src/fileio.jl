function read_10x(tenx_dir::String; 
        version::String="v2", 
        min_gene::Union{Float64, Int64} = 0.0, 
        min_cell::Union{Float64, Int64}=0.0
    )
    if version === "v2"
        counts = MatrixMarket.mmread(tenx_dir * "/matrix.mtx")
        counts = Matrix(counts)
        counts = DataFrame(counts, :auto)
        cells = CSV.File(tenx_dir * "/barcodes.tsv", header = false) |> DataFrame
        cells = string.(cells.Column1)
        genes = CSV.File(tenx_dir * "/genes.tsv", header = false) |> DataFrame
        genes = string.(genes.Column2)
        gene_kept = rowSums(counts, 1:nrow(counts), 1:ncol(counts)) .> min_gene
        genes = genes[gene_kept]
        cell_kept = colSums(counts, 1:nrow(counts), 1:ncol(counts)) .> min_cell
        cells = cells[cell_kept]
        counts = counts[gene_kept, cell_kept]
        rename!(counts, cells)
        counts.gene=genes
        counts = unique(counts, :gene)
        select!(counts, :gene, Not(:gene))
        return counts
    end
end
