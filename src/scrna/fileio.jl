function read_10x(tenx_dir::String; 
    version::String="v2", 
    min_gene::Real = 0.0, 
    min_cell::Real = 0.0
)
    if version === "v2"
        counts = MatrixMarket.mmread(tenx_dir * "/matrix.mtx")
        #counts = Matrix(counts)
        cells = CSV.File(tenx_dir * "/barcodes.tsv", header = false) |> DataFrame
        cells = string.(cells.Column1)
        genes = CSV.File(tenx_dir * "/genes.tsv", header = false) |> DataFrame
        genes = string.(genes.Column2)
        gene_kept = (vec ∘ collect)(rowSum(counts).> min_cell)
        genes = genes[gene_kept]
        cell_kept = (vec ∘ collect)(colSum(counts) .> min_gene)
        cells = cells[cell_kept]
        counts = counts[gene_kept, cell_kept]
        gene_kept, gene_removed = check_duplicates(genes)
        gene_removed = collect(values(gene_removed))
        counts = counts[Not(gene_removed), :]
        rawcount = RawCountObject(counts, cells, gene_kept)
        return rawcount
    end
end

function SaveObj(sc_obj::AbstractCellScope; key::String = "CSObject", filename::String = "cs_obj.jld2")
    JLD2.save(filename, key, sc_obj)
end

