abstract type AbstractCount end

mutable struct RawCountObject <: AbstractCount
    count_mtx::Matrix{Float64}
    cell_name::Vector{String}
    gene_name::Vector{String}
    function RawCountObject(count_mtx, cell_name, gene_name)
        if length(cell_name) !== size(count_mtx)[2]
            error("Total cell number mismatches between count matrix and cell ID!")
        end
        if length(gene_name) !== size(count_mtx)[1]
            error("Total gene number mismatches between count matrix and gene names!")
        end
        raw_count = new(count_mtx, cell_name, gene_name)
        return raw_count
    end
end

mutable struct NormCountObject <: AbstractCount
    gene_name::Vector{String}
    cell_name::Vector{String}
    count_mtx::Matrix{Float64}
    scale_factor::Union{Int64, Float64}
    norm_method::String
    pseudocount::Union{Int64, Float64}
    function NormCountObject(count_mtx, cell_name, gene_name, scale_factor, norm_method, pseudocount)
        if length(cell_name) !== size(count_mtx)[2]
            error("Total cell number mismatches between count matrix and cell ID!")
        end
        if length(gene_name) !== size(count_mtx)[1]
            error("Total gene number mismatches between count matrix and gene names!")
        end
        normcountobj = new(count_mtx, cell_name, gene_name, scale_factor, norm_method, pseudocount)
        return normcountobj
    end
end

