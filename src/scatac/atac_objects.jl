
abstract type AbstractATAC <: AbstractCellScope end

mutable struct GeneActivityObject <: AbstractCount
    peak_anno::Union{DataFrame, Nothing}
    count_mtx::Union{AbstractMatrix{<:Real}, Nothing}
    cell_name::Union{Vector{String}, Nothing}
    gene_name::Union{Vector{String}, Nothing}
    GeneActivityObject(peak_anno, count_mtx, cell_name, gene_name) = new(peak_anno, count_mtx, cell_name, gene_name)
end

mutable struct FragmentObject <: AbstractATAC
    fragment::Union{DataFrame, Nothing}
    genecode::Union{DataFrame, Nothing}
    FragmentObject(fragment_data, gene_code) = new(fragment_data, gene_code)
end

mutable struct scATACObject <: AbstractATAC
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    varPeak::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    peakAnno::Union{DataFrame, Nothing}
    fragmentData::Union{FragmentObject, Nothing}
    activityData::Union{GeneActivityObject, Nothing}
    undefinedData::Union{UndefinedObject, Nothing}
    function scATACObject(raw_count::RawCountObject; 
            meta_data::Union{DataFrame, Nothing} = nothing,
            min_gene::Int64 = 0,
            min_cell::Int64 = 0,
            prefix::Union{String, Nothing} = nothing,
            postfix::Union{String, Nothing} = nothing)
        count_mat = raw_count.count_mtx
        genes = raw_count.gene_name
        cells = raw_count.cell_name
        gene_kept = (vec ∘ collect)(rowSum(count_mat).> min_cell)
        genes = genes[gene_kept]
        cell_kept = (vec ∘ collect)(colSum(count_mat) .> min_gene)
        cells = cells[cell_kept]
        count_mat = count_mat[gene_kept, cell_kept]
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = raw_count.cell_name, nFeatures=nFeatures, nGenes=nGenes)
        end
        if prefix !== nothing
            println("Adding prefix " * prefix * " to all cells...")
            cellnames = prefix * "_" .* raw_count.cell_name
        end
        if postfix !== nothing
            println("Adding postfix " * postfix * " to all cells...")
            cellnames = raw_count.cell_name .* "_" .* postfix
        end
        count_obj = RawCountObject(count_mat, cells, genes)
        scATAC_obj = new(count_obj)
        scATAC_obj.metaData = meta_data
        return scATAC_obj
    end
end
