# AbstractCellScope: toppest level type
abstract type AbstractCellScope end

# AbstractSingleCell: type for scRNA-seq
abstract type AbstractSingleCell <: AbstractCellScope end
abstract type AbstractCount <: AbstractSingleCell end
abstract type AbstractDimReduction <: AbstractSingleCell end

# scRNA-seq object constructors
mutable struct RawCountObject <: AbstractCount
    count_mtx::AbstractMatrix{<:Real}
    cell_name::Vector{String}
    gene_name::Vector{String}
    function RawCountObject(count_mtx, cell_name, gene_name)
        if length(cell_name) !== size(count_mtx)[2]
            error("The total number of cells in the count matrix and the cell ID do not match!")
        end
        if length(gene_name) !== size(count_mtx)[1]
            error("The total number of genes in the count matrix and the gene name do not match!")
        end
        raw_count = new(count_mtx, cell_name, gene_name)
        return raw_count
    end
end

mutable struct NormCountObject <: AbstractCount
    count_mtx::AbstractMatrix{<:Real}
    cell_name::Vector{String}
    gene_name::Vector{String}
    scale_factor::Real
    norm_method::String
    pseudocount::Real
    function NormCountObject(count_mtx, cell_name, gene_name, scale_factor, norm_method, pseudocount)
        if length(cell_name) !== size(count_mtx)[2]
            error("The total number of cells in the count matrix and the cell ID do not match!")
        end
        if length(gene_name) !== size(count_mtx)[1]
            error("The total number of genes in the count matrix and the gene name do not match!")
        end
        normcountobj = new(count_mtx, cell_name, gene_name, scale_factor, norm_method, pseudocount)
        return normcountobj
    end
end

mutable struct ScaleCountObject <: AbstractCount
    count_mtx::AbstractMatrix{<:Real}
    cell_name::Vector{String}
    gene_name::Vector{String}
    do_scale::Bool
    do_center::Bool
    scale_max::Real
    function ScaleCountObject(count_mtx, cell_name, gene_name, do_scale, do_center,scale_max)
        if length(cell_name) !== size(count_mtx)[2]
            error("The total number of cells in the count matrix and the cell ID do not match!")
        end
        if length(gene_name) !== size(count_mtx)[1]
            error("The total number of genes in the count matrix and the gene name do not match!")
        end
        scale_obj = new(count_mtx, cell_name, gene_name, do_scale, do_center,scale_max)
        return scale_obj
    end
end

mutable struct PCAObject <: AbstractDimReduction
    cell_embedding::AbstractMatrix{<:Real}
    pca_model::Union{PCA{Float64}, Nothing}
    percent_var::Union{Vector{Float64}, Nothing}
    key::Union{String, Nothing}
    method::Union{Symbol, Nothing}
    pratio::Union{Real, Nothing}
    maxoutdim::Union{Int64, Nothing}
    PCAObject(cell_embedding, pca_model, percent_var, key, method, pratio, maxoutdim) = new(cell_embedding, pca_model, percent_var, key, method, pratio, maxoutdim)
end

mutable struct tSNEObject <: AbstractDimReduction
    cell_embedding::AbstractMatrix{<:Real}
    key::Union{String, Nothing}
    ndim::Union{Int64, Nothing}
    reduce_dims::Union{Int64, Nothing}
    max_iter::Union{Int64, Nothing}
    perplexit::Union{Int64, Nothing}
    tSNEObject(cell_embedding, key, ndim, reduce_dims, max_iter, perplexit) = new(cell_embedding, key, ndim, reduce_dims, max_iter, perplexit)
end

mutable struct UMAPObject <: AbstractDimReduction
    cell_embedding::AbstractMatrix{<:Real}
    key::Union{String, Nothing}
    n_components::Union{Int64, Nothing}
    n_dimensions::Union{Int64, Nothing}
    n_neighbors::Union{Int64, Nothing}
    metric::Union{String, Nothing}
    min_dist::Union{Real, Nothing}
    knn_data::Union{AbstractMatrix{<:Real}, Nothing}
    UMAPObject(cell_embedding, key, n_components, n_dimensions, n_neighbors, metric, min_dist, knn_data) = new(cell_embedding, key, n_components, n_dimensions, n_neighbors, metric, min_dist, knn_data)    
end

mutable struct ReductionObject <: AbstractDimReduction
    pca::Union{PCAObject, Nothing}
    tsne::Union{tSNEObject, Nothing}
    umap::Union{UMAPObject, Nothing}
    ReductionObject(pca, tsne, umap) = new(pca, tsne, umap)
end

mutable struct ClusteringObject <: AbstractSingleCell
    clustering::DataFrame
    metric::Union{String, Nothing}
    adj_mat::Union{AbstractMatrix{<:Real}, Nothing}
    ledein_res::Union{NamedTuple{(:quality, :partition), Tuple{Float64, Vector{Vector{Int64}}}}, Nothing}
    resolution::Union{Real, Nothing}
    ClusteringObject(clustering, metric, adj_mat, ledein_res, resolution)=new(clustering, metric, adj_mat, ledein_res, resolution)
end

mutable struct VariableGeneObject <: AbstractSingleCell
    var_gene::Vector{String}
    vst_data::DataFrame
    VariableGeneObject(var_gene, vst_data) = new(var_gene, vst_data)
end

mutable struct UndefinedObject <: AbstractSingleCell
    uns1::Any
    uns2::Any
    uns3::Any
    uns4::Any
    uns5::Any
    UndefinedObject(uns1, uns2, uns3, uns4, uns5) = new(uns1, uns2, uns3, uns4, uns5)
end

mutable struct scRNAObject <: AbstractSingleCell
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    undefinedData::Union{UndefinedObject, Nothing}
    function scRNAObject(raw_count::RawCountObject; 
            meta_data::Union{DataFrame, Nothing} = nothing,
            #min_gene::Int64 = 0,
            #min_cell::Int64 = 0,
            prefix::Union{String, Nothing} = nothing,
            postfix::Union{String, Nothing} = nothing)
        count_mat = raw_count.count_mtx
        genes = raw_count.gene_name
        cells = raw_count.cell_name
        #gene_kept = (vec ∘ collect)(rowSum(count_mat).> min_cell)
        #genes = genes[gene_kept]
        #cell_kept = (vec ∘ collect)(colSum(count_mat) .> min_gene)
        #cells = cells[cell_kept]
        #count_mat = count_mat[gene_kept, cell_kept]
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = raw_count.cell_name, nFeatures=nFeatures, nGenes = nGenes)
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
        scRNA_obj = new(count_obj)
        scRNA_obj.metaData = meta_data
        return scRNA_obj
    end
end