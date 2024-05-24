# AbstractCellScope: toppest level type
abstract type AbstractCellScope end

# AbstractSingleCell: type for scRNA-seq
abstract type AbstractSingleCell <: AbstractCellScope end
abstract type AbstractCount <: AbstractSingleCell end
abstract type AbstractDimReduction <: AbstractSingleCell end
abstract type AbstractHarmony <: AbstractSingleCell end

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
            min_gene::Int64 = 0,
            min_cell::Int64 = 0,
            prefix::Union{String, Nothing} = nothing,
            postfix::Union{String, Nothing} = nothing)
        count_mat = raw_count.count_mtx
        genes = raw_count.gene_name
        cells = raw_count.cell_name
        count_mat, genes, cells = subset_matrix(count_mat, genes, cells, min_gene, min_cell)
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

mutable struct HarmonyObject <: AbstractHarmony
    Z_corr::Array{Float64, 2}
    Z_orig::Array{Float64, 2}
    Z_cos::Array{Float64, 2}
    Phi::Array{Float64, 2}
    Phi_moe::Array{Float64, 2}
    N::Int
    Pr_b::Vector{Float64}
    B::Int
    d::Int
    window_size::Int
    epsilon_kmeans::Float64
    epsilon_harmony::Float64
    lamb::Array{Float64, 2}
    sigma::Vector{Float64}
    sigma_prior::Vector{Float64}
    block_size::Float64
    K::Int
    max_iter_harmony::Int
    max_iter_kmeans::Int
    verbose::Bool
    theta::Vector{Float64}
    objective_harmony::Vector{Float64}
    objective_kmeans::Vector{Float64}
    objective_kmeans_dist::Vector{Float64}
    objective_kmeans_entropy::Vector{Float64}
    objective_kmeans_cross::Vector{Float64}
    kmeans_rounds::Vector{Int}
    scale_dist::Array{Float64, 2}
    dist_mat::Array{Float64, 2}
    O::Array{Float64, 2}
    E::Array{Float64, 2}
    W::Array{Float64, 2}
    Phi_Rk::Array{Float64, 2}
    Y::Array{Float64, 2}
    R::Array{Float64, 2}

    function HarmonyObject(
        data_mat::Array{Float64, 2},
        meta_data::DataFrame,
        vars_use;
        theta::Union{Nothing, Float64, Int, Array{Float64}} = nothing,
        lamb::Union{Nothing, Float64, Int, Array{Float64}} = nothing,
        sigma::Union{Nothing, Float64, Int} = 0.1,
        nclust::Union{Nothing, Int} = nothing,
        tau::Float64 = 0.0,
        block_size::Float64 = 0.05,
        max_iter_harmony::Int = 10,
        max_iter_kmeans::Int = 20,
        epsilon_cluster::Float64 = 1e-5,
        epsilon_harmony::Float64 = 1e-4,
        verbose::Bool = true,
        random_state::Int = 123
    )
        N = size(meta_data, 1)
        if size(data_mat, 2) != N
            data_mat = transpose(data_mat)
        end
        @assert size(data_mat, 2) == N "The number of cells in the metadata must match the number of cells in the count matrix."
        if nclust === nothing
            nclust = min(round(Int, N / 30.0), 100)
        end
        K = nclust
        if typeof(sigma) == Float64 && nclust > 1
            sigma = fill(sigma, nclust)
        end       
        phi = get_dummies(meta_data, [vars_use])
        phi = phi[:, size(meta_data)[2]:end]
        phi = transpose(Matrix(phi))
        phi_n = length(unique(meta_data[!, vars_use]))
        if theta === nothing
            theta = fill(2.0, phi_n)
        else
            theta = fill(theta, phi_n)
        end
        @assert length(theta) == phi_n 
        if lamb === nothing
            lamb = fill(1.0, phi_n)
        else
            lamb = fill(lamb, phi_n)
        end
        @assert length(lamb) == phi_n
        N_b = sum(phi, dims=2)
        Pr_b = vec(N_b / N)
        if tau > 0
            theta = theta .* (1 .- exp.(-(N_b / (nclust * tau)) .^ 2))
        end
        lamb_mat = diagm(0 => vcat(0.0, lamb))
        phi_moe = vcat(ones(1, N), phi)
        Random.seed!(random_state)                                        
        Z_corr = deepcopy(data_mat)
        Z_orig = deepcopy(data_mat)
        Z_cos = Z_orig' ./ maximum(Z_orig', dims=2)
        Z_cos = Matrix(Z_cos')
        Z_cos = Z_cos ./ mapslices(x -> norm(x, 2), Z_cos, dims=1)   
        N = size(Z_corr, 2)
        B = size(phi, 1)
        d = size(Z_corr, 1)
        window_size = 3
        objective_harmony = Float64[]
        objective_kmeans = Float64[]
        objective_kmeans_dist = Float64[]
        objective_kmeans_entropy = Float64[]
        objective_kmeans_cross = Float64[]
        kmeans_rounds = Int[]
        scale_dist = zeros(K, N)
        dist_mat = zeros(K, N)
        O = zeros(K, B)
        E = zeros(K, B)
        W = zeros(B + 1, d)
        Phi_Rk = zeros(B + 1, N)
        harmony_obj = new(Z_corr, Z_orig, Z_cos, phi, phi_moe, N, Pr_b, B, d, window_size,
                epsilon_cluster, epsilon_harmony, lamb_mat, sigma, sigma, block_size, K,
                max_iter_harmony, max_iter_kmeans, verbose, theta,
                objective_harmony, objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy,
                objective_kmeans_cross, kmeans_rounds, scale_dist, dist_mat, O, E, W, Phi_Rk,
                zeros(Float64, 0, 0), zeros(Float64, 0, 0))
        init_cluster!(harmony_obj, random_state)
        harmonize!(harmony_obj, max_iter_harmony, verbose)
        return harmony_obj
    end
end