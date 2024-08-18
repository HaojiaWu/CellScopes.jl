
abstract type AbstractAncillaryObject <: AbstractCellScope end
abstract type AbstractHarmony <: AbstractSingleCell end

mutable struct AncillaryObject <: AbstractAncillaryObject
    layerData::Union{Layers, Nothing}
    spmetaData::Union{SpaMetaObj, DataFrame, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    imageData::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}, VisiumImgObject, Nothing}
    polygonData::Union{Array{Array{Float64, 2}, 1}, Nothing}
    alterImgData::Union{AlterHDImgObject, Nothing}
    defaultData::Union{String, Nothing}
    function AncillaryObject(sp_obj)
        ancillary_obj = fill_nothing(AncillaryObject)
        for field in (:layerData, :spmetaData, :polyCount, :polynormCount, :imputeData, :polygonData, :alterImgData, :defaultData)
            if isdefined(sp_obj, field)
                setfield!(ancillary_obj, field, getfield(sp_obj, field))
            end
        end
        return ancillary_obj
    end
end

mutable struct AncillaryObjects <: AbstractAncillaryObject
    ancillaryObjs::Dict{String, AncillaryObject}
    AncillaryObjects() = new(Dict{String, AncillaryObject}())
end

mutable struct IntegratedObject <: AbstractCellScope
    ancillaryObjs::Union{AncillaryObjects, Nothing}
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    dataType::Union{String, Nothing}
end

function CreateIntegratedObject(obj_list;
    sample_names::Union{Vector{String}, Nothing}=nothing, 
    min_gene::Union{Int64, Float64, Nothing} = nothing,
    min_cell::Union{Int64, Float64, Nothing} = nothing
)
    int_obj = fill_nothing(IntegratedObject)
    if isa(sample_names, Nothing)
        sample_names = string.("dataset_", 1:length(obj_list))
    end
    int_obj.dataType = get_type(obj_list)
    # prepare ancillary objects
    seq_types = [scRNAObject, SlideseqObject, scATACObject, VisiumObject]
    if any(x -> typeof(x) in seq_types, obj_list)
        ancillary_objs = nothing
    else
        ancillary_objs = AncillaryObjects()
        for i in 1:length(obj_list)
            ancillary_obj = AncillaryObject(obj_list[i])
            ancillary_objs.ancillaryObjs[sample_names[i]] = ancillary_obj
        end
    end
    int_obj.ancillaryObjs = ancillary_objs

    # merge raw counts
    ct_objs = [getfield(ct_obj, :rawCount) for ct_obj in obj_list]
    counts_all = [getfield(raw_count, :count_mtx) for raw_count in ct_objs]
    genes_all = [getfield(genes, :gene_name) for genes in ct_objs]
    cells_all = [getfield(cells, :cell_name) for cells in ct_objs]
    cells_all = vcat(cells_all...)
    merged_mtx, genes_all2 = merge_matrices(counts_all, genes_all)
    if !isa(min_gene, Nothing) && !isa(min_cell, Nothing)
        count_mat, genes, cells = subset_matrix(merged_mtx, genes_all2, cells_all, min_gene, min_cell)
        merged_ct = RawCountObject(count_mat, cells, genes)
    else
        merged_ct = RawCountObject(merged_mtx, cells_all, genes_all2)
    end
    int_obj.rawCount = merged_ct

    # merge metadata
    meta_data = [getfield(meta, :metaData) for meta in obj_list]
    new_meta = Vector{DataFrame}()
    for i in 1:length(meta_data)
        meta1 = meta_data[i]
        meta1.dataset .= sample_names[i]
        meta1 = meta1[!, [:Cell_id, :dataset]]
        push!(new_meta, meta1)
    end
    combined_meta = vcat(new_meta...)
    nFeatures = vec(colSum(merged_ct.count_mtx))
    nGenes = vec(sum(x->x>0, merged_ct.count_mtx, dims=1))
    combined_meta.nFeatures = nFeatures
    combined_meta.nGenes = nGenes
    int_obj.metaData = combined_meta
    return int_obj
    println("IntegratedObject was successfully created!")
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

mutable struct PairedSpObject <: AbstractCellScope
    vsObj::Union{VisiumHDObject, Nothing}
    xnObj::Union{XeniumObject, Nothing}
    vsMat::Union{Matrix{Float64}, Nothing} 
    xnMat::Union{Matrix{Float64}, Nothing} 
    PairedSpObject(vsObj, xnObj, vsMat, xnMat)=new(vsObj, xnObj, vsMat, xnMat)
end

mutable struct PairedObject <: AbstractCellScope
    pairedData::Union{PairedSpObject, Nothing}
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polygonData::Array{Array{Float64, 2}, 1}
    function PairedObject(paired_obj::PairedSpObject, counts::RawCountObject;        
        prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, meta_data::Union{DataFrame, Nothing} = nothing,
        min_gene::Int64=0, min_cell::Int64=0)
        cell_name = counts.cell_name
        if isa(prefix, String)
            println("Adding prefix " * prefix * " to all cells...")
            counts.cell_name = prefix * "_" .* counts.cell_name
        end
        if isa(postfix, String)
            println("Adding postfix " * postfix * " to all cells...")
            counts.cell_name = counts.cell_name .* "_" .* postfix
        end
        count_mat = counts.count_mtx
        gene_name = counts.gene_name
        cell_name = counts.cell_name
        if min_gene > 0 || min_cell > 0
            count_mat, gene_name, cell_name = subset_matrix(count_mat, gene_name, cell_name, min_gene, min_cell)
        end
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cell_name, nFeatures=nFeatures, nGenes = nGenes)
        end
        counts = RawCountObject(count_mat, cell_name, gene_name)
        pairedObj = new(paired_obj, counts)
        pairedObj.metaData = meta_data
        return pairedObj
        println("PairedObject was successfully created!")
    end
end