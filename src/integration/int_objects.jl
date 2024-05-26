
abstract type AbstractAncillaryObject <: AbstractCellScope end

mutable struct AncillaryObject <: AbstractAncillaryObject
    spmetaData::Union{SpaMetaObj, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    imageData::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}, Nothing}
    polygonData::Union{Array{Array{Float64, 2}, 1}, Nothing}
    function AncillaryObject(sp_obj)
        ancillary_obj = fill_nothing(AncillaryObject)
        for field in (:spmetaData, :polyCount, :polynormCount, :imputeData, :polygonData)
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

    # prepare ancillary objects
    seq_types = (scRNAObject, SlideseqObject, scATACObject, VisiumObject)
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
    new_meta = []
    for i in 1:length(meta_data)
        meta1 = meta_data[i]
        meta1.dataset .= sample_names[i]
        meta1 = meta1[!, [:Cell_id, :dataset]]
        new_meta[i] = meta1
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