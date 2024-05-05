# Types for high resolution full transcriptome spatial data

abstract type AbstractSpaFullObj <: AbstractCellScope end
abstract type AbstractLayers <: AbstractSpaFullObj end
abstract type AbstractPositions <: AbstractSpaFullObj end
abstract type AbstractHDImages <: AbstractSpaFullObj end

mutable struct Layer <: AbstractLayers
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{DataFrame, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polygonData::Union{Array{Array{Float64, 2}, 1}, Nothing}
    jsonParameters::Union{Dict{String, Any}, Nothing}
    function Layer(raw_count; 
        min_gene::Int64 = 0,
        min_cell::Int64 = 0,
        meta_data::Union{DataFrame, Nothing} = nothing,
        prefix::Union{String, Nothing} = nothing,
        postfix::Union{String, Nothing} = nothing)
        count_mat = raw_count.count_mtx
        genes = raw_count.gene_name
        cells = raw_count.cell_name
        count_mat, genes, cells = subset_matrix(count_mat, genes, cells, min_gene, min_cell)
        if isa(prefix, String)
            println("Adding prefix " * prefix * " to all cells...")
            cellnames = prefix * "_" .* raw_count.cell_name
        end
        if isa(postfix, String)
            println("Adding postfix " * postfix * " to all cells...")
            cellnames = raw_count.cell_name .* "_" .* postfix
        end
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cells, nFeatures=nFeatures, nGenes = nGenes)
        end
        count_obj = RawCountObject(count_mat, cells, genes)
        layer_obj = new(count_obj, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
        layer_obj.metaData = meta_data
        return layer_obj
    end
end

mutable struct Layers <: AbstractLayers
    layers::Dict{String, Layer}
    Layers() = new(Dict{String, Layer}())
end

mutable struct Positions <: AbstractPositions
    Positions() = new(Dict{String, DataFrame}())
end

mutable struct AlterImages <: AbstractHDImages
    AlterImages() = new(Dict{String, Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}, Nothing}}())
end

mutable struct AlterHDImgObject <: AbstractHDImages
    imgData::Union{AlterImages, Nothing}
    posData::Union{Positions, Nothing}
    AlterHDImgObject = new(imageData, posData)
end

mutable struct VisiumHDObject <: AbstractSpaFullObj
    layerData::Union{Layers, Nothing}
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{DataFrame, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    imageData::Union{VisiumImgObject, Nothing}
    alterImgData::Union{AlterHDImgObject, Nothing}
    polygonData::Union{Array{Array{Float64, 2}, 1}, Nothing}
    defaultData::Union{String, Nothing}
    function VisiumHDObject(layer_data::Layers; 
            default_bin::String="8_um",
            meta_data::Union{DataFrame, Nothing} = nothing,
            sp_meta::Union{DataFrame, Nothing} = nothing,
            min_gene::Int64 = 0,
            min_cell::Int64 = 0,
            prefix::Union{String, Nothing} = nothing,
            postfix::Union{String, Nothing} = nothing)
        hd_obj = new(layer_data, nothing, nothing, nothing,nothing,nothing,nothing,nothing,nothing,nothing, nothing, nothing,default_bin)
        raw_count = hd_obj.layerData.layers[default_bin].rawCount
        meta_data = hd_obj.layerData.layers[default_bin].metaData
        sp_meta = hd_obj.layerData.layers[default_bin].spmetaData
        count_mat = raw_count.count_mtx
        genes = raw_count.gene_name
        cells = raw_count.cell_name
        count_mat, genes, cells = subset_matrix(count_mat, genes, cells, min_gene, min_cell)
        if prefix == String
            println("Adding prefix " * prefix * " to all cells...")
            cellnames = prefix * "_" .* raw_count.cell_name
        end
        if postfix == String
            println("Adding postfix " * postfix * " to all cells...")
            cellnames = raw_count.cell_name .* "_" .* postfix
        end
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cells, nFeatures=nFeatures, nGenes = nGenes)
            hd_obj.layerData.layers[default_bin].metaData = meta_data
        end
        count_obj = RawCountObject(count_mat, cells, genes)
        hd_obj.layerData.layers[default_bin].rawCount = count_obj
        if sp_meta !== nothing
            sp_meta = filter(:barcode => ∈(Set(meta_data.Cell_id)), sp_meta)
            hd_obj.layerData.layers[default_bin].spmetaData = sp_meta
        end
        image_obj = VisiumImgObject(nothing, nothing, nothing, nothing, nothing, nothing)
        hd_obj.imageData = image_obj
        update_data!(hd_obj)
        return hd_obj
        println("VisiumHDObject was successfully created!")
    end
end

function update_data!(obj::VisiumHDObject)
    layer = get(obj.layerData.layers, obj.defaultData, nothing)
    if layer !== nothing
        obj.rawCount = layer.rawCount
        obj.normCount = layer.normCount
        obj.scaleCount = layer.scaleCount
        obj.metaData = layer.metaData
        obj.spmetaData = layer.spmetaData
        obj.varGene = layer.varGene
        obj.dimReduction = layer.dimReduction
        obj.clustData = layer.clustData
        obj.polygonData = layer.polygonData
        obj.imageData.jsonParameters  = layer.jsonParameters
    end
end

function Base.setproperty!(obj::VisiumHDObject, sym::Symbol, value)
    setfield!(obj, sym, value)
    if sym == :defaultData
        update_data!(obj)
    end
end

