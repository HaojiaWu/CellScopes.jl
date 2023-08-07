# AbstractSpaObj: type for spatial data

abstract type AbstractSpaObj <: AbstractCellScope end
abstract type AbstractImagingObj <: AbstractSpaObj end
abstract type AbstractSequencingObj <: AbstractSpaObj end

mutable struct SpaCountObj <: AbstractCount
    count_mtx::AbstractMatrix{<:Real}
    cell_name::Vector{String}
    gene_name::Vector{String}
    function SpaCountObj(count_mtx, cell_name, gene_name)
        if length(cell_name) !== size(count_mtx)[2]
            error("The total number of cells in the count matrix and the cell ID do not match!")
        end
        if length(gene_name) !== size(count_mtx)[1]
            error("The total number of genes in the count matrix and the gene name do not match!")
        end
        spa_count = new(count_mtx, cell_name, gene_name)
        return spa_count
    end
end

mutable struct SpaImputeObj <: AbstractImagingObj
    tgCount::Union{SpaCountObj, Nothing}
    spageCount::Union{SpaCountObj, Nothing}
    gimviCount::Union{SpaCountObj, Nothing}
    function SpaImputeObj(imp_type::String; 
        imp_data::Union{SpaCountObj, Nothing}=nothing
        )
        impute_obj = new(nothing, nothing, nothing)
        if imp_type === "tangram"
            impute_obj.tgCount = imp_data
        elseif imp_type === "SpaGE"
            impute_obj.spageCount = imp_data
        elseif imp_type === "gimVI"
            impute_obj.gimviCount = imp_data
        else
            error("imp_type can only be \"tangram\", \"SpaGE\" and \"gimVI\"")
        end
        return impute_obj
    end
end

mutable struct SpaMetaObj <: AbstractImagingObj
    cell::Union{DataFrame, Nothing}
    molecule::Union{DataFrame, Nothing}
    polygon::Union{DataFrame, Nothing}
    SpaMetaObj(cell, molecule, polygon) = new(cell, molecule, polygon)
end

mutable struct SpaCoordObj <: AbstractImagingObj
    cellCoord::Union{DataFrame, Nothing}
    molCoord::Union{DataFrame, Nothing}
    polygonCoord::Union{DataFrame, Nothing}
    otherCoord::Union{DataFrame, Nothing}
    SpaCoordObj(cellCoord, molCoord, polygonCoord, otherCoord) = new(cellCoord, molCoord, polygonCoord, otherCoord)
end

mutable struct ImagingSpatialObject <: AbstractImagingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    coordData::Union{SpaCoordObj, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    imageData::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}}
    polygonData::Union{Vector{Matrix{Float64}}, Nothing}
    function ImagingSpatialObject(molecule_data::DataFrame, cell_data::DataFrame, counts::RawCountObject; 
        prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, meta_data::Union{DataFrame, Nothing} = nothing,
        min_gene::Int64=0, min_cell::Int64=0, x_col::Union{String, Symbol} = "x", 
        y_col::Union{String, Symbol} = "y", cell_col::Union{String, Symbol} = "cell")
        if prefix !== nothing
            println("Adding prefix " * prefix * " to all cells...")
            counts.cell_name = prefix * "_" .* counts.cell_name
            molecule_data[!, cell_col] = prefix * "_" .* molecule_data[!, cell_col]
            cell_data[!, cell_col] = prefix * "_" .* cell_data[!, cell_col]
        end
        if postfix !== nothing
            println("Adding postfix " * postfix * " to all cells...")
            counts.cell_name = counts.cell_name .* "_" .* postfix
            molecule_data[!, cell_col] = molecule_data[!, cell_col] .* "_" .* postfix
            cell_data[!, cell_col] = cell_data[!, cell_col] .* "_" .* postfix
        end
        count_mat = counts.count_mtx
        gene_name = counts.gene_name
        cell_name = counts.cell_name
        count_mat, gene_name, cell_name = subset_matrix(count_mat, gene_name, cell_name, min_gene, min_cell)
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cell_name, nFeatures=nFeatures, nGenes = nGenes)
        end
        counts = RawCountObject(count_mat, cell_name, gene_name)
        cell_check = check_vec(cell_name, cell_data[!, cell_col])
        cell_data = cell_data[cell_check, :]
        mol_check = check_vec(cell_name, molecule_data[!, cell_col])
        molecule_data = molecule_data[mol_check, :]
        spObj=new(counts)
        meta = SpaMetaObj(cell_data, molecule_data, nothing)
        spObj.spmetaData = meta
        cell_coord = cell_data[!, [x_col, y_col]]
        mol_coord = molecule_data[!, [x_col, y_col]]
        coord = SpaCoordObj(cell_coord, mol_coord, nothing, nothing)
        spObj.coordData = coord
        spObj.metaData = meta_data
        spObj.polygonData = nothing
        return spObj
        println("ImagingSpatialObject was successfully created!")
    end
end

mutable struct CartanaObject <: AbstractImagingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    coordData::Union{SpaCoordObj, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    imageData::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}}
    polygonData::Union{Vector{Matrix{Float64}}, Nothing}
    function CartanaObject(molecule_data::DataFrame, cell_data::DataFrame, counts::RawCountObject; 
        prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, meta_data::Union{DataFrame, Nothing} = nothing,
        min_gene::Int64=0, min_cell::Int64=0, x_col::Union{String, Symbol} = "x", 
        y_col::Union{String, Symbol} = "y", cell_col::Union{String, Symbol} = "cell")
        if prefix !== nothing
            println("Adding prefix " * prefix * " to all cells...")
            counts.cell_name = prefix * "_" .* counts.cell_name
            molecule_data[!, cell_col] = prefix * "_" .* molecule_data[!, cell_col]
            cell_data[!, cell_col] = prefix * "_" .* cell_data[!, cell_col]
        end
        if postfix !== nothing
            println("Adding postfix " * postfix * " to all cells...")
            counts.cell_name = counts.cell_name .* "_" .* postfix
            molecule_data[!, cell_col] = molecule_data[!, cell_col] .* "_" .* postfix
            cell_data[!, cell_col] = cell_data[!, cell_col] .* "_" .* postfix
        end
        count_mat = counts.count_mtx
        gene_name = counts.gene_name
        cell_name = counts.cell_name
        count_mat, gene_name, cell_name = subset_matrix(count_mat, gene_name, cell_name, min_gene, min_cell)
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cell_name, nFeatures=nFeatures, nGenes = nGenes)
        end
        counts = RawCountObject(count_mat, cell_name, gene_name)
        cell_check = check_vec(cell_name, cell_data[!, cell_col])
        cell_data = cell_data[cell_check, :]
        mol_check = check_vec(cell_name, molecule_data[!, cell_col])
        molecule_data = molecule_data[mol_check, :]
        spObj=new(counts)
        meta = SpaMetaObj(cell_data, molecule_data, nothing)
        spObj.spmetaData = meta
        cell_coord = cell_data[!, [x_col, y_col]]
        mol_coord = molecule_data[!, [x_col, y_col]]
        coord = SpaCoordObj(cell_coord, mol_coord, nothing, nothing)
        spObj.coordData = coord
        spObj.metaData = meta_data
        spObj.polygonData = nothing
        return spObj
        println("CartanaObject was successfully created!")
    end
end

mutable struct XeniumObject <: AbstractImagingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    coordData::Union{SpaCoordObj, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    imageData::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}}
    polygonData::Array{Array{Float64, 2}, 1}

    function XeniumObject(molecule_data::DataFrame, cell_data::DataFrame, counts::RawCountObject, poly_data::Array{Array{Float64, 2}, 1}, umap_obj::UMAPObject; 
        prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, meta_data::Union{DataFrame, Nothing} = nothing,
        min_gene::Int64=0, min_cell::Int64=0, x_col::Union{String, Symbol} = "x", 
        y_col::Union{String, Symbol} = "y", cell_col::Union{String, Symbol} = "cell")
        if isa(prefix, String)
            println("Adding prefix " * prefix * " to all cells...")
            counts.cell_name = prefix * "_" .* counts.cell_name
            molecule_data[!, cell_col] = prefix * "_" .* molecule_data[!, cell_col]
            cell_data[!, cell_col] = prefix * "_" .* cell_data[!, cell_col]
        end
        if isa(postfix, String)
            println("Adding postfix " * postfix * " to all cells...")
            counts.cell_name = counts.cell_name .* "_" .* postfix
            molecule_data[!, cell_col] = molecule_data[!, cell_col] .* "_" .* postfix
            cell_data[!, cell_col] = cell_data[!, cell_col] .* "_" .* postfix
        end
        count_mat = counts.count_mtx
        gene_name = counts.gene_name
        cell_name = counts.cell_name
        if min_gene > 0 | min_cell > 0
            count_mat, gene_name, cell_name = subset_matrix(count_mat, gene_name, cell_name, min_gene, min_cell)
            cell_check = check_vec(cell_name, cell_data[!, cell_col])
            cell_data = cell_data[cell_check, :]
            mol_check = check_vec(cell_name, molecule_data[!, cell_col])
            molecule_data = molecule_data[mol_check, :]
        end
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cell_name, nFeatures=nFeatures, nGenes = nGenes)
        end
        counts = RawCountObject(count_mat, cell_name, gene_name)
        spObj = new(counts)
        polygon_df = DataFrame(polygon_number = 1:length(poly_data), mapped_cell = cell_data.cell, cluster=cell_data.cluster)
        meta = SpaMetaObj(cell_data, molecule_data, polygon_df)
        spObj.spmetaData = meta
        cell_coord = cell_data[!, [x_col, y_col]]
        mol_coord = molecule_data[!, [x_col, y_col]]
        coord = SpaCoordObj(cell_coord, mol_coord, nothing, nothing)
        spObj.coordData = coord
        spObj.metaData = meta_data
        spObj.polygonData = poly_data
        reduct_obj = ReductionObject(nothing, nothing, umap_obj)
        spObj.dimReduction = reduct_obj
        spObj = normalize_object(spObj)
        spObj.polynormCount = spObj.normCount
        return spObj
        println("XeniumObject was successfully created!")
    end
end

mutable struct VisiumImgObject <: AbstractImagingObj
    highresImage::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}, Nothing}
    lowresImage::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}, Nothing}
    fullresImage::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}, Nothing}
    detectedTissue::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}, Nothing}
    alignedImage::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}, Nothing}
    jsonParameters::Union{Dict{String, Any}, Nothing}
    VisiumImgObject(highresImage, lowresImage, fullresImage, detectedTissue, alignedImage, jsonParameters) = new(highresImage, lowresImage, fullresImage, detectedTissue, alignedImage, jsonParameters)
end

mutable struct VisiumObject <: AbstractSequencingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{DataFrame, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    imageData::Union{VisiumImgObject, Nothing}
    function VisiumObject(raw_count::RawCountObject; 
            meta_data::Union{DataFrame, Nothing} = nothing,
            sp_meta::Union{DataFrame, Nothing} = nothing,
            min_gene::Int64 = 0,
            min_cell::Int64 = 0,
            prefix::Union{String, Nothing} = nothing,
            postfix::Union{String, Nothing} = nothing)
        count_mat = raw_count.count_mtx
        genes = raw_count.gene_name
        cells = raw_count.cell_name
        count_mat, genes, cells = subset_matrix(count_mat, genes, cells, min_gene, min_cell)
        if prefix !== nothing
            println("Adding prefix " * prefix * " to all cells...")
            cellnames = prefix * "_" .* raw_count.cell_name
        end
        if postfix !== nothing
            println("Adding postfix " * postfix * " to all cells...")
            cellnames = raw_count.cell_name .* "_" .* postfix
        end
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cells, nFeatures=nFeatures, nGenes = nGenes)
        end
        count_obj = RawCountObject(count_mat, cells, genes)
        visium_obj = new(count_obj)
        if sp_meta !== nothing
            sp_meta = filter(:cell => ∈(Set(meta_data.Cell_id)), sp_meta)
            visium_obj.spmetaData = sp_meta
        end
        visium_obj.metaData = meta_data
        return visium_obj
        println("VisiumObject was successfully created!")
    end
end

function add_impdata(impute_obj::SpaImputeObj, imp_type::String, imp_data::SpaCountObj)
        if imp_type === "tangram"
            impute_obj.tgCount = imp_data
        elseif imp_type === "SpaGE"
            impute_obj.spageCount = imp_data
        elseif imp_type === "gimVI"
            impute_obj.gimviCount = imp_data
        else
            error("imp_type can only be \"tangram\", \"SpaGE\" and \"gimVI\"")
        end
        return impute_obj
end

mutable struct MerfishObject <: AbstractImagingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    coordData::Union{SpaCoordObj, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    imageData::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}}
    polygonData::Array{Array{Float64, 2}, 1}

    function MerfishObject(molecule_data::DataFrame, cell_data::DataFrame, counts::RawCountObject, poly_data::Array{Array{Float64, 2}, 1}; 
        prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, meta_data::Union{DataFrame, Nothing} = nothing,
        min_gene::Int64=0, min_cell::Int64=0, x_col::Union{String, Symbol} = "x", 
        y_col::Union{String, Symbol} = "y", cell_col::Union{String, Symbol} = "cell")
        if prefix !== nothing
            println("Adding prefix " * prefix * " to all cells...")
            counts.cell_name = prefix * "_" .* counts.cell_name
            cell_data[!, cell_col] = prefix * "_" .* cell_data[!, cell_col]
        end
        if postfix !== nothing
            println("Adding postfix " * postfix * " to all cells...")
            counts.cell_name = counts.cell_name .* "_" .* postfix
            cell_data[!, cell_col] = cell_data[!, cell_col] .* "_" .* postfix
        end
        count_mat = raw_count.count_mtx
        gene_name = raw_count.gene_name
        cell_name = raw_count.cell_name
        count_mat, gene_name, cell_name = subset_matrix(count_mat, gene_name, cell_name, min_gene, min_cell)
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cell_name, nFeatures=nFeatures, nGenes = nGenes)
        end
        counts = RawCountObject(count_mat, cell_name, gene_name)
        cell_check = check_vec(cell_name, cell_data[!, cell_col])
        cell_data = cell_data[cell_check, :]
        spObj = new(counts)
        polygon_df = DataFrame(polygon_number = 1:length(poly_data), mapped_cell = cell_data.cell)
        meta = SpaMetaObj(cell_data, molecule_data, polygon_df)
        spObj.spmetaData = meta
        cell_coord = cell_data[!, [x_col, y_col]]
        mol_coord = molecule_data[!, [x_col, y_col]]
        coord = SpaCoordObj(cell_coord, mol_coord, nothing, nothing)
        spObj.coordData = coord
        spObj.metaData = meta_data
        spObj.polygonData = poly_data
        spObj = normalize_object(spObj)
        spObj.polynormCount = spObj.normCount
        replace!(spObj.polynormCount.count_mtx, NaN=>0)
        return spObj
        println("MerfishObject was successfully created!")
    end
end

mutable struct SlideseqObject <: AbstractSequencingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{DataFrame, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    function SlideseqObject(raw_count::RawCountObject; 
            meta_data::Union{DataFrame, Nothing} = nothing,
            sp_meta::Union{DataFrame, Nothing} = nothing,
            min_gene::Int64 = 0,
            min_cell::Int64 = 0,
            prefix::Union{String, Nothing} = nothing,
            postfix::Union{String, Nothing} = nothing)
        count_mat = raw_count.count_mtx
        genes = raw_count.gene_name
        cells = raw_count.cell_name
        count_mat, genes, cells = subset_matrix(count_mat, genes, cells, min_gene, min_cell)
        if prefix !== nothing
            println("Adding prefix " * prefix * " to all cells...")
            cellnames = prefix * "_" .* raw_count.cell_name
        end
        if postfix !== nothing
            println("Adding postfix " * postfix * " to all cells...")
            cellnames = raw_count.cell_name .* "_" .* postfix
        end
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cells, nFeatures=nFeatures, nGenes = nGenes)
        end
        count_obj = RawCountObject(count_mat, cells, genes)
        slideseq_obj = new(count_obj)
        if sp_meta !== nothing
            sp_meta = filter(:cell => ∈(Set(meta_data.Cell_id)), sp_meta)
            slideseq_obj.spmetaData = sp_meta
        end
        slideseq_obj.metaData = meta_data
        return slideseq_obj
        println("SlideseqObject was successfully created!")
    end
end

mutable struct STARmapObject <: AbstractImagingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    coordData::Union{SpaCoordObj, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    polygonData::Union{Vector{Matrix{Float64}}, Nothing}
    function STARmapObject(molecule_data::DataFrame, cell_data::DataFrame, counts::RawCountObject; 
        prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, meta_data::Union{DataFrame, Nothing} = nothing,
        min_gene::Int64=0, min_cell::Int64=0, x_col::Union{String, Symbol} = "x", 
        y_col::Union{String, Symbol} = "y", cell_col::Union{String, Symbol} = "cell")
        if prefix !== nothing
            println("Adding prefix " * prefix * " to all cells...")
            counts.cell_name = prefix * "_" .* counts.cell_name
            molecule_data[!, cell_col] = prefix * "_" .* molecule_data[!, cell_col]
            cell_data[!, cell_col] = prefix * "_" .* cell_data[!, cell_col]
        end
        if postfix !== nothing
            println("Adding postfix " * postfix * " to all cells...")
            counts.cell_name = counts.cell_name .* "_" .* postfix
            molecule_data[!, cell_col] = molecule_data[!, cell_col] .* "_" .* postfix
            cell_data[!, cell_col] = cell_data[!, cell_col] .* "_" .* postfix
        end
        count_mat = raw_count.count_mtx
        gene_name = raw_count.gene_name
        cell_name = raw_count.cell_name
        count_mat, gene_name, cell_name = subset_matrix(count_mat, gene_name, cell_name, min_gene, min_cell)
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cell_name, nFeatures=nFeatures, nGenes = nGenes)
        end
        counts = RawCountObject(count_mat, cell_name, gene_name)
        cell_check = check_vec(cell_name, cell_data[!, cell_col])
        cell_data = cell_data[cell_check, :]
        mol_check = check_vec(cell_name, molecule_data[!, cell_col])
        molecule_data = molecule_data[mol_check, :]
        spObj=new(counts)
        meta = SpaMetaObj(cell_data, molecule_data, nothing)
        spObj.spmetaData = meta
        cell_coord = cell_data[!, [x_col, y_col]]
        mol_coord = molecule_data[!, [x_col, y_col]]
        coord = SpaCoordObj(cell_coord, mol_coord, nothing, nothing)
        spObj.coordData = coord
        spObj.metaData = meta_data
        spObj.polygonData = nothing
        return spObj
        println("STARmapObject was successfully created!")
    end
end

mutable struct seqFishObject <: AbstractImagingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{DataFrame, Nothing}
    spmetaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    coordData::Union{SpaCoordObj, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    polygonData::Union{Vector{Matrix{Float64}}, Nothing}
    function seqFishObject(molecule_data::DataFrame, cell_data::DataFrame, counts::RawCountObject; 
        prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, meta_data::Union{DataFrame, Nothing} = nothing,
        min_gene::Int64=0, min_cell::Int64=0, x_col::Union{String, Symbol} = "x", 
        y_col::Union{String, Symbol} = "y", cell_col::Union{String, Symbol} = "cell")
        if prefix !== nothing
            println("Adding prefix " * prefix * " to all cells...")
            counts.cell_name = prefix * "_" .* counts.cell_name
            molecule_data[!, cell_col] = prefix * "_" .* molecule_data[!, cell_col]
            cell_data[!, cell_col] = prefix * "_" .* cell_data[!, cell_col]
        end
        if postfix !== nothing
            println("Adding postfix " * postfix * " to all cells...")
            counts.cell_name = counts.cell_name .* "_" .* postfix
            molecule_data[!, cell_col] = molecule_data[!, cell_col] .* "_" .* postfix
            cell_data[!, cell_col] = cell_data[!, cell_col] .* "_" .* postfix
        end
        count_mat = raw_count.count_mtx
        gene_name = raw_count.gene_name
        cell_name = raw_count.cell_name
        count_mat, gene_name, cell_name = subset_matrix(count_mat, gene_name, cell_name, min_gene, min_cell)
        if isa(meta_data, Nothing)
            nFeatures = vec(colSum(count_mat))
            nGenes = vec(sum(x->x>0, count_mat, dims=1))
            meta_data = DataFrame(Cell_id = cell_name, nFeatures=nFeatures, nGenes = nGenes)
        end
        counts = RawCountObject(count_mat, cell_name, gene_name)
        cell_check = check_vec(cell_name, cell_data[!, cell_col])
        cell_data = cell_data[cell_check, :]
        mol_check = check_vec(cell_name, molecule_data[!, cell_col])
        molecule_data = molecule_data[mol_check, :]
        spObj=new(counts)
        meta = SpaMetaObj(cell_data, molecule_data, nothing)
        spObj.spmetaData = meta
        cell_coord = cell_data[!, [x_col, y_col]]
        mol_coord = molecule_data[!, [x_col, y_col]]
        coord = SpaCoordObj(cell_coord, mol_coord, nothing, nothing)
        spObj.coordData = coord
        spObj.metaData = meta_data
        spObj.polygonData = nothing
        return spObj
        println("seqFishObject was successfully created!")
    end
end
