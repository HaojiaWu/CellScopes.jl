
function get_object_group(obj_type::String)
    if obj_type == "Sequencing"
        objs = Union{VisiumObject, SlideseqObject, StereoSeqObject, VisiumHDObject}
    elseif obj_type == "Imaging"
        objs = Union{ImagingSpatialObject, CartanaObject, XeniumObject, MerfishObject, STARmapObject, seqFishObject, StereoSeqObject, CosMxObject}
    elseif obj_type == "Spatial"
        objs = Union{ImagingSpatialObject, CartanaObject, VisiumObject, VisiumHDObject, XeniumObject, MerfishObject, SlideseqObject, STARmapObject, seqFishObject, StereoSeqObject, CosMxObject}
    elseif obj_type == "Spatial2"
        objs = Union{ImagingSpatialObject, CartanaObject, VisiumObject, XeniumObject, MerfishObject, SlideseqObject, STARmapObject, seqFishObject, StereoSeqObject, CosMxObject}
    elseif obj_type == "All"
        objs = Union{scRNAObject, scATACObject, ImagingSpatialObject, CartanaObject, VisiumObject, VisiumHDObject, XeniumObject, MerfishObject, SlideseqObject, STARmapObject, seqFishObject, StereoSeqObject, CosMxObject, IntegratedObject}
    else
        error("Obj_type can only be Sequencing, Imaging, Spatial or All")
    end
    return objs
end

function Base.show(io::IO, sc_obj::scRNAObject)
    println(io, "scRNAObject in CellScopes.jl")
    println("Genes x Cells = ", length(sc_obj.rawCount.gene_name), " x ", length(sc_obj.rawCount.cell_name))
    println("Available data:")
    all_field = fieldnames(scRNAObject)
    reduce_field = fieldnames(ReductionObject)
    for i in all_field
        if isdefined(sc_obj, i)
            matchtype(getfield(sc_obj, i))
        end
    end

    if isdefined(sc_obj, :dimReduction)
        for i in reduce_field
            if isdefined(sc_obj.dimReduction, i)
                matchtype(getfield(sc_obj.dimReduction, i))
            end 
        end
    end
    println("All fields:")
    [println("- ", string(i)) for i in fieldnames(typeof(sc_obj))]
end

function Base.show(io::IO, sc_obj::scATACObject)
    println(io, "scATACObject in CellScopes.jl")
    println("Peaks x Cells = ", length(sc_obj.rawCount.gene_name), " x ", length(sc_obj.rawCount.cell_name))
    println("Available data:")
    all_field = fieldnames(scATACObject)
    reduce_field = fieldnames(ReductionObject)
    for i in all_field
        if isdefined(sc_obj, i)
            matchtype(getfield(sc_obj, i))
        end
    end

    if isdefined(sc_obj, :peakAnno)
        println("- Peak data")
    end

    if isdefined(sc_obj, :dimReduction)
        for i in reduce_field
            if isdefined(sc_obj.dimReduction, i)
                matchtype(getfield(sc_obj.dimReduction, i))
            end 
        end
    end

    println("All fields:")
    [println("- ", string(i)) for i in fieldnames(typeof(sc_obj))]
end

function Base.show(io::IO, count::AbstractCount)
    println(io, string(typeof(count)))
    println("Genes x Cells = ", length(count.gene_name), " x ", length(count.cell_name))
    println("All fields:")
    [println("- ", string(i)) for i in fieldnames(typeof(count))]
end

function Base.show(io::IO, obj_type::Union{AbstractDimReduction, ClusteringObject, VariableGeneObject, UndefinedObject, SpaImputeObj, VisiumImgObject, SpaCountObj, SpaMetaObj, SpaImageObj, Layer, Layers, Positions, AlterImages, AlterHDImgObject, HarmonyObject, AncillaryObject, AncillaryObjects})
    println(io, string(typeof(obj_type)))
    println("All fields:")
    [println("- ", string(i)) for i in fieldnames(typeof(obj_type))]
end

matchtype(field) = @match field begin
    raw::RawCountObject => println("- Raw count")
    norm::NormCountObject => println("- Normalized count")
    scale::ScaleCountObject => println("- Scaled count")
    meta::DataFrame => println("- Metadata")
    vargene::VariableGeneObject => println("- Variable genes")
    reduce::ReductionObject => println("- Dimensional reduction")
    pca::PCAObject => println("- PCA data")
    tsne::tSNEObject => println("- tSNE data")
    umap::UMAPObject => println("- UMAP data")
    cluster::ClusteringObject => println("- Clustering data")
    activity::GeneActivityObject => println("- Gene activity data")
    imputation::SpaImputeObj => println("- Imputed data")
    spMeta::SpaMetaObj => println("- Spatial metadata")
    spCoord::SpaCoordObj => println("- Spatial coordiantes")
    fragment::FragmentObject => println("- Fragment data")
    imgObject::SpaImageObj => println("- Image data")
    vsmImgObject::VisiumImgObject => println("- Visium image data")
    layer::Layer => println("- Layer")
    layers::Layers => println("- Layers data")
    positions::Positions => println("- Positions data")
    ancObject::AncillaryObject => println("- Ancillary object")
    ancObjects::AncillaryObjects => println("- Ancillary data")
    alterimg::AlterImages => println("- Processed image data")
    alterimgobj::AlterHDImgObject => println("- Alternative image data")
    polygonData::Polygons => println("- Polygon data")
    uns::UndefinedObject => println("- Undefined slot")
    emptyObj::Nothing => print("")
end

function Base.show(io::IO, sp_obj::get_object_group("All"))
    if isa(sp_obj, ImagingSpatialObject)
        println(io, "ImagingSpatialObject in CellScopes.jl")
    elseif isa(sp_obj, CartanaObject)
        println(io, "CartanaObject in CellScopes.jl")
    elseif isa(sp_obj, XeniumObject)
        println(io, "XeniumObject in CellScopes.jl")
    elseif isa(sp_obj, MerfishObject)
        println(io, "MerfishObject in CellScopes.jl")
    elseif isa(sp_obj, STARmapObject)
        println(io, "STARmapObject in CellScopes.jl")
    elseif isa(sp_obj, seqFishObject)
        println(io, "seqFishObject in CellScopes.jl")
    elseif isa(sp_obj, SlideseqObject)
        println(io, "SlideseqObject in CellScopes.jl")
    elseif isa(sp_obj, StereoSeqObject)
        println(io, "StereoSeqObject in CellScopes.jl")
    elseif isa(sp_obj, VisiumHDObject)
        println(io, "VisiumHDObject in CellScopes.jl")
    elseif isa(sp_obj, IntegratedObject)
        println(io, "IntegratedObject in CellScopes.jl")
    else
        println(io, "VisiumObject in CellScopes.jl")
    end
    println("Genes x Cells = ", size(sp_obj.rawCount.count_mtx)[1], " x ", size(sp_obj.rawCount.count_mtx)[2])
    println("Available data:")
    all_field = fieldnames(typeof(sp_obj))
    for i in all_field
        if isdefined(sp_obj, i)
            println("- ", string(i))
        end
    end
    println("All fields:")
    [println("- ", string(i)) for i in fieldnames(typeof(sp_obj))]
end
