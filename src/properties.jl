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

function Base.show(io::IO, obj_type::Union{AbstractDimReduction, ClusteringObject, VariableGeneObject, UndefinedObject, SpaImputeObj, VisiumImgObject, SpaCountObj,SpaMetaObj})
    println(io, string(typeof(obj_type)))
    println("All fields:")
    [println("- ", string(i)) for i in fieldnames(typeof(obj_type))]
end

matchtype(field) = @match field begin
    raw::RawCountObject => println("- Raw count")
    norm::NormCountObject => println("- Normalized count")
    scale::ScaleCountObject => println("- Scaled count")
    vargene::VariableGeneObject => println("- Variable genes")
    pca::PCAObject => println("- PCA data")
    tsne::tSNEObject => println("- tSNE data")
    umap::UMAPObject => println("- UMAP data")
    cluster::ClusteringObject => println("- Clustering data")
    activity::GeneActivityObject => println("- Gene activity data")
    imputation::SpaImputeObj => println("- Imputed data")
    spMeta::SpaMetaObj => println("- Spatial metadata")
    spCoord::SpaCoordObj => println("- Spatial coordiantes")
    fragment::FragmentObject => println("- Fragment data")
    uns::UndefinedObject => println("- Undefined slot")
end

function Base.show(io::IO, sp_obj::Union{CartanaObject, VisiumObject, XeniumObject, MerfishObject, SlideseqObject, starMapObject, seqFishObject})
    if isa(sp_obj, CartanaObject)
        println(io, "CartanaObject in CellScopes.jl")
    elseif isa(sp_obj, XeniumObject)
        println(io, "XeniumObject in CellScopes.jl")
    elseif isa(sp_obj, MerfishObject)
        println(io, "MerfishObject in CellScopes.jl")
    elseif isa(sp_obj, starMapObject)
        println(io, "starMapObject in CellScopes.jl")
    elseif isa(sp_obj, seqFishObject)
        println(io, "seqFishObject in CellScopes.jl")
    elseif isa(sp_obj, SlideseqObject)
        println(io, "SlideseqObject in CellScopes.jl")
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
