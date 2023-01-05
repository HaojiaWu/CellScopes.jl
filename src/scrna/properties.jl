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
    println("Available fields:")
    [println("- ", string(i)) for i in fieldnames(typeof(sc_obj))]
end

function Base.show(io::IO, count::AbstractCount)
    println(io, string(typeof(count)))
    println("Genes x Cells = ", length(count.gene_name), " x ", length(count.cell_name))
    println("Available fields:")
    [println("- ", string(i)) for i in fieldnames(typeof(count))]
end

function Base.show(io::IO, obj_type::Union{AbstractDimReduction, ClusteringObject, VariableGeneObject, UndefinedObject})
    println(io, string(typeof(obj_type)))
    println("Available fields:")
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
    uns::UndefinedObject => println("- Undefined slot")
end
