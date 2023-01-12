# AbstractSpaObj: type for spatial data

abstract type AbstractSpaObj <: AbstractCellScope end
abstract type AbstractImagingObj <: AbstractSpaObj end
abstract type AbstractSequencingObj <: AbstractSpaObj end

mutable struct SpaCountObj <: AbstractSpaObj
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

mutable struct CartanaObject <: AbstractImagingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    polyCount::Union{RawCountObject, Nothing}
    polynormCount::Union{NormCountObject, Nothing}
    coordData::Union{SpaCoordObj, Nothing}
    imputeData::Union{SpaImputeObj, Nothing}
    imageData::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}}
    polygonData::Array{Array{Float64, 2}, 1}
    function CartanaObject(molecule_data::DataFrame, cell_data::DataFrame, counts::RawCountObject; 
        prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, 
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
        gene_kept = (vec ∘ collect)(rowSum(count_mat).>= min_gene)
        gene_name = gene_name[gene_kept]
        cell_kept = (vec ∘ collect)(colSum(count_mat) .>= min_cell)
        cell_name = cell_name[cell_kept]
        count_mat = count_mat[gene_kept, cell_kept]
        counts = RawCountObject(count_mat, cell_name, gene_name)
        cell_check = Folds.collect(x in cell_name for x in cell_data[!, cell_col])
        cell_data = cell_data[cell_check, :]
        cell_check = Folds.collect(x in cell_name for x in molecule_data[!, cell_col])
        molecule_data = molecule_data[cell_check, :]
        spObj=new(counts)
        meta = SpaMetaObj(cell_data, nothing, nothing)
        meta.molecule = molecule_data
        spObj.metaData = meta
        cell_coord = cell_data[!, [x_col, y_col]]
        mol_coord = molecule_data[!, [x_col, y_col]]
        coord = SpaCoordObj(cell_coord, mol_coord, nothing, nothing)
        spObj.dimData = coord
        println("CartanaObject was successfully created!")
        return spObj
    end
end

mutable struct VisiumObject <: AbstractSequencingObj
    rawCount::Union{RawCountObject, Nothing}
    normCount::Union{NormCountObject, Nothing}
    scaleCount::Union{ScaleCountObject, Nothing}
    metaData::Union{SpaMetaObj, Nothing}
    varGene::Union{VariableGeneObject, Nothing}
    dimReduction::Union{ReductionObject, Nothing}
    clustData::Union{ClusteringObject, Nothing}
    coordData::Union{SpaCoordObj, Nothing}
    imageData::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}}
    function VisiumObject(raw_count::RawCountObject; 
            meta_data::Union{DataFrame, Nothing} = nothing,
            min_gene::Int64 = 0,
            min_cell::Int64 = 0,
            prefix::Union{String, Nothing} = nothing,
            postfix::Union{String, Nothing} = nothing)
        count_mat = raw_count.count_mtx
        genes = raw_count.gene_name
        cells = raw_count.cell_name
        gene_kept = (vec ∘ collect)(rowSum(count_mat).> min_gene)
        genes = genes[gene_kept]
        cell_kept = (vec ∘ collect)(colSum(count_mat) .> min_cell)
        cells = cells[cell_kept]
        count_mat = count_mat[gene_kept, cell_kept]
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
        visium_obj = new(count_obj)
        meta = SpaMetaObj(meta_data)
        visium_obj.metaData = meta
        println("VisiumObject was successfully created!")
        return visium_obj
    end
end

function add_impdata(impute_obj::SpaImputeObj, imp_type::String, imp_data::DataFrame)
        if imp_type === "tangram"
            impute_obj.tg_data = imp_data
        elseif imp_type === "SpaGE"
            impute_obj.spage_data = imp_data
        elseif imp_type === "gimVI"
            impute_obj.gimvi_data = imp_data
        else
            error("imp_type can only be \"tangram\", \"SpaGE\" and \"gimVI\"")
        end
        return impute_obj
end

function normalizeData(sp::AbstractSpaObj)
    if isa(sp.norm_counts, DataFrame)
        println("Your data has been normalized. No need to normalize again.")
    else
        orig_count=deepcopy(sp.counts)
        for celln in DataFrames.names(orig_count)[2:end]
            if celln !==String
                celln = string(celln)
            end
            orig_count[!, celln] = orig_count[!, celln] ./ sum(orig_count[!, celln])
        end
        new_df=Matrix(orig_count[!,2:end])
        dt = StatsBase.fit(UnitRangeTransform, new_df, dims=2)
        new_df=StatsBase.transform(dt, new_df)
        new_df=DataFrame(new_df,:auto)
        replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
        new_df = map(replace_nan, eachcol(new_df))
        new_df = DataFrame(new_df, :auto)
        DataFrames.rename!(new_df, DataFrames.names(orig_count)[2:end])
        new_df[!,:gene]=orig_count[!,:gene]
        total_cell=length(names(new_df))
        new_df=new_df[!,[total_cell; collect(1:total_cell-1)]]
        sp.norm_counts = new_df
    end
    return sp
end

function normalizeData(sp::DataFrame)
    orig_count=deepcopy(sp)
    for celln in DataFrames.names(orig_count)[2:end]
        if celln !==String
            celln = string(celln)
        end
        orig_count[!, celln] = orig_count[!, celln] ./ sum(orig_count[!, celln])
    end
    new_df=Matrix(orig_count[!,2:end])
    dt = StatsBase.fit(UnitRangeTransform, new_df, dims=2)
    new_df=StatsBase.transform(dt, new_df)
    new_df=DataFrame(new_df,:auto)
    replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
    new_df = map(replace_nan, eachcol(new_df))
    new_df = DataFrame(new_df, :auto)
    DataFrames.rename!(new_df, DataFrames.names(orig_count)[2:end])
    new_df[!,:gene]=orig_count[!,:gene]
    total_cell=length(names(new_df))
    new_df=new_df[!,[total_cell; collect(1:total_cell-1)]]
    return new_df
end

function polygons_cell_mapping(sp::AbstractSpaObj)
    polygons=sp.polygons
    center_df=DataFrame()
    cell=Int[]
    X=Float64[]
    Y=Float64[]
    for (i,p) in enumerate(polygons)
        x=mean(p[:,1])
        y=mean(p[:,2])
        push!(cell,i)
        push!(X,x)
        push!(Y,y)
    end
    center_df.cell=cell
    center_df.X=X
    center_df.Y=Y
    cell_coord=sp.cell_coord
    poly_mtx=Matrix(center_df[:,2:3])
    ref_mtx=Matrix(cell_coord[:,2:3])
    ref_coord=[]
    for i in eachrow(ref_mtx)
        ref_coord=push!(ref_coord,i)
    end
    query_coord=[]
    for i in eachrow(poly_mtx)
        query_coord=push!(query_coord,i)
    end
    n=length(query_coord)
    query_dist=[]
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    Threads.@threads for j in 1:length(query_coord)
        ref_dist=[]
        for i in 1:length(ref_coord)
            ref_dist=push!(ref_dist, Distances.euclidean(query_coord[j],ref_coord[i]))
        end
        query_dist=push!(query_dist,findmin(ref_dist)[2])
        next!(p)
    end
    println("Polygons have been mapped to the cells!")
    my_annot=DataFrame(polygon_number=collect(1:length(query_dist)), mapped_cell=query_dist);
    from=sp.cells.cell
    to=sp.cells.cluster
    anno=mapvalues(my_annot, :mapped_cell, :cluster,from, to)
    sp.poly_anno=anno
    return sp
end

function generate_polygon_counts(sp::AbstractSpaObj)
    coord_molecules=sp.molecules
    if isa(sp.poly_anno, Nothing)
        error("Please run polygons_cell_mapping first!")
    end
    anno=sp.poly_anno
    anno=rename!(anno, :polygon_number => :number)
    anno=rename!(anno, :mapped_cell => :cell)
    anno.cell=string.(anno.cell)
    coord_molecules.cell=string.(coord_molecules.cell)
    join_df=innerjoin(coord_molecules[!,[:gene,:cell]],anno[!,[:number,:cell]], on=:cell)
    join_df.cell=parse.(Int64, join_df.cell)
    join_df.gene=string.(join_df.gene)
    gdf = groupby(join_df, :gene);
    new_df=DataFrame()
    for (i, df) in enumerate(gdf)
        map_dict=countmap(df.number)
        anno1=DataFrame(cell=collect(keys(map_dict)),count=collect(values(map_dict)))
        anno1.gene.=gdf[i].gene[1]
        new_df=[new_df;anno1]
    end
    final_df=unstack(new_df, :cell, :gene, :count)
    final_df .= ifelse.(isequal.(final_df, missing), 0, final_df)
    final_df=mapcols(ByRow(Int64), final_df)
    sp.poly_counts=final_df
    println("Poly_counts was added to SpaData!")
    norm_df=mapcols(ByRow(Float64), final_df[!,2:end])
    norm_df=Matrix(norm_df)
    norm_df=norm_df'
    dt = StatsBase.fit(UnitRangeTransform, norm_df, dims=2)
    norm_df=StatsBase.transform(dt, norm_df)
    norm_df=norm_df'
    norm_df=DataFrame(norm_df,:auto)
    DataFrames.rename!(norm_df, DataFrames.names(final_df)[2:end])
    norm_df.cell=final_df.cell
    sort!(norm_df, :cell)
    sp.poly_norm=norm_df
    println("Poly_counts was normalized!")
    return sp
end

function subset_SpaObj(sp::AbstractSpaObj, cell_col::Union{String, Symbol}, subset_names::Union{Vector{String}, Vector{Int64}})
    spObj=deepcopy(sp)
    barcodes = deepcopy(subset_names)
    barcodes2 = [["gene"]; barcodes]
    spObj.norm_counts=spObj.norm_counts[!, barcodes2]
    spObj.counts=spObj.counts[!, barcodes2]
    spObj.cells=filter(cell_col => x -> x in subset_names, spObj.cells)
    if isa(sp, SpaObj)
        spObj.cell_names=subset_names
        spObj.molecules=filter(cell_col => x -> x in subset_names, spObj.molecules)
        spObj.cell_coord=filter(cell_col => x -> x in subset_names, spObj.cell_coord)
        spObj.mol_coord=filter(cell_col => x -> x in subset_names, spObj.mol_coord)    
    end
    return spObj
end

function subset_SpaObj_cluster(sp::AbstractSpaObj, cluster_col::Union{String, Symbol}, cell_col::Union{String, Symbol}, subset_names::Union{Vector{String}, Vector{Int64}})
    spObj=deepcopy(sp)
    cluster_names = deepcopy(subset_names)
    cell_coord = spObj.cells
    cell_coord = filter(cluster_col => x-> x in cluster_names, cell_coord)
    barcodes = cell_coord[!, cell_col]
    barcodes2 = [["gene"]; barcodes]
    spObj.norm_counts=spObj.norm_counts[!, barcodes2]
    spObj.counts=spObj.counts[!, barcodes2]
    spObj.cells=filter(cell_col => x -> x in barcodes, spObj.cells)
     if isa(sp, SpaObj)
        spObj.cell_names=barcodes
        spObj.molecules=filter(cell_col => x -> x in barcodes, spObj.molecules)
        spObj.cell_coord=filter(cell_col => x -> x in barcodes, spObj.cell_coord)
        spObj.mol_coord=filter(cell_col => x -> x in barcodes, spObj.mol_coord)
    end
    return spObj
end