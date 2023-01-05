abstract type AbstractSpaObj <: AbstractCellScope end

mutable struct imputeObj <: AbstractSpaObj
    tg_data::Union{DataFrame, Nothing}
    spage_data::Union{DataFrame, Nothing}
    gimvi_data::Union{DataFrame, Nothing}
    function imputeObj(imp_type::String; 
        imp_data::Union{DataFrame, Nothing}=nothing
        )
        impute_obj = new(nothing, nothing, nothing)
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
end

mutable struct SpaObj <: AbstractSpaObj
    molecules::DataFrame
    cells::DataFrame
    counts::DataFrame
    norm_counts::Union{DataFrame, Nothing}
    mol_coord::DataFrame
    cell_coord::DataFrame
    gene_names::Union{Vector, String}
    cell_names::Union{Vector, String}
    image::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}}
    polygons::Array{Array{Float64, 2}, 1}
    poly_anno::Union{DataFrame, Nothing}
    poly_counts::Union{DataFrame, Nothing}
    poly_norm::Union{DataFrame, Nothing}
    imp_data::Union{imputeObj, Nothing}
    imp_meta::Union{DataFrame, Nothing}
    sc_data::Union{DataFrame, Nothing}
    sc_meta::Union{DataFrame, Nothing}
    function SpaObj(molecules::DataFrame, cells::DataFrame, counts::DataFrame; 
            cell_prefix::Union{String, Nothing}=nothing, normalize::Bool=true)
        spObj=new(molecules, cells, counts)
        spObj.norm_counts=nothing
        if normalize
            println("Normalizing data...")
            orig_count = deepcopy(counts)
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
            DataFrames.rename!(new_df, DataFrames.names(orig_count)[2:end])
            new_df[!,:gene]=orig_count[!,:gene]
            total_cell=length(names(new_df))
            new_df=new_df[!,[total_cell; collect(1:total_cell-1)]]
            spObj.norm_counts = new_df
            println("Normalization done!")
        end
        spObj.gene_names=unique(molecules.gene)
        cellnames=unique(cells.cell)
        if cell_prefix !== nothing
            println("Adding prefix " * cell_prefix * " to all cells...")
            cellnames=cell_prefix * "_" .* names(counts)[2:end]
            cellnames2=deepcopy(cellnames)
            rename!(counts, append!(["gene"], cellnames))
            if typeof(cells.cell[1]) !== String
                cells[!, :cell2]=string.(cells.cell)
            end
            cells.cell2=cell_prefix * "_" .* cells.cell2
            if typeof(molecules.cell[1]) !== String
                molecules[!,:cell2]=string.(molecules.cell)
            end
            molecules.cell2=cell_prefix * "_" .* molecules.cell2
            println("Cell names with " * cell_prefix * " added!")
            if spObj.norm_counts !==nothing
               rename!(spObj.norm_counts, append!(["gene"], cellnames2))
            end
        end
        spObj.cell_names=cellnames
        spObj.mol_coord=spObj.molecules[!, [:gene, :cell,:cell2, :x, :y]]
        spObj.cell_coord=spObj.cells[!, [:cell,:cell2, :x, :y]]
        println("SpaObj was succussfully created!")
        return spObj
    end
end

mutable struct visiumObj <: AbstractSpaObj
    cells::DataFrame
    counts::DataFrame
    norm_counts::Union{DataFrame, Nothing}
    cell_coord::DataFrame
    gene_names::Union{Vector, String}
    cell_names::Union{Vector, String}
    image::Union{Matrix{RGB{N0f8}},Matrix{Gray{N0f8}}}
    imp_data::Union{DataFrame, Nothing}
    imp_meta::Union{DataFrame, Nothing}
    sc_data::Union{DataFrame, Nothing}
    sc_meta::Union{DataFrame, Nothing}
    function visiumObj(visium_dir::String; data_process::String="Seurat", cell_prefix::Union{String, Nothing}=nothing, 
        project_name::String="seuratObj", n_var_genes::Int64=1000,n_pca::Int64=20, res::Union{Float64, Int64}=0.6, min_dist::Union{Float64, Int64}=0.3)
        if data_process==="Seurat"
            println("Please make sure Seurat and its dependencies have been installed in your R environment.")
            visium = run_seurat(visium_dir; data_type="visium", cell_prefix=cell_prefix, n_pca=n_pca, project_name = project_name,
                                n_var_genes=n_var_genes,res=res,min_dist=min_dist)
            cells, raw_counts = retrieve_seurat_data(visium; data_type="visium")
            original_counts = deepcopy(raw_counts)
            norm_counts = normalizeData(original_counts)
            vsmObj = new(cells, raw_counts, norm_counts)
            println("visiumObj was succussfully created!")
            return vsmObj
        elseif data_process==="Scanpy"
            println("Please make sure Scanpy and its dependencies have been installed in your Python environment.")
            visium = run_scanpy(visium_dir; data_type="visium", cell_prefix=cell_prefix, n_pca=n_pca,
                                n_var_genes=n_var_genes,res=res,min_dist=min_dist)
            cells, raw_counts = retrieve_scanpy_data(visium; data_type="visium")
            original_counts = deepcopy(raw_counts)
            norm_counts = normalizeData(original_counts)
            vsmObj = new(cells, raw_counts, norm_counts)
            println("visiumObj was succussfully created!")
            return vsmObj
        else
            error("data_process can only be \"Seurat\" or \"Scanpy\" so far.")
        end
    end
end

mutable struct scRNAObj <: AbstractSpaObj
    cells::DataFrame
    counts::DataFrame
    norm_counts::Union{DataFrame, Nothing}
    function scRNAObj(file_name::String; data_process::String="Seurat", file_format::String="rda")
        if data_process==="Seurat"
            seu = read_seurat_data(file_name; file_format=file_format)
            cells, raw_counts, norm_counts = retrieve_seurat_data(seu; data_type="scRNA")
            scObj = new(cells, raw_counts, norm_counts)
            println("scRNAObj was succussfully created!")
            return scObj
        elseif data_process==="Scanpy"
            adata = read_scanpy_data(file_name)
            cells, raw_counts, norm_counts = retrieve_scanpy_data(adata; data_type="scRNA")
            scObj = new(cells, raw_counts, norm_counts)
            println("scRNAObj was succussfully created!")
            return scObj
        else
            error("data_process can only be \"Seurat\" or \"Scanpy\" so far.")
        end
    end
end

function add_impdata(impute_obj::imputeObj, imp_type::String, imp_data::DataFrame)
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

function normalizeData(sp::SpaObj)
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

function polygons_cell_mapping(sp::SpaObj)
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

function generate_polygon_counts(sp::SpaObj)
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

function subset_SpaObj(sp::Union{SpaObj, visiumObj, scRNAObj}, cell_col::Union{String, Symbol}, subset_names::Union{Vector{String}, Vector{Int64}})
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

function subset_SpaObj_cluster(sp::Union{SpaObj, visiumObj, scRNAObj}, cluster_col::Union{String, Symbol}, cell_col::Union{String, Symbol}, subset_names::Union{Vector{String}, Vector{Int64}})
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