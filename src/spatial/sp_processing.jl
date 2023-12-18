function polygons_cell_mapping(sp::AbstractImagingObj; anno::Union{String, Symbol}="cluster")
    polygons = sp.polygonData
    n_polygons = length(polygons)
    cell = Vector{Int}(undef, n_polygons)
    X = Vector{Float64}(undef, n_polygons)
    Y = Vector{Float64}(undef, n_polygons)
    for (i, p) in enumerate(polygons)
        X[i], Y[i] = mean(p[:,1]), mean(p[:,2])
        cell[i] = i
    end
    center_df = DataFrame(cell=cell, X=X, Y=Y)
    cell_coord = sp.spmetaData.cell

    poly_mtx = Matrix(center_df[!, 2:3])
    ref_mtx = Matrix(cell_coord[:, ["x","y"]])

    ref_coord = [row for row in eachrow(ref_mtx)]
    query_coord = [row for row in eachrow(poly_mtx)]
    query_dist = Vector{Int}(undef, length(query_coord))
    p = Progress(length(query_coord), dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    for (j, q) in enumerate(query_coord)
        ref_dists = [Distances.euclidean(q, r) for r in ref_coord]
        query_dist[j] = argmin(ref_dists)
        next!(p)
    end

    println("Polygons have been mapped to the cells!")

    my_annot = DataFrame(polygon_number=collect(1:length(query_dist)), mapped_cell=query_dist)
    from = sp.spmetaData.cell.cell
    to = sp.spmetaData.cell[!, anno]

    prefix = Base.split(sp.spmetaData.cell.cell[1], "_")
    if length(prefix) > 1
        my_annot.mapped_cell .= prefix[1] .* "_" .* string.(my_annot.mapped_cell)
    end
    
    anno = map_values(my_annot, :mapped_cell, :cluster, from, to)
    sp.spmetaData.polygon = anno

    return sp
end

#=
function polygons_cell_mapping(sp::AbstractImagingObj; anno::Union{String, Symbol}="cluster")
    polygons=sp.polygonData
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
    cell_coord=sp.spmetaData.cell
    poly_mtx = Matrix(center_df[!,2:3])
    ref_mtx = Matrix(cell_coord[:,["x","y"]])
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
    for j in 1:length(query_coord)
        ref_dist=[Distances.euclidean(query_coord[j],ref_coord[i]) for i in 1:length(ref_coord)]
        query_dist=push!(query_dist,findmin(ref_dist)[2])
        next!(p)
    end
    println("Polygons have been mapped to the cells!")
    my_annot=DataFrame(polygon_number=collect(1:length(query_dist)), mapped_cell=query_dist);
    from=sp.spmetaData.cell.cell
    to=sp.spmetaData.cell[!, anno]
    prefix = split(sp.spmetaData.cell.cell[1], "_")
    if length(prefix) > 1
        my_annot.mapped_cell = prefix[1] * "_" .* string.(my_annot.mapped_cell);
    end
    anno=map_values(my_annot, :mapped_cell, :cluster,from, to)
    sp.spmetaData.polygon=anno
    return sp
end
=#

function generate_polygon_counts(sp::AbstractImagingObj)
    coord_molecules = deepcopy(sp.spmetaData.molecule)
    if isa(sp.spmetaData.polygon, Nothing)
        error("Please run polygons_cell_mapping first!")
    end
    anno = deepcopy(sp.spmetaData.polygon)
    coord_molecules.cell = string.(coord_molecules.cell)
    anno = rename!(anno, :polygon_number => :number)
    anno = rename!(anno, :mapped_cell => :cell)
    anno.cell = string.(anno.cell)
    join_df = innerjoin(coord_molecules[!,[:gene,:cell]],anno[!,[:number,:cell]], on=:cell)
    join_df.gene = string.(join_df.gene)
    gdf = groupby(join_df, :gene)
    new_df = []
    p = Progress(length(gdf), dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    if isa(sp, StereoSeqObject)
        for df in gdf
            anno1 = DataFrames.combine(groupby(df, :number), :MIDCount => sum => :count)
            anno1.gene .= df.gene[1]
            push!(new_df, anno1)
            next!(p)
        end
        new_df = vcat(new_df...)
        rename!(new_df, :number => :cell)    
    else
        for df in gdf
            map_dict = countmap(df.number)
            anno1 = DataFrame(cell=collect(keys(map_dict)), count=collect(values(map_dict)))
            anno1.gene .= df.gene[1]
            push!(new_df, anno1)
            next!(p)
        end
        new_df = vcat(new_df...) 
    end
    final_df = unstack(new_df, :cell, :gene, :count)
    final_df .= ifelse.(isequal.(final_df, missing), 0, final_df)
    final_df = mapcols(ByRow(Int64), final_df);
    sort!(final_df, :cell)
    prefix = Base.split(sp.spmetaData.cell.cell[1], "_")
    if length(prefix) > 1
        final_df.cell = prefix[1] * "_" .* string.(final_df.cell)
    end
    cellnames = final_df.cell
    genenames = names(final_df)[2:end]
    final_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(final_df[!,2:end])')
    raw_count = RawCountObject(final_df, cellnames, genenames)
    norm_count = normalize_object(raw_count)
    sp.polyCount = raw_count
    sp.polynormCount = norm_count
    println("Polygon counts was normalized!")
    return sp
end

function process_scs_file(file_path, add_x, add_y, prefix)
    df = CSV.File(file_path; delim='\t', header=false) |> DataFrame
    rename!(df, [:coord, :cell])
    coords = Base.split.(df[:, :coord], ':')
    df.x = parse.(Int, getindex.(coords, 1)) .+ add_x
    df.y = parse.(Int, getindex.(coords, 2)) .+ add_y
    df.cell = prefix .* "_" .* string.(df.cell)
    return df
end

function process_scs_directory(directory)
    filenames = glob("spot2cell_*.txt", directory)
    all_dataframes = DataFrame[]
    p = Progress(length(filenames), desc="Processing Files")
    for filename in filenames
        numbers = [parse(Int, m.match) for m in eachmatch(r"\d+", basename(filename))]
        add_x, add_y = numbers[2], numbers[3]
        file_path = joinpath(directory, filename)
        prefix = join(numbers[2:4], ':')
        df = process_scs_file(file_path, add_x, add_y, prefix)
        push!(all_dataframes, df)
        next!(p)
    end
    return vcat(all_dataframes...)
end

function create_stereoseq_scs(scs_results, spot_coord; prefix="sp", min_gene=0, min_cell=0)
    @info("1. Reading input files...")
    final_dataframe = process_scs_directory(scs_results)
    orig_cord = CSV.File(spot_coord; delim='\t') |> DataFrame
    final_dataframe[!, :spot_loc] = [string(i) * "_" * string(j) for (i, j) in zip(final_dataframe.x, final_dataframe.y)]
    orig_cord[!, :spot_loc] = [string(i) * "_" * string(j) for (i, j) in zip(orig_cord.x, orig_cord.y)]
    orig_cord2 = filter(:spot_loc => ∈(Set(final_dataframe.spot_loc)), orig_cord)
    new_coord = innerjoin(orig_cord2, final_dataframe[:, [:spot_loc, :cell]], on=:spot_loc)
    gdf = groupby(new_coord, :cell)
    new_df = []
    @info("2. Creating count matrix...")
    p = Progress(length(gdf), dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    for df in gdf
        anno1 = DataFrames.combine(groupby(df, :geneID), :MIDCount => sum => :count)
        anno1.cell .= df.cell[1]
        push!(new_df, anno1)
        next!(p)
    end
    new_df = vcat(new_df...)
    @info("3. Constructing StereoSeq object...")
    final_df = unstack(new_df, :cell, :geneID, :count)
    final_df .= ifelse.(isequal.(final_df, missing), 0, final_df)
    cellnames = final_df.cell
    final_df = mapcols(ByRow(Int64), final_df[!, 2:end])
    genenames = names(final_df)
    final_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(final_df)')
    raw_count = RawCountObject(final_df, cellnames, genenames)
    cell_coord = DataFrames.combine(groupby(new_coord, :cell), :x => mean => :x, :y => mean => :y)
    rename!(new_coord,:geneID => :gene)
    sp = StereoSeqObject(new_coord, cell_coord, raw_count;
    prefix = prefix, min_gene = min_gene, min_cell = min_cell)
    return sp
end
