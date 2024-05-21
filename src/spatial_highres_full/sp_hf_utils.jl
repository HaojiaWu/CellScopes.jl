function set_default_layer(obj::VisiumHDObject; layer_slot::Union{String, Nothing}=nothing)
    layer = get(obj.layerData.layers, obj.defaultData, nothing)
    if layer !== nothing
        layer.rawCount = obj.rawCount
        layer.normCount = obj.normCount
        layer.scaleCount = obj.scaleCount
        layer.metaData = obj.metaData
        layer.spmetaData = obj.spmetaData
        layer.varGene = obj.varGene
        layer.dimReduction = obj.dimReduction
        layer.clustData = obj.clustData
        if !isa(obj.alterImgData, Nothing)
            layer.posData = obj.alterImgData.posData
            layer.polygonData.polygons = merge(Dict("original" => obj.polygonData), obj.alterImgData.polyData.polygons)
        end
    end
    obj.defaultData = layer_slot
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
        obj.polygonData = layer.polygonData.polygons["original"]
        obj.imageData.jsonParameters  = layer.jsonParameters
        if !isa(obj.alterImgData, Nothing)
            obj.alterImgData.posData = layer.posData
            obj.alterImgData.polyData.polygons = Dict(key => value for (key, value) in layer.polygonData.polygons if key != "original")
        end
    end
    return obj
end

default_layer(obj::VisiumHDObject)=obj.defaultData
list_layers(obj::VisiumHDObject) = println(keys(obj.layerData.layers))

function set_python_environment(python_path::Union{String, Nothing})
    if python_path !== Nothing
        ENV["PYTHON"] = python_path
    end
    Pkg.build("PyCall")
    @eval begin
        using PyCall
    end
end

function read_parquet_py(parquet_file)
    pd = pyimport("pandas")
    parquet_df = pd.read_parquet(parquet_file)
    parquet_df = pd_to_df(parquet_df)
    return parquet_df
end

function read_parquet(parquet_file)
    db = DuckDB.open(":memory:")
    con = DuckDB.connect(db)
    res = DuckDB.execute(con, "SELECT * FROM '$parquet_file';")
    pos_data = DataFrame(res)
    return pos_data
end

function read_hd_pos(pos_file)
    positions = read_parquet(pos_file)
    rename!(positions, ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"])
    positions.pxl_col_in_fullres = Float64.(positions.pxl_col_in_fullres)
    positions.pxl_row_in_fullres = Float64.(positions.pxl_row_in_fullres)
    positions.array_col = Int64.(positions.array_col)
    positions.array_row = Int64.(positions.array_row)
    positions.barcode = string.(positions.barcode)
    return positions
end

function create_categorical_colormap(values, alpha)
    unique_values = sort(unique(values)) 
    colors = [RGBA(1, 1, 1, 1)]  
    append!(colors, [RGBA(c.r, c.g, c.b, alpha) for c in distinguishable_colors(length(unique_values) - (0 in unique_values ? 1 : 0))])
    color_map = Dict(zip(unique_values, colors))
    return color_map
end

function assign_integers(strings::Vector{String})
    string_to_int = Dict{String, Int}()
    counter = 1  
    output = Int[]  

    function is_integer(s::String)
        try
            parse(Int, s)
            return true
        catch
            return false
        end
    end

    for s in strings
        if is_integer(s)
            parsed_int = parse(Int, s)
            push!(output, parsed_int)  
            string_to_int[s] = parsed_int  
        else
            if !haskey(string_to_int, s)
                string_to_int[s] = counter
                counter += 1
            end
            push!(output, string_to_int[s])  
        end
    end

    return string_to_int, output
end

function convert_df_hm(df; x_col="x", y_col="y", anno = "cluster")
    df[!, anno] = assign_integers(df[!, anno])
    df[!, x_col]=Int.(round.(df[!, x_col]))
    df[!, y_col]=Int.(round.(df[!, y_col]))
    max_x = maximum(df[!, x_col])
    max_y = maximum(df[!, y_col])
    matrix1 = zeros(Int, max_x + 1, max_y + 1)  
    for row in eachrow(df)
        matrix1[row[x_col] + 1, row[y_col] + 1] = row[anno]
    end
    return matrix1
end

function get_key(d::Dict, value)
    for (k, v) in d
        if v == value
            return k
        end
    end
    return nothing 
end

function color_to_rgba(color_name::String, alpha::Float64 = 1.0)
    color = parse(Colorant, color_name)
    rgba = RGBA(color.r, color.g, color.b, alpha)
    return rgba
end

function compute_corner_points(df::DataFrame, width::Real; x_col="x", y_col="y", cell= "barcode")
    n = nrow(df)
    df.ID = collect(1:nrow(df))
    corners = DataFrame(ID = Int[], barcode = String[], new_x = Float64[], new_y = Float64[])
    for i in 1:n
        x, y = df[i, x_col], df[i, y_col]
        half_width = width / 2
        push!(corners, [df[i, :ID], df[i, cell], x - half_width, y + half_width])  
        push!(corners, [df[i, :ID], df[i, cell], x - half_width, y - half_width]) 
        push!(corners, [df[i, :ID], df[i, cell], x + half_width, y - half_width]) 
        push!(corners, [df[i, :ID], df[i, cell], x + half_width, y + half_width])  
    end
    return corners
end

function process_hd_coordinates(img, pos, scale_factor; return_img=true)
    max_x = maximum(pos.pxl_row_in_fullres) + 1
    transformed_pos = DataFrames.transform(pos, [:pxl_row_in_fullres, :pxl_col_in_fullres] => ByRow((x, y) -> (y, -x + max_x)) => [:pxl_row_in_fullres, :pxl_col_in_fullres]);
    h1, w1 = size(img)
    df = DataFrame(x = repeat(1:w1, inner = h1),
                   y = repeat(1:h1, outer = w1),
                   color = vec(img))
    transformed_pos.pxl_row_in_fullres1 = transformed_pos.pxl_row_in_fullres .* scale_factor
    transformed_pos.pxl_col_in_fullres1 = transformed_pos.pxl_col_in_fullres .* scale_factor
    df1=df[!, ["x","y"]]
    df1.datatype .= "HE"
    df2=transformed_pos[!, ["pxl_col_in_fullres1","pxl_row_in_fullres1"]]
    rename!(df2, ["x","y"])
    df2.datatype .= "cell"
    df3=[df1; df2]
    df3.y = invert_y_axis(df3.y)
    new_df1 = df3[df3.datatype .== "HE",:]
    new_df1.color = df.color
    new_df1.x=Int64.(round.(new_df1.x))
    new_df1.y=Int64.(round.(new_df1.y))
    max_y = maximum(new_df1.y)
    max_x = maximum(new_df1.x)
    if max_y < h1
        new_df1.y .+= 1
        max_y += 1
    elseif max_y > h1
        new_df1.y .-= 1
        max_y -= 1
    end
    if max_x < w1
        new_df1.x .+= 1
        max_x += 1
    elseif max_x > w1
        new_df1.x .-= 1
        max_x -= 1
    end
    new_img = fill(RGBA(1, 1, 1, 1), max_x, max_y)
    indices = CartesianIndex.(new_df1.x, new_df1.y)
    new_img[indices] = new_df1.color
    pos.pxl_row_in_fullres1 = pos.pxl_row_in_fullres .* scale_factor
    pos.pxl_col_in_fullres1 = pos.pxl_col_in_fullres .* scale_factor
    df4=pos[!, ["pxl_col_in_fullres1","pxl_row_in_fullres1"]]
    rename!(df4, ["x","y"])
    df4.datatype .= "cell"
    df5=[df1; df4];
    df5.y = invert_y_axis(df5.y)
    new_pos = df5[df5.datatype .== "cell",:]
    new_pos.cell = pos.barcode
    new_pos = new_pos[!, [:cell, :x, :y]]
    new_pos.x = new_pos.x ./ scale_factor
    new_pos.y = new_pos.y ./ scale_factor
    if return_img
        return new_img, new_pos
    else
        return new_pos
    end
end

function update_coordinates_hd(sp::VisiumHDObject)
    alter_imgdata = AlterImages()
    pos_data = Positions()
    poly_data = Polygons()
    layer_slot = sp.defaultData
    px_width = parse(Int, Base.split(layer_slot, "_")[1]) / sp.imageData.jsonParameters["microns_per_pixel"]
    low_res = deepcopy(sp.imageData.lowresImage)
    sp_meta = deepcopy(sp.spmetaData)
    scale_factor = get_vs_sf(sp; img_res = "low")
    img1, pos1 = process_hd_coordinates(low_res, sp_meta, scale_factor)
    alter_imgdata.imgs["low"] = img1
    pos_data.positions["low_pos"] = pos1
    poly1 = create_polygon(pos1, px_width; x_col="x", y_col="y", cell_col = "cell")
    poly_data.polygons["low_poly"] = poly1
    sp_meta = deepcopy(sp.spmetaData)
    hi_res = deepcopy(sp.imageData.highresImage)
    scale_factor = get_vs_sf(sp; img_res = "high")
    img2, pos2 = process_hd_coordinates(hi_res, sp_meta, scale_factor; return_img=true)
    poly2 = create_polygon(pos2, px_width; x_col="x", y_col="y", cell_col = "cell")
    alter_imgdata.imgs["high"] = img2
    pos_data.positions["high_pos"] = pos2
    poly_data.polygons["high_poly"] = poly2
    if !isa(sp.imageData.fullresImage, Nothing)
        sp_meta = deepcopy(sp.spmetaData)
        full_res = deepcopy(sp.imageData.fullresImage)
        scale_factor = get_vs_sf(sp; img_res = "full")
        img3, pos3 = process_hd_coordinates(full_res, sp_meta, scale_factor)
        poly3 = create_polygon(pos3, px_width; x_col="x", y_col="y", cell_col = "cell")
        alter_imgdata.imgs["full"] = img3
        pos_data.positions["full_pos"] = pos3
        poly_data.polygons["full_poly"] = poly3
    end
    sp.alterImgData = AlterHDImgObject(alter_imgdata, pos_data, poly_data)
    return sp
end

function starts_with(s, patterns)
    for pattern in patterns
        if startswith(s, pattern)
            return true
        end
    end
    return false
end