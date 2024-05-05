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
        obj.imageData.jsonParameters  = layer.jsonParameters
        if !isa(obj.alterImgData, Nothing)
            obj.alterImgData.posData = layer.posData
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

function create_image(df)
    max_y = maximum(df.y)
    max_x = maximum(df.x)
    new_img = fill(RGBA(1, 1, 1, 1), max_x, max_y)
    x_coords = df.x
    y_coords = df.y
    colors = df.color
    indices = CartesianIndex.(df.x, df.y)
    new_img[indices] = df.color
    return new_img
end

function process_hd_coordinates(img, pos, scale_factor; return_img=true)
    h1, w1 = size(img)
    df = DataFrame(x = repeat(1:w1, inner = h1),
               y = repeat(1:h1, outer = w1),
               color = vec(img))
    pos.pxl_row_in_fullres = pos.pxl_row_in_fullres .* scale_factor
    pos.pxl_col_in_fullres = pos.pxl_col_in_fullres .* scale_factor
    df2=df[!, ["x","y"]]
    df2.datatype .= "HE"
    df3=pos[!, ["pxl_col_in_fullres","pxl_row_in_fullres"]]
    rename!(df3, ["x","y"])
    df3.datatype .= "cell"
    df4=[df2; df3]
    df4.y = invert_y_axis(df4.y)
    df_grp = groupby(df4, :datatype)
    new_pos = DataFrame(cell = pos.barcode, x=df_grp[2].x ./ scale_factor, y=df_grp[2].y ./ scale_factor)
    if return_img
        df2.color = df.color
        new_img = create_image(df2)
        return new_img, new_pos
    else
        return new_pos
    end
end

function update_coordinates_hd(sp::VisiumHDObject)
    low_res = deepcopy(sp.imageData.lowresImage)
    sp_meta = deepcopy(sp.spmetaData)
    scale_factor = get_vs_sf(sp; img_res = "low")
    img1, pos1 = process_hd_coordinates(low_res, sp_meta, scale_factor)
    alter_imgdata = AlterHDImgObject(nothing, nothing)
    alter_imgdata.imgData["low"] = img1
    alter_imgdata.posData["low_pos"] = pos1
    sp_meta = deepcopy(sp.spmetaData)
    hi_res = deepcopy(sp.imageData.highresImage)
    scale_factor = get_vs_sf(sp; img_res = "high")
    img2, pos2 = process_hd_coordinates(hi_res, sp_meta, scale_factor; return_img=true)
    alter_imgdata.imgData["high"] = img2
    alter_imgdata.posData["high_pos"] = pos2
    if !isa(sp.imageData.fullresImage, Nothing)
        sp_meta = deepcopy(sp.spmetaData)
        full_res = deepcopy(sp.imageData.fullresImage)
        scale_factor = get_vs_sf(sp; img_res = "full")
        img3, pos3 = process_hd_coordinates(full_res, sp_meta, scale_factor)
        alter_imgdata.imgData["full"] = img3
        alter_imgdata.posData["full_pos"] = pos3
    end
    sp.alterImgData = alter_imgdata
    return sp
end
