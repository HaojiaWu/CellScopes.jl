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
    end
    return obj
end

function default_layer(obj::VisiumHDObject)
    return obj.defaultData
end

function set_python_environment(python_path::Union{String, Nothing})
    if python_path !== Nothing
        ENV["PYTHON"] = python_path
    end
    Pkg.build("PyCall")
    @eval begin
        using PyCall
    end
end

function read_parquet(parquet_file; python_path::Union{String, Nothing} = nothing)
    set_python_environment(python_path)
    pd = pyimport("pandas")
    parquet_df = pd.read_parquet(parquet_file)
    parquet_df = pd_to_df(parquet_df)
    return parquet_df
end

function read_hd_pos(pos_file; python_path::Union{String, Nothing} = nothing)
    positions = read_parquet(pos_file; python_path=python_path)
    rename!(positions, ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"])
    if !isa(positions.pxl_col_in_fullres, Vector{<:Real})
        positions.pxl_col_in_fullres = parse.(Int64, positions.pxl_col_in_fullres)
    end
    if !isa(positions.pxl_row_in_fullres, Vector{<:Real})
        positions.pxl_row_in_fullres = parse.(Int64, positions.pxl_row_in_fullres)
    end
    if !isa(positions.array_col, Vector{<:Real})
        positions.array_col = parse.(Int64, positions.array_col)
    end
    if !isa(positions.array_row, Vector{<:Real})
        positions.array_row = parse.(Int64, positions.array_row)
    end
    positions.barcode = string.(positions.barcode)
    return positions
end

