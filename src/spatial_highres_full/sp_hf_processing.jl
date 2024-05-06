function convert_image_data(sp; layer_slot = "8_um")
    sp = set_default_layer(sp; layer_slot = layer_slot)
    sp = update_coordinates_hd(sp)
    return sp
end

function create_polygon(sp; layer_slot = "8_um", img_res = "low", x_col="x", y_col="y", cell_col = "cell")
        pos = sp.alterImgData.posData.positions[img_res * "_pos"]
        px_width = parse(Int, split(img_res, "_")[1]) / hd.imageData.jsonParameters["microns_per_pixel"]
        corner_coordinates = compute_corner_points(pos, px_width; cell = cell_col, x_col = x_col, y_col = y_col)
        seg = DataFrame(x = corner_coordinates.new_x, y = corner_coordinates.new_y, cell_id = corner_coordinates.barcode)
        grouped = groupby(seg, :cell_id)
        cell_ids = unique(seg.cell_id)
        poly = Vector{Matrix{Float64}}(undef, length(cell_ids))
        n = length(cell_ids)
        println("Formatting cell polygons...")
        p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
        for idx in 1:length(cell_ids)
            cell_data = grouped[idx]
            cell_1 = Matrix(cell_data[!, 1:2])
            poly[idx] = cell_1
            next!(p)
        end
        return poly
end