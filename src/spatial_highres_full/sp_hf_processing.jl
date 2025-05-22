function convert_image_data(sp; layer_slot = "8_um")
    sp = set_default_layer(sp; layer_slot = layer_slot)
    sp = update_coordinates_hd(sp)
    return sp
end

function create_polygon(pos::DataFrame, px_width; x_col="x", y_col="y", cell_col = "cell")
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

function parse_molecule(hd_dir)
    data = HDF5.h5read(hd_dir * "/molecule_info.h5","count")
    data = Int64.(data)
    feature_ID = HDF5.h5read(hd_dir * "/molecule_info.h5","feature_idx")
    feature_ID = Int64.(feature_ID) .+ 1
    gene_name = HDF5.h5read(hd_dir * "/molecule_info.h5","features/name")
    barcode_id = HDF5.h5read(hd_dir * "/molecule_info.h5","barcode_idx")
    barcode_id = Int64.(barcode_id) .+ 1
    barcode = HDF5.h5read(hd_dir * "/molecule_info.h5","barcodes")
    barcode_all = [barcode[i] for i in barcode_id]
    feature_all = [gene_name[i] for i in feature_ID]
    umi_type = HDF5.h5read(hd_dir * "/molecule_info.h5","umi_type")
    umi_type = Int64.(umi_type)
    transcript = DataFrame(:gene => feature_all, :count => data, :barcode2 => barcode_all, :umi_type => umi_type)
    pos = read_parquet(hd_dir * "/binned_outputs/square_002um/spatial/tissue_positions.parquet")
    filter!(:in_tissue => !=(0), pos)
    pos[:, :barcode2] = [replace(s, r"-\d+" => "") for s in pos.barcode]
    filter!(:barcode2 => ∈(Set(pos.barcode2)), transcript)
    transcript = leftjoin(transcript, pos, on = :barcode2)
    transcript = transcript[!, [:barcode, :gene, :pxl_row_in_fullres, :pxl_col_in_fullres, :count, :umi_type]]
    return transcript
end