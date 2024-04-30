function read_10x(tenx_dir::String; 
    version::String="v2", 
    min_gene::Int64 = 0, 
    min_cell::Int64 = 0
)
    if version === "v2"
        counts = MatrixMarket.mmread(tenx_dir * "/matrix.mtx")
        cells = CSV.File(tenx_dir * "/barcodes.tsv", header = false) |> DataFrame
        genes = CSV.File(tenx_dir * "/genes.tsv", header = false) |> DataFrame
    elseif version === "v3"
        gene_file = tenx_dir * "/features.tsv.gz"
        cell_file = tenx_dir * "/barcodes.tsv.gz"
        count_file = tenx_dir * "/matrix.mtx.gz"
        genes = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", gene_file); header_exists=false))
        cells = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", cell_file); header_exists=false))
        counts = MatrixMarket.mmread(gunzip(count_file))
    else
        error("version can only be v2 or v3!")
    end
    cells = string.(cells.Column1)
    genes = string.(genes.Column2)
    gene_kept = (vec ∘ collect)(rowSum(counts).> min_cell)
    genes = genes[gene_kept]
    cell_kept = (vec ∘ collect)(colSum(counts) .> min_gene)
    cells = cells[cell_kept]
    gene_indices = findall(gene_kept)
    cell_indices = findall(cell_kept)
    counts = counts[gene_indices, cell_indices]
    gene_kept, gene_removed = check_duplicates(genes)
    gene_removed = collect(values(gene_removed))
    total_num = collect(1:length(genes))
    gene_num = total_num[Not(gene_removed)]
    to_keep = [i ∈ gene_num for i in total_num]
    counts = counts[findall(to_keep), :]
    rawcount = RawCountObject(counts, cells, gene_kept)
    return rawcount
end

function save(sc_obj::AbstractCellScope; filename::String = "cs_obj.jld2")
    JLD2.save(filename, "key", sc_obj)
end

function load(;filename::String = "cs_obj.jld2")
    cs = JLD2.load(filename)
    cs = cs["key"]
    return cs
end

function read_visium(visium_dir::String; 
    min_gene::Int64 = 0, 
    min_cell::Int64 = 0
)
    # locate all files
    highres_image_file = visium_dir * "/spatial/tissue_hires_image.png"
    lowres_image_file = visium_dir * "/spatial/tissue_lowres_image.png"
    fullres_image_file = visium_dir * "/spatial/tissue_fullres_image.png"
    detected_tissue = visium_dir * "/spatial/detected_tissue_image.jpg"
    aligned_image_file = visium_dir * "/spatial/aligned_fiducials.jpg"
    position_file = visium_dir * "/spatial/tissue_positions_list.csv"
    if !isfile(position_file)
        position_file = visium_dir * "/spatial/tissue_positions.csv"
    end
    json_file = visium_dir * "/spatial/scalefactors_json.json"
    gene_file = visium_dir * "/filtered_feature_bc_matrix/features.tsv.gz"
    cell_file = visium_dir * "/filtered_feature_bc_matrix/barcodes.tsv.gz"
    count_file = visium_dir * "/filtered_feature_bc_matrix/matrix.mtx.gz"
    # prepare counts
    genes = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", gene_file); header_exists=false))
    genes = string.(genes.Column2)
    cells = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", cell_file); header_exists=false))
    cells = string.(cells.Column1)
    counts = MatrixMarket.mmread(gunzip(count_file))
    gene_kept = (vec ∘ collect)(rowSum(counts).> min_cell)
    genes = genes[gene_kept]
    cell_kept = (vec ∘ collect)(colSum(counts) .> min_gene)
    cells = cells[cell_kept]
    gene_indices = findall(gene_kept)
    cell_indices = findall(cell_kept)
    counts = counts[gene_indices, cell_indices]
    gene_kept, gene_removed = check_duplicates(genes)
    gene_removed = collect(values(gene_removed))
    total_num = collect(1:length(genes))
    gene_num = total_num[Not(gene_removed)]
    to_keep = [i ∈ gene_num for i in total_num]
    counts = counts[findall(to_keep), :]
    rawcount = RawCountObject(counts, cells, gene_kept)
    # prepare spatial info
    positions = DataFrame(CSV.File(position_file, header=true))
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
    high_img = FileIO.load(highres_image_file)
    low_img = FileIO.load(lowres_image_file)
    if isfile(fullres_image_file)
        full_img = FileIO.load(fullres_image_file)
        full_img = convert(Matrix{RGB{N0f8}}, full_img)
    else
        full_img = nothing
    end
    tissue_img = FileIO.load(detected_tissue)
    aligned_img = FileIO.load(aligned_image_file)
    json_data = JSON.parsefile(json_file)
    image_obj = VisiumImgObject(high_img, low_img, full_img, tissue_img, aligned_img, json_data)
    # create visiumobject
    vsm_obj = VisiumObject(rawcount)
    vsm_obj.spmetaData = positions
    vsm_obj.imageData = image_obj
    vsm_obj.spmetaData = vsm_obj.spmetaData[indexin(vsm_obj.metaData.Cell_id, vsm_obj.spmetaData.barcode),:]
    return vsm_obj
end

function read_xenium(xenium_dir::String; prefix = "xenium", min_gene::Int64 = 0, min_cell::Int64 = 0, version="1.1")
    gene_file = xenium_dir * "/cell_feature_matrix/features.tsv.gz"
    cell_file = xenium_dir * "/cell_feature_matrix/barcodes.tsv.gz"
    count_file = xenium_dir * "/cell_feature_matrix/matrix.mtx.gz"
    cell_meta = xenium_dir * "/cells.csv.gz"
    transcript_meta = xenium_dir * "/transcripts.csv.gz"
    seg_file = xenium_dir * "/cell_boundaries.csv.gz"
    cluster_file = xenium_dir * "/analysis/clustering/gene_expression_graphclust/clusters.csv"
    umap_file = xenium_dir * "analysis/umap/gene_expression_2_components/projection.csv"
    xenium_umap =  DataFrame(CSV.File(umap_file))
    genes = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", gene_file); header_exists=false))
    cells = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", cell_file); header_exists=false))
    clustering =  DataFrame(CSV.File(cluster_file))
    cell_ids_filter = clustering.Barcode
    seg = DataFrame(CSVFiles.load(CSVFiles.File(format"CSV", seg_file); header_exists=true))
    cell_ids = unique(seg.cell_id)
    cell_kept = check_vec(cell_ids_filter, cell_ids)
    count_molecules =  DataFrame(CSV.File(transcript_meta))
    count_cells =  DataFrame(CSV.File(cell_meta))
    counts = MatrixMarket.mmread(gunzip(count_file))
    counts = counts[:, cell_kept]
    seg = filter(:cell_id => ∈(Set(clustering.Barcode)), seg)
    cells = filter(:Column1 => ∈(Set(clustering.Barcode)), cells)
    count_molecules = filter(:cell_id => ∈(Set(clustering.Barcode)), count_molecules)
    count_cells = filter(:cell_id => ∈(Set(clustering.Barcode)), count_cells)
    grouped = groupby(seg, :cell_id)
    cell_ids = unique(seg.cell_id)
    poly = Vector{Matrix{Float64}}(undef, length(cell_ids))
    n = length(cell_ids)
    println("Formatting cell polygons...")
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    for idx in 1:length(cell_ids)
        cell_data = grouped[idx]
        cell_1 = Matrix(cell_data[!, 2:end])
        poly[idx] = cell_1
        next!(p)
    end
    println("Cell polygons formatted!")
    xenium_umap.Barcode = string.(xenium_umap.Barcode)
    xenium_umap = Matrix(xenium_umap[!, 2:end])
    umap_obj = UMAPObject(xenium_umap, "UMAP", 2, nothing, nothing, nothing, nothing, nothing)
    genes = string.(genes.Column2)
    blank_gene = Grep.grep("BLANK", genes)
    neg_gene = Grep.grep("NegControl", genes)
    antisense_gene = Grep.grep("antisense", genes)
    unassigned_gene = Grep.grep("Unassigned", genes)
    gene_rm = [blank_gene; neg_gene; antisense_gene; unassigned_gene]
    cells = string.(cells.Column1)
    rename!(count_molecules, :cell_id => :cell, :feature_name => :gene, :x_location => :x, :y_location => :y, :z_location => :z)
    rename!(count_cells, :cell_id => :cell, :x_centroid => :x, :y_centroid => :y)
    if version == 1.0
        count_molecules.cell[count_molecules.cell.==-1] .= 0
    end
    raw_count = RawCountObject(counts, cells, genes)
    genes2 = setdiff(genes, gene_rm)
    println("Filtering cells and genes...")
    raw_count = subset_count(raw_count; genes=genes2)
    println("Cells and genes filtered!")
    count_molecules.cell = string.(count_molecules.cell)
    count_molecules = filter(:gene => ∈(Set(genes2)), count_molecules)
    count_cells.cell = string.(count_cells.cell)
    count_cells = filter(:cell => ∈(Set(string.(clustering.Barcode))), count_cells)
    count_cells.cluster = clustering.Cluster
    spObj = XeniumObject(count_molecules, count_cells, raw_count;
            prefix = prefix, min_gene = min_gene, min_cell = min_gene)
    clustering.cell = prefix .*  "_" .* string.(clustering.Barcode)
    all_cells = spObj.rawCount.cell_name
    clustering = filter(:cell=> ∈(Set(all_cells)), clustering)
    spObj.spmetaData.cell.cluster = string.(clustering.Cluster)
    spObj.metaData.cluster = string.(clustering.Cluster)
    spObj.spmetaData.molecule = map_values(spObj.spmetaData.molecule, :cell, :cluster, clustering.cell, clustering.Cluster)
    spObj.spmetaData.molecule.cluster = string.(spObj.spmetaData.molecule.cluster)
    polygon_df = DataFrame(polygon_number = 1:length(poly), mapped_cell = count_cells.cell, cluster=count_cells.cluster)
    spObj.spmetaData.polygon = polygon_df
    spObj.polygonData = poly
    reduct_obj = ReductionObject(nothing, nothing, umap_obj)
    spObj.dimReduction = reduct_obj
    spObj = normalize_object(spObj)
    spObj.polynormCount = spObj.normCount
    return spObj
end

function read_atac_count(atac_path::String; 
    min_peak::Int64 = 0, 
    min_cell::Int64 = 0
)
    peak_loc = atac_path * "/filtered_peak_bc_matrix/peaks.bed"
    cell_file = atac_path * "/filtered_peak_bc_matrix/barcodes.tsv"
    count_file = atac_path * "/filtered_peak_bc_matrix/matrix.mtx"
    peaks = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", peak_loc); header_exists=false))
    peak_names = string.(peaks.Column1) .* "_" .* string.(peaks.Column2) .* "_" .* string.(peaks.Column3)
    cells = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", cell_file); header_exists=false))
    cells = cells.Column1
    counts = MatrixMarket.mmread(count_file);
    gene_kept = (vec ∘ collect)(rowSum(counts).> min_cell)
    peak_names = peak_names[gene_kept]
    cell_kept = (vec ∘ collect)(colSum(counts) .> min_peak)
    cells = cells[cell_kept]
    counts = counts[gene_kept, cell_kept]
    rawcount = RawCountObject(counts, cells, peak_names)
    return rawcount
end

function read_atac(atac_path; min_peak::Int64=0, min_cell::Int64=0)
    println("This step reads all information directly from cellranger-atac output for downstream analysis. It may take 10 - 15 mins to complete as certain files (e.g. the fragment file) can be large in size.")
    println("1/3 Reading peak count data...")
    raw_count = read_atac_count(atac_path; min_peak=min_peak, min_cell=min_cell)
    atac_obj = scATACObject(raw_count)
    println("1/3 Peak count was loaded!")
    println("2/3 Reading the peak annotation file...")
    peak_anno_file = atac_path * "/peak_annotation.tsv"
    peak_anno = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", peak_anno_file);header_exists=true))
    rename!(peak_anno, "end" => "stop") # end is a Julia key word.
    peak_anno.peak_names = string.(peak_anno.chrom) .* "_" .* string.(peak_anno.start) .* "_" .* string.(peak_anno.stop)
    cells = colnames(atac_obj)
    println("2/3 Peak annotation was loaded!")
    println("3/3 Reading the fragment file...")
    fragment_file = atac_path * "/fragments.tsv.gz"
    fragments = CSV.read(fragment_file, DataFrame; delim = "\t", comment = "#", header =false)
    println("3/3 Fragments were loaded!")
    fragments = filter(:Column4 => ∈(Set(cells)), fragments)
    atac_obj.peakAnno = peak_anno
    frag_data = FragmentObject(fragments, nothing)
    atac_obj.fragmentData = frag_data
    println("scATACObject was successfully constructed!")
    return atac_obj
end

function read_merfish(merfish_dir::String; prefix = "merfish", min_gene::Int64 = 0, min_cell::Int64 = 0)
    cell_meta = merfish_dir * "/cell_metadata.csv"
    transcript_meta = merfish_dir * "/detected_transcripts.csv"
    count_file = merfish_dir * "/cell_by_gene.csv"
    println("1.Loading transcript file...")
    count_molecules = CSV.read(transcript_meta, DataFrame)
    println("1.Transcript file loaded...")
    println("2.Reading cell polygons data...")
    cell_boundary_path = merfish_dir * "/" * "cell_boundaries/"
    all_file = Glob.glob("*.hdf5", cell_boundary_path)
    first_char = length(cell_boundary_path) + length("feature_data_") + 1
    seg = DataFrame()
    n = length(all_file)
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    for i in 1:length(all_file)
        last_char = length(all_file[i]) - length(".hdf5")
        fov_id = parse(Int, all_file[i][first_char:last_char])
        fid = h5open(all_file[i], "r")
        featuredata = fid["featuredata"]
        all_keys = keys(featuredata)
        df1 = DataFrame()
        for j in 1:length(all_keys)
            coord = read(featuredata[all_keys[j]]["zIndex_3"]["p_0"]["coordinates"])
            coord2d = coord[:,:]'
            coord2d = DataFrame(coord2d, :auto)
            rename!(coord2d, :x1 => :x, :x2 => :y)
            coord2d.cell .= all_keys[j]
            coord2d.fov .= fov_id
            df1 = [df1; coord2d]
        end
        seg = [seg; df1]
        close(fid)
        next!(p)
    end
    println("2.Cell polygons loaded!")
    grouped = groupby(seg, :cell)
    cell_ids = unique(seg.cell)
    poly = Vector{Matrix{Float64}}(undef, length(cell_ids))
    n = length(cell_ids)
    println("3.Formatting cell polygons...")
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    for idx in 1:length(cell_ids)
        cell_data = grouped[idx]
        cell_1 = Matrix(cell_data[!, 1:2])
        poly[idx] = cell_1
        next!(p)
    end
    println("3.Cell polygons formatted!")
    cells = unique(seg.cell)
    count_cells = CSV.read(cell_meta, DataFrame; types=Dict(1=>String))
    rename!(count_cells, :Column1 => :cell, :center_x => :x, :center_y => :y)
    count_cells = filter(:cell => ∈(Set(cells)),  count_cells)
    counts = CSV.read(count_file, DataFrame;  types=Dict(1=>String))
    counts = counts[indexin(cells, counts.cell),:]
    count_cells = count_cells[indexin(cells, count_cells.cell),:]
    cells = counts.cell
    genes = names(counts)[2:end]
    gene_rm = Grep.grep("Blank", genes)
    genes2 = setdiff(genes, gene_rm)
    counts = counts[!, genes2]
    counts = convert(SparseMatrixCSC{Int64, Int64},Matrix(counts))
    counts = counts'
    raw_count = RawCountObject(counts, cells, genes2)
    count_molecules = filter(:gene => ∈(Set(genes2)), count_molecules)
    spObj = MerfishObject(count_molecules, count_cells, raw_count, poly;
            prefix = prefix, min_gene = min_gene, min_cell = min_gene)          
    return spObj

end

function read_slideseq(bead_coord_file, count_file; min_gene::Int64 = 0, min_cell::Int64 = 0)
    loc = CSV.read(bead_coord_file, DataFrame; delim=",")
    rename!(loc, ["cell","x","y"])
    counts = CSV.read(count_file, DataFrame; delim="\t")
    gene_name = counts.GENE
    cell_name = string.(names(counts)[2:end])
    counts = counts[!, 2:end]
    counts = convert(SparseMatrixCSC{Int64, Int64},Matrix(counts))
    raw_count = RawCountObject(counts, cell_name, gene_name)
    slide = SlideseqObject(raw_count; sp_meta=loc, min_gene=min_gene, min_cell=min_cell)
    return slide
end

function read_10x_h5(h5file_path)
    data = h5read(h5file_path, "matrix/data")
    indices = h5read(h5file_path, "matrix/indices") .+ 1
    indptr = h5read(h5file_path, "matrix/indptr") .+ 1
    shape = h5read(h5file_path, "matrix/shape")
    genes = h5read(h5file_path, "matrix/features/name")
    cells = h5read(h5file_path, "matrix/barcodes")
    count_mtx = SparseMatrixCSC(shape[1], shape[2], indptr, indices, data)
    count_mtx = sparse(Matrix(count_mtx))
    raw_ct = RawCountObject(count_mtx, cells, genes)
    return raw_ct
end

function read_baysor_loom(loom_path)
    counts = h5read(loom_path, "matrix")
    counts = counts'
    genes = h5read(loom_path, "row_attrs/Name")
    cells = h5read(loom_path, "col_attrs/CellID")
    cells = string.(cells)
    raw_ct = RawCountObject(counts, cells, genes)
    return raw_ct
end

function read_scanpy_loom(loom_path)
    counts = h5read(loom_path, "matrix")
    counts = counts'
    genes = h5read(loom_path, "row_attrs/var_names")
    cells = h5read(loom_path, "col_attrs/obs_names")
    cells = string.(cells)
    raw_ct = RawCountObject(counts, cells, genes)
    return raw_ct
end

function read_seurat_loom(loom_path)
    counts = h5read(loom_path, "matrix")
    counts = counts'
    genes = h5read(loom_path, "row_attrs/Gene")
    cells = h5read(loom_path, "col_attrs/CellID")
    cells = string.(cells)
    raw_ct = RawCountObject(counts, cells, genes)
    return raw_ct
end

function read_baysor(baysor_output::String; tech::String="newtech", 
    prefix::Union{String, Nothing}=nothing, 
    postfix::Union{String, Nothing}=nothing, 
    meta_data::Union{DataFrame, Nothing} = nothing,
    min_gene::Int64=0, min_cell::Int64=0, x_col::Union{String, Symbol} = "x", 
    y_col::Union{String, Symbol} = "y", cell_col::Union{String, Symbol} = "cell")
molecules = CSV.read(baysor_output * "/segmentation.csv", DataFrame)
count_df = CSV.read(baysor_output * "/segmentation_counts.tsv", DataFrame)
cells = CSV.read(baysor_output * "/segmentation_cell_stats.csv", DataFrame)
toml_info = TOML.parsefile(baysor_output * "/segmentation_config.toml")
scale_value = toml_info["Data"]["scale"]
grid_step = scale_value / 7
bandwidth= scale_value / 10
println("Generating cell boundary polygons...")
poly = baysor_boundary_polygons(molecules, molecules[!, cell_col]; grid_step=grid_step, bandwidth=bandwidth)
println("Cell polygons were created!")
gene_name = count_df.gene
cell_name = string.(names(count_df)[2:end])
count_df = count_df[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = RawCountObject(count_df, cell_name, gene_name)
molecules.cell = string.(molecules.cell)
cells.cell = string.(cells.cell)
println("Constructing spatial object...")
if tech == "newtech"
    spObj = ImagingSpatialObject(molecules, cells, raw_count; 
                    prefix=prefix, 
                    postfix=postfix, 
                    meta_data= meta_data,
                    min_gene=min_gene, 
                    min_cell=min_cell, 
                    x_col = x_col, 
                    y_col = y_col, 
                    cell_col = cell_col)
elseif tech == "cartana"
    spObj = CartanaObject(molecules, cells, raw_count; 
                            prefix=prefix, 
                            postfix=postfix, 
                            meta_data= meta_data,
                            min_gene=min_gene, 
                            min_cell=min_cell, 
                            x_col = x_col, 
                            y_col = y_col, 
                            cell_col = cell_col)
elseif tech == "xenium"
    spObj = XeniumObject(molecules, cells, raw_count; 
                            prefix=prefix, 
                            postfix=postfix, 
                            meta_data= meta_data,
                            min_gene=min_gene, 
                            min_cell=min_cell, 
                            x_col = x_col, 
                            y_col = y_col, 
                            cell_col = cell_col)
elseif tech == "merfish"
    spObj = MerfishObject(molecules, cells, raw_count; 
                            prefix=prefix, 
                            postfix=postfix, 
                            meta_data= meta_data,
                            min_gene=min_gene, 
                            min_cell=min_cell, 
                            x_col = x_col, 
                            y_col = y_col, 
                            cell_col = cell_col)
elseif tech == "starmap"
    spObj = STARmapObject(molecules, cells, raw_count; 
                            prefix=prefix, 
                            postfix=postfix, 
                            meta_data= meta_data,
                            min_gene=min_gene, 
                            min_cell=min_cell, 
                            x_col = x_col, 
                            y_col = y_col, 
                            cell_col = cell_col)
elseif tech == "seqfish"
    spObj = seqFishObject(molecules, cells, raw_count; 
                            prefix=prefix, 
                            postfix=postfix, 
                            meta_data= meta_data,
                            min_gene=min_gene, 
                            min_cell=min_cell, 
                            x_col = x_col, 
                            y_col = y_col, 
                            cell_col = cell_col)
elseif tech == "stereoseq"
    spObj = StereoSeqObject(molecules, cells, raw_count; 
                            prefix=prefix, 
                            postfix=postfix, 
                            meta_data= meta_data,
                            min_gene=min_gene, 
                            min_cell=min_cell, 
                            x_col = x_col, 
                            y_col = y_col, 
                            cell_col = cell_col)
elseif tech == "cosmx"
    spObj = CosMxObject(molecules, cells, raw_count; 
                            prefix=prefix, 
                            postfix=postfix, 
                            meta_data= meta_data,
                            min_gene=min_gene, 
                            min_cell=min_cell, 
                            x_col = x_col, 
                            y_col = y_col, 
                            cell_col = cell_col)
else
    error("Please pass the correct spatial technique name to the 'tech' parameter. It can be 'newtech', 'xenium', 'cartana', 'merfish', 'seqfish', 'starmap', 'stereoseq', or 'cosmx', etc. ")
end
return spObj
end

function read_cosmx(cosmx_dir; fov::Union{Nothing, Int64}=nothing, input_type = "fullview",
    prefix::Union{String, Nothing}=nothing, postfix::Union{String, Nothing}=nothing, 
    meta_data::Union{DataFrame, Nothing} = nothing,
    min_gene::Int64=0, min_cell::Int64=0, x_col::Union{String, Symbol} = "x", 
    y_col::Union{String, Symbol} = "y", cell_col::Union{String, Symbol} = "cell")
if input_type=="fullview"
    transcripts = read_cosmx_transcript(cosmx_dir; fov = nothing)
    cell_coord = DataFrames.combine(groupby(transcripts, [:CellId, :fov]), :x => mean => :x, :y => mean => :y)
    gene_counts = make_ct_from_tx(transcripts)
    count_mtx = make_ct_from_df(gene_counts)
    rename!(transcripts, :target => :gene, :CellId => :cell)
    rename!(cell_coord, :CellId => :cell)
    fov_img = read_cosmx_image(cosmx_dir; fov=nothing)
    fov_positions = CSV.read(cosmx_dir * "/RunSummary/latest.fovs.csv", DataFrame; header = false)
    col_num = length(StatsBase.countmap(fov_positions.Column2))
    row_num = collect(values(StatsBase.countmap(fov_positions.Column2)))[1]
    grid_data = DataFrame(:x=>generate_repeating_pattern(row_num, col_num), :y=>generate_alternating_pattern(row_num, col_num))
    cosmx_obj = CosMxObject(transcripts, cell_coord, count_mtx; prefix = prefix, postfix=postfix,
                            meta_data=meta_data, min_gene=min_gene, min_cell=min_cell,x_col=x_col,
                            y_col=y_col, cell_col=cell_col)
    img_obj = SpaImageObj(nothing, nothing, fov_img)
    cosmx_obj.imgObject = img_obj
    cosmx_obj.gridData = grid_data
elseif input_type == "singleview"
    if isa(fov, Nothing)
        error("Please select a fov to analyze!")
    end
    transcripts = read_cosmx_transcript(cosmx_dir; fov = fov)
    cell_coord = DataFrames.combine(groupby(transcripts, [:CellId, :fov]), :x => mean => :x, :y => mean => :y)
    gene_counts = make_ct_from_tx(transcripts)
    count_mtx = make_ct_from_df(gene_counts)
    rename!(transcripts, :target => :gene, :CellId => :cell)
    rename!(cell_coord, :CellId => :cell)
    fov_img = read_cosmx_image(cosmx_dir; fov=fov)
    cosmx_obj = CosMxObject(transcripts, cell_coord, count_mtx; prefix = prefix, postfix=postfix,
                            meta_data=meta_data, min_gene=min_gene, min_cell=min_cell,x_col=x_col,
                            y_col=y_col, cell_col=cell_col)
    img_obj = SpaImageObj(nothing, nothing, fov_img)
    cosmx_obj.imgObject = img_obj
else
    error("input_type can only be \"fullview\" or \"singleview\"!")
end

return cosmx_obj
end

function read_visiumHD(hd_dir::String; 
    min_genes::Union{Vector{Int64}, Tuple{Int64} }= [0, 0, 0], 
    min_cells::Union{Vector{Int64}, Tuple{Int64} }= [0, 0, 0],
    prefix::Union{String, Nothing} = nothing,
    postfix::Union{String, Nothing} = nothing,
    default_bin = "8_um"
)
    layers = Layers()
    println("1. loading 2um binned data...")
    tenx_dir = hd_dir * "/binned_outputs/square_002um/filtered_feature_bc_matrix"
    pos_file = hd_dir * "/binned_outputs/square_002um/spatial/tissue_positions.parquet"
    json_file = hd_dir * "/binned_outputs/square_002um/spatial/scalefactors_json.json"
    counts = read_10x(tenx_dir; version ="v3", min_gene=min_genes[1], min_cell=min_cells[1])
    layer1 = Layer(counts; prefix = prefix, postfix = postfix)
    all_cells = layer1.rawCount.cell_name
    pos1 = read_hd_pos(pos_file)
    pos1 = filter(:barcode => ∈(Set(all_cells)), pos1)
    pos1 = reorder(pos1, "barcode", all_cells)
    rename!(pos1, :barcode => :cell)
    layer1.spmetaData = pos1
    json1 = JSON.parsefile(json_file)
    layer1.jsonParameters = json1
    layers.layers["2_um"] = layer1
    println("1. 2um binned data loaded!")
    println("2. loading 8um binned data...")
    tenx_dir = hd_dir * "/binned_outputs/square_008um/filtered_feature_bc_matrix"
    pos_file = hd_dir * "/binned_outputs/square_008um/spatial/tissue_positions.parquet"
    json_file = hd_dir * "/binned_outputs/square_008um/spatial/scalefactors_json.json"
    cluster_file = hd_dir * "/binned_outputs/square_008um/analysis/clustering/gene_expression_graphclust/clusters.csv"
    umap_file = hd_dir * "/binned_outputs/square_008um/analysis/umap/gene_expression_2_components/projection.csv"
    counts = read_10x(tenx_dir; version ="v3", min_gene=min_genes[2], min_cell=min_cells[2])
    layer2 = Layer(counts; prefix = prefix, postfix = postfix)
    all_cells = layer2.rawCount.cell_name
    pos2 = read_hd_pos(pos_file)
    pos2 = filter(:barcode => ∈(Set(all_cells)), pos2)
    pos2 = reorder(pos2, "barcode", all_cells)
    layer2.spmetaData = pos2
    json2 = JSON.parsefile(json_file)
    layer2.jsonParameters = json2
    hd_umap =  DataFrame(CSV.File(umap_file))
    clustering =  DataFrame(CSV.File(cluster_file))
    rename!(clustering, ["cell", "cluster"]) 
    if isa(prefix, String)
        clustering.cell = prefix .*  "_" .* string.(clustering.cell)
    end
    if isa(postfix, String)
        clustering.cell = string.(clustering.cell) .* "_" .* postfix
    end
    clustering = filter(:cell=> ∈(Set(all_cells)), clustering)
    layer2.spmetaData.cluster = string.(clustering.cluster)
    layer2.metaData.cluster = string.(clustering.cluster)
    if isa(prefix, String)
        hd_umap.Barcode = prefix .*  "_" .* string.(hd_umap.Barcode)
    end
    if isa(postfix, String)
        hd_umap.Barcode = string.(hd_umap.Barcode) .* "_" .* postfix
    end
    hd_umap = filter(:Barcode => ∈(Set(all_cells)), hd_umap)
    hd_umap = Matrix(hd_umap[!, 2:end])
    umap_obj = UMAPObject(hd_umap, "UMAP", 2, nothing, nothing, nothing, nothing, nothing)
    reduct_obj = ReductionObject(nothing, nothing, umap_obj)
    layer2.dimReduction = reduct_obj
    layers.layers["8_um"] = layer2
    println("2. 8um binned data loaded!")
    println("3. loading 16um binned data...")
    tenx_dir = hd_dir * "/binned_outputs/square_016um/filtered_feature_bc_matrix"
    pos_file = hd_dir * "/binned_outputs/square_016um/spatial/tissue_positions.parquet"
    json_file = hd_dir * "/binned_outputs/square_016um/spatial/scalefactors_json.json"
    cluster_file = hd_dir * "/binned_outputs/square_016um/analysis/clustering/gene_expression_graphclust/clusters.csv"
    umap_file = hd_dir * "/binned_outputs/square_016um/analysis/umap/gene_expression_2_components/projection.csv"
    counts = read_10x(tenx_dir; version ="v3", min_gene=min_genes[3], min_cell=min_cells[3])
    layer3 = Layer(counts; prefix = prefix, postfix = postfix)
    all_cells = layer3.rawCount.cell_name
    pos3 = read_hd_pos(pos_file)
    pos3 = filter(:barcode => ∈(Set(all_cells)), pos3)
    pos3 = reorder(pos3, "barcode", all_cells)
    layer3.spmetaData = pos3
    json3 = JSON.parsefile(json_file)
    layer3.jsonParameters = json3
    hd_umap =  DataFrame(CSV.File(umap_file))
    clustering =  DataFrame(CSV.File(cluster_file))
    rename!(clustering, ["cell", "cluster"]) 
    if isa(prefix, String)
        clustering.cell = prefix .*  "_" .* string.(clustering.cell)
    end
    if isa(postfix, String)
        clustering.cell = string.(clustering.cell) .* "_" .* postfix
    end
    clustering = filter(:cell=> ∈(Set(all_cells)), clustering)
    layer3.spmetaData.cluster = string.(clustering.cluster)
    layer3.metaData.cluster = string.(clustering.cluster)
    if isa(prefix, String)
        hd_umap.Barcode = prefix .*  "_" .* string.(hd_umap.Barcode)
    end
    if isa(postfix, String)
        hd_umap.Barcode = string.(hd_umap.Barcode) .* "_" .* postfix
    end
    hd_umap = filter(:Barcode => ∈(Set(all_cells)), hd_umap)
    hd_umap = Matrix(hd_umap[!, 2:end])
    umap_obj = UMAPObject(hd_umap, "UMAP", 2, nothing, nothing, nothing, nothing, nothing)
    reduct_obj = ReductionObject(nothing, nothing, umap_obj)
    layer3.dimReduction = reduct_obj
    layers.layers["16_um"] = layer3
    println("3. 16um binned data loaded!")
    highres_image_file = hd_dir * "/spatial/tissue_hires_image.png"
    lowres_image_file = hd_dir * "/spatial/tissue_lowres_image.png"
    fullres_image_file = hd_dir * "/spatial/tissue_fullres_image.png"
    detected_tissue = hd_dir * "/spatial/detected_tissue_image.jpg"
    aligned_image_file = hd_dir * "/spatial/aligned_fiducials.jpg"
    high_img = FileIO.load(highres_image_file)
    low_img = FileIO.load(lowres_image_file)
    if isfile(fullres_image_file)
        full_img = FileIO.load(fullres_image_file)
        full_img = convert(Matrix{RGB{N0f8}}, full_img)
    else
        full_img = nothing
    end
    tissue_img = FileIO.load(detected_tissue)
    aligned_img = FileIO.load(aligned_image_file)
    hd_obj = VisiumHDObject(layers; default_bin=default_bin)
    image_obj = VisiumImgObject(high_img, low_img, full_img, tissue_img, aligned_img, layers.layers[default_bin].jsonParameters)
    hd_obj.imageData = image_obj
    return hd_obj
end