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

function load(filename::String)
    cs = JLD2.load(filename)
    cs = cs["key"]
    return cs
end

function save_cs(filename::String, obj)
    open(filename, "w") do io
        serialize(io, obj)
    end
end

function load_cs(filename::String)
    open(filename, "r") do io
        return deserialize(io)
    end
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

function read_xenium(xenium_dir::String; prefix = nothing, min_gene::Int64 = 0, min_cell::Int64 = 0, version="1.1")
    gene_file = xenium_dir * "/cell_feature_matrix/features.tsv.gz"
    cell_file = xenium_dir * "/cell_feature_matrix/barcodes.tsv.gz"
    count_file = xenium_dir * "/cell_feature_matrix/matrix.mtx.gz"
    cell_meta = xenium_dir * "/cells.csv.gz"
    transcript_meta = xenium_dir * "/transcripts.parquet"
    seg_file = xenium_dir * "/cell_boundaries.parquet"
    cluster_file = xenium_dir * "/analysis/clustering/gene_expression_graphclust/clusters.csv"
    umap_file = xenium_dir * "analysis/umap/gene_expression_2_components/projection.csv"
    xenium_umap =  DataFrame(CSV.File(umap_file))
    genes = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", gene_file); header_exists=false))
    cells = DataFrame(CSVFiles.load(CSVFiles.File(format"TSV", cell_file); header_exists=false))
    clustering =  DataFrame(CSV.File(cluster_file))
    cell_ids_filter = clustering.Barcode
    seg = read_parquet(seg_file)
    seg.cell_id = String.(seg.cell_id)
    cell_ids = unique(seg.cell_id)
    cell_kept = check_vec(cell_ids_filter, cell_ids)
    count_molecules = read_parquet(transcript_meta)
    count_cells =  DataFrame(CSV.File(cell_meta))
    counts = MatrixMarket.mmread(gunzip(count_file))
    counts = counts[:, cell_kept]
    seg = filter(:cell_id => ∈(Set(clustering.Barcode)), seg)
    cells = filter(:Column1 => ∈(Set(clustering.Barcode)), cells)
    count_molecules = filter(:cell_id => ∈(Set(clustering.Barcode)), count_molecules)
    count_cells = filter(:cell_id => ∈(Set(clustering.Barcode)), count_cells)
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
    println("\033[1;34mFiltering cells and genes...\033[0m")
    raw_count = subset_count(raw_count; genes=genes2)
    println("Cells and genes filtered!")
    count_molecules.cell = string.(count_molecules.cell)
    count_molecules = filter(:gene => ∈(Set(genes2)), count_molecules)
    count_cells.cell = string.(count_cells.cell)
    count_cells = filter(:cell => ∈(Set(string.(clustering.Barcode))), count_cells)
    count_cells.cluster = clustering.Cluster
    spObj = XeniumObject(count_molecules, count_cells, raw_count;
            prefix = prefix, min_gene = min_gene, min_cell = min_gene)
    clustering.cell = string.(clustering.Barcode)
    if !isa(prefix, Nothing)
        clustering.cell = prefix .*  "_" .* string.(clustering.cell)
        seg.cell_id = prefix .*  "_" .* string.(seg.cell_id)
    end
    all_cells = spObj.rawCount.cell_name
    clustering = filter(:cell=> ∈(Set(all_cells)), clustering)
    spObj.spmetaData.cell.cluster = string.(clustering.Cluster)
    spObj.metaData.cluster = string.(clustering.Cluster)
    spObj.spmetaData.molecule = map_values(spObj.spmetaData.molecule, :cell, :cluster, clustering.cell, clustering.Cluster)
    spObj.spmetaData.molecule.cluster = string.(spObj.spmetaData.molecule.cluster)

    seg = filter(:cell_id=> ∈(Set(all_cells)), seg)
    grouped = groupby(seg, :cell_id)
    cell_ids = unique(seg.cell_id)
    poly = Vector{Matrix{Float64}}(undef, length(cell_ids))
    n = length(cell_ids)
    println("\033[1;34mFormatting cell polygons...\033[0m")
    p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:blue)
    for idx in 1:length(cell_ids)
        cell_data = grouped[idx]
        cell_1 = Matrix(cell_data[!, 2:end])
        poly[idx] = cell_1
        next!(p)
    end
    println("Cell polygons formatted!")
    polygon_df = DataFrame(polygon_number = 1:length(poly), mapped_cell = all_cells, cluster=spObj.spmetaData.cell.cluster)
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
patterns = [r"^NegControlProbe_", r"^Unassigned", r"^antisense_", r"^NegControlCodeword_", r"^BLANK_"]
if is_match(String.(count_df.gene), patterns)
    count_df = @chain count_df begin
        @rsubset !match_patterns(:gene, patterns)
    end

    molecules = @chain molecules begin
        @rsubset !match_patterns(:gene, patterns)
    end
end
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
spObj.polygonData = poly
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

function read_cellseg_geojson(geojson_file)
    geojson = JSON.parsefile(geojson_file)
    features = geojson["features"]
    rows = []
    for (i, f) in enumerate(features)
        geom = f["geometry"]
        gtype = geom["type"]
        props = f["properties"]

        cell_int = haskey(props, "cell_id") ? props["cell_id"] : i
        cell_id = "cellid_" * lpad(string(cell_int), 9, '0') * "-1"
        cluster = get(get(props, "classification", Dict()), "name", missing)

        if gtype == "Polygon"
            for ring in geom["coordinates"]
                for pt in ring
                    push!(rows, (cell_id = cell_id, x = pt[2], y = pt[1], cluster = cluster))
                end
            end
        elseif gtype == "MultiPolygon"
            for poly in geom["coordinates"]
                for ring in poly
                    for pt in ring
                        push!(rows, (cell_id = cell_id, x = pt[2], y = pt[1], cluster = cluster))
                    end
                end
            end
        elseif gtype == "Point"
            pt = geom["coordinates"]
            push!(rows, (cell_id = cell_id, x = pt[2], y = pt[1], cluster = cluster))
        else
            @warn "Unhandled geometry type: $gtype"
        end
    end

    df = DataFrame(rows)
    return df
end

function read_layers(hd_dir; 
    min_gene::Int64 = 0,
    min_cell::Int64 = 0,
    seg_type::String = "bin",
    bin_size::Int64 = 8,
    prefix::Union{String, Nothing}=nothing, 
    postfix::Union{String, Nothing}=nothing
    )
    json_file = hd_dir * "/spatial/scalefactors_json.json"
    json = JSON.parsefile(json_file)
    if seg_type == "bin"
        tenx_dir = hd_dir * "/filtered_feature_bc_matrix"
        pos_file = hd_dir * "/spatial/tissue_positions.parquet"
        pos = read_hd_pos(pos_file)
        px_width = bin_size/json["microns_per_pixel"]
        corner_coordinates = compute_corner_points(pos, px_width; cell = "barcode", x_col = "pxl_row_in_fullres", y_col = "pxl_col_in_fullres")
        seg = DataFrame(x = corner_coordinates.new_x, y = corner_coordinates.new_y, cell_id = corner_coordinates.barcode)
    elseif seg_type == "cell"
        tenx_dir = hd_dir * "/filtered_feature_cell_matrix"
        cellseg_geojson = hd_dir * "/graphclust_annotated_cell_segmentations.geojson"
        seg = read_cellseg_geojson(cellseg_geojson)
        seg = seg[:, [2:ncol(seg); 1]]
        pos = DataFrames.combine(groupby(seg, :cell_id)) do subdf
                (
                    pxl_row_in_fullres = mean(subdf.x),
                    pxl_col_in_fullres = mean(subdf.y)
                )
            end
        rename!(pos, :cell_id => :barcode)
        pos.ID = collect(1:nrow(pos))
    else
        error(""""seg_type" can only be "bin" or "cell"!""")
    end
    counts = read_10x(tenx_dir; version ="v3", min_gene=min_gene, min_cell=min_cell)
    all_cells = counts.cell_name
    pos = filter(:barcode => ∈(Set(all_cells)), pos)
    pos =  pos[(pos[!, :pxl_row_in_fullres] .> 0) .& (pos[!, :pxl_col_in_fullres] .> 0), :]
    if isa(prefix, String)
        pos.barcode = prefix .*  "_" .* string.(pos.barcode)
        counts.cell_name = prefix .*  "_" .* string.(counts.cell_name)
    end
    if isa(postfix, String)
        pos.barcode = string.(pos.barcode ) .* "_" .* postfix
        counts.cell_name = string.(counts.cell_name ) .* "_" .* postfix
    end
    all_cells = pos[!, :barcode]
    counts = subset_count(counts; cells = all_cells)
    layer = Layer(counts)
    all_cells = layer.rawCount.cell_name
    pos = filter(:barcode => ∈(Set(all_cells)), pos)
    pos = reorder(pos, "barcode", all_cells)
    layer.spmetaData = pos
    layer.jsonParameters = json
    if bin_size !== 2
        cluster_file = hd_dir * "/analysis/clustering/gene_expression_graphclust/clusters.csv"
        umap_file = hd_dir * "/analysis/umap/gene_expression_2_components/projection.csv"
        hd_umap =  DataFrame(CSV.File(umap_file))
        clustering =  DataFrame(CSV.File(cluster_file))
        rename!(clustering, ["cell", "cluster"]) 
        if isa(prefix, String)
            clustering.cell = prefix .*  "_" .* string.(clustering.cell)
            hd_umap.Barcode = prefix .*  "_" .* string.(hd_umap.Barcode)
        end
        if isa(postfix, String)
            clustering.cell = string.(clustering.cell) .* "_" .* postfix
            hd_umap.Barcode = string.(hd_umap.Barcode) .* "_" .* postfix
        end
        clustering = filter(:cell=> ∈(Set(all_cells)), clustering)
        layer.spmetaData.cluster = string.(clustering.cluster)
        layer.metaData.cluster = string.(clustering.cluster)
        hd_umap = filter(:Barcode => ∈(Set(all_cells)), hd_umap)
        hd_umap = Matrix(hd_umap[!, 2:end])
        umap_obj = UMAPObject(hd_umap, "UMAP", 2, nothing, nothing, nothing, nothing, nothing)
        reduct_obj = ReductionObject(nothing, nothing, umap_obj)
        layer.dimReduction = reduct_obj
    end
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
    poly_data = Polygons()
    poly_data.polygons["original"] = poly
    layer.polygonData = poly_data
    return layer
end

function read_visiumHD(hd_dir::String; 
    min_genes::Union{Nothing, Vector{Int}, Tuple{Int}} = nothing, 
    min_cells::Union{Nothing, Vector{Int}, Tuple{Int}} = nothing,
    prefix::Union{String, Nothing} = nothing,
    postfix::Union{String, Nothing} = nothing,
    default_bin = "8_um"
)
    if !isdir(hd_dir * "/segmented_outputs")
        min_genes === nothing && (min_genes = [0, 0, 0])
        min_cells === nothing && (min_cells = [0, 0, 0])
    else
        min_genes === nothing && (min_genes = [0, 0, 0, 0])
        min_cells === nothing && (min_cells = [0, 0, 0, 0])
    end
    layers = Layers()
    println("\033[1;34m1. loading 2um binned data...\033[0m")
    bin_dir = hd_dir * "/binned_outputs/square_002um"
    layer1 = read_layers(bin_dir; min_gene = min_genes[1], min_cell =  min_cells[1], prefix = prefix, postfix = postfix, seg_type = "bin", bin_size = 2)
    layers.layers["2_um"] = layer1
    println("1. 2um binned data loaded!")
    println("\033[1;34m2. loading 8um binned data...\033[0m")
    bin_dir = hd_dir * "/binned_outputs/square_008um"
    layer2 = read_layers(bin_dir; min_gene = min_genes[2], min_cell =  min_cells[2], prefix = prefix, postfix = postfix, seg_type = "bin", bin_size = 8)
    layers.layers["8_um"] = layer2
    println("2. 8um binned data loaded!")
    println("\033[1;34m3. loading 16um binned data...\033[0m")
    bin_dir = hd_dir * "/binned_outputs/square_016um"
    layer3 = read_layers(bin_dir; min_gene = min_genes[3], min_cell =  min_cells[3], prefix = prefix, postfix = postfix, seg_type = "bin", bin_size = 16)
    layers.layers["16_um"] = layer3
    println("3. 16um binned data loaded!")
    if isdir(hd_dir * "/segmented_outputs")
        println("\033[1;34m4. loading cell seg data...\033[0m")
        cell_dir = hd_dir * "/segmented_outputs"
        layer4 = read_layers(cell_dir; min_gene = min_genes[4], min_cell =  min_cells[4], prefix = prefix, postfix = postfix, seg_type = "cell")
        layers.layers["cell"] = layer4
        println("4. Cell seg data loaded!")
    end
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

function read_paired_data(xn_dir, vs_dir, xn_img_path, vs_img_path; 
    vs_mat::Union{Matrix{Float64}, Nothing} = nothing,
    xn_mat::Union{Matrix{Float64}, Nothing} = nothing,
    kwargs...
)
    @info("Step I. Reading visiumHD data...")
    hd_obj = read_visiumHD(vs_dir)
    @info("Step II. Reading Xenium data...")
    xn_obj = read_xenium(xn_dir)
    @info("Step III. Pairing both data...")
    cell_coord = deepcopy(xn_obj.spmetaData.cell)
    mol_coord = deepcopy(xn_obj.spmetaData.molecule)
    cell_coord.x ./= 0.2125
    cell_coord.y ./= 0.2125
    mol_coord.x ./= 0.2125
    mol_coord.y ./= 0.2125
    vs_img = FileIO.load(vs_img_path)
    vs_img = convert(Matrix{RGB{N0f8}}, vs_img)
    xn_img = FileIO.load(xn_img_path)
    xn_img = convert(Matrix{RGB{N0f8}}, xn_img)
    println("Generating cell count...")
    cell_counts, molecule_data = generate_hd_segcount(xn_dir, vs_dir; t_mat = vs_mat, img_lims=size(vs_img))
    println("Transforming coordinates...")
    if !isa(vs_mat, Nothing)
        inv_vs_mat = inv(vs_mat)
        cell_coord = transform_coord(cell_coord, inv_vs_mat; x_old = :x, y_old = :y, x_new=:y, y_new = :x)
        mol_coord = transform_coord(mol_coord, inv_vs_mat; x_old = :x, y_old = :y, x_new=:y, y_new = :x)
        poly = reformat_polygons(xn_dir, vs_mat)
        xn_obj.polygonData = poly
        polygon_df = DataFrame(polygon_number = 1:length(poly), mapped_cell = cell_coord.cell, cluster=cell_coord.cluster)
        xn_obj.spmetaData.polygon = polygon_df
    end

    if !isa(xn_mat, Nothing) 
        if isa(vs_mat, Nothing)
            xn_mat = xn_mat
        else
            xn_mat = inv(vs_mat) * xn_mat
        end
        h1, w1 = size(xn_img)
        new_df = DataFrame(x = repeat(1:w1, inner = h1),
           y = repeat(1:h1, outer = w1),
           color = vec(xn_img))
        new_df = transform_coord(new_df, xn_mat; x_old = :x, y_old = :y, x_new=:new_y, y_new = :new_x)
        new_df.new_x = Int64.(round.(new_df.new_x))
        new_df.new_y = Int64.(round.(new_df.new_y))
        new_df = new_df[(new_df.new_x .> 0) .&& (new_df.new_y .> 0), :]
        max_y = maximum(new_df.new_y)
        max_x = maximum(new_df.new_x)
        new_img = fill(RGB{N0f8}(1.0, 1.0, 1.0), max_x, max_y)
        indices = CartesianIndex.(new_df.new_x, new_df.new_y)
        new_img[indices] = new_df.color
    end
    new_img = smoothe_img(new_img)
    cell_coord_xn = deepcopy(cell_coord)
    cell_coord_vs = deepcopy(hd_obj.spmetaData)
    coord_xn2,coord_vs2,img_xn2,img_vs2, x_lims, y_lims= crop_img_coord(cell_coord_xn, cell_coord_vs, new_img, vs_img)
    mol_coord[!, "x"] .-= x_lims[1]
    mol_coord[!, "y"] .-= y_lims[1]
    xn_obj.polygonData = [m .- [x_lims[1] y_lims[1]] for m in xn_obj.polygonData]
    xn_obj.spmetaData.cell = coord_xn2
    xn_obj.spmetaData.molecule = mol_coord
    xn_obj.imageData = img_xn2
    hd_obj.spmetaData = coord_vs2
    hd_obj.layerData.layers["8_um"].spmetaData = coord_vs2
    hd_obj.imageData.fullresImage = img_vs2
    paired_sp_obj = PairedSpObject(hd_obj, xn_obj, vs_mat, xn_mat)
    paired_obj = PairedObject(paired_sp_obj, cell_counts; kwargs...)
    cell_kept = cell_counts.cell_name
    cell_data = deepcopy(xn_obj.spmetaData.cell)
    cell_data = filter(:cell => ∈(Set(cell_kept)), cell_data)
    poly_data = deepcopy(xn_obj.spmetaData.polygon)
    poly_data = filter(:mapped_cell => ∈(cell_kept), poly_data)
    poly = xn_obj.polygonData[poly_data.polygon_number]
    poly_data = DataFrame(polygon_number = 1:length(poly), mapped_cell = cell_data.cell, cluster=cell_data.cluster)
    paired_obj.polygonData = poly
    meta = SpaMetaObj(cell_data, molecule_data, poly_data)
    paired_obj.spmetaData = meta
    @info("All done!")
    return paired_obj
end

function read_joint_data(xn_dir, vs_dir, img_path;
    t_mat::Union{Matrix{Float64}, Nothing} = nothing,
    kwargs...
    )
    @info("Step I. Reading visiumHD data...")
    hd_obj = read_visiumHD(vs_dir)
    @info("Step II. Reading Xenium data...")
    xn_obj = read_xenium(xn_dir)
    @info("Step III. Aligning both data...")
    vs_cell = deepcopy(hd_obj.spmetaData)
    x_lims = [minimum(vs_cell.pxl_row_in_fullres), maximum(vs_cell.pxl_row_in_fullres)]
    y_lims = [minimum(vs_cell.pxl_col_in_fullres), maximum(vs_cell.pxl_col_in_fullres)]
    cell_coord = deepcopy(xn_obj.spmetaData.cell)
    mol_coord = deepcopy(xn_obj.spmetaData.molecule)
    cell_coord.x ./= 0.2125
    cell_coord.y./= 0.2125
    mol_coord.x ./= 0.2125
    mol_coord.y ./= 0.2125
    img = FileIO.load(img_path)
    img = convert(Matrix{RGB{N0f8}}, img)
    println("\033[1;34mGenerating cell count...\033[0m")
    cell_counts, molecule_data = generate_hd_segcount(xn_dir, vs_dir; t_mat = t_mat, img_lims=size(img))
    println("\033[1;34mTransforming coordinates...\033[0m")
    inv_mat = inv(t_mat)
    cell_coord = transform_coord(cell_coord, inv_mat; x_old = :x, y_old = :y, x_new=:y, y_new = :x)
    mol_coord = transform_coord(mol_coord, inv_mat; x_old = :x, y_old = :y, x_new=:y, y_new = :x)
    poly = reformat_polygons(xn_dir, t_mat)
    polygon_df = DataFrame(polygon_number = 1:length(poly), mapped_cell = cell_coord.cell, cluster=cell_coord.cluster)
    cell_coord[!, :x] .-= x_lims[1]
    cell_coord[!, :y] .-= y_lims[1]
    mol_coord[!, :x] .-= x_lims[1]
    mol_coord[!, :y] .-= y_lims[1]
    poly = [m .- [x_lims[1] y_lims[1]] for m in poly]
    xn_obj.polygonData = deepcopy(poly)
    xn_obj.spmetaData.polygon = deepcopy(polygon_df)
    xn_obj.spmetaData.cell = deepcopy(cell_coord)
    xn_obj.spmetaData.molecule = deepcopy(mol_coord)
    println("\033[1;34mTrimming images...\033[0m")
    img_df, colors = img_to_df(img)
    img_df.colors = colors
    img_df = img_df[(img_df[!, :x] .> x_lims[1]) .& (img_df[!, :x] .< x_lims[2]) .&
                    (img_df[!, :y] .> y_lims[1]) .& (img_df[!, :y] .< y_lims[2]), :]
    img_new = df_to_img(img_df.x, img_df.y, img_df.colors)
    img_new = smoothe_img(img_new)
    xn_obj.imageData = img_new
    hd_obj.imageData.fullresImage = img_new
    paired_sp_obj = PairedSpObject(hd_obj, xn_obj, nothing, t_mat)
    paired_obj = PairedObject(paired_sp_obj, cell_counts; kwargs...)
    cell_kept = cell_counts.cell_name
    cell_data = deepcopy(xn_obj.spmetaData.cell)
    cell_data = filter(:cell => ∈(Set(cell_kept)), cell_data)
    poly_data = deepcopy(xn_obj.spmetaData.polygon)
    poly_data = filter(:mapped_cell => ∈(cell_kept), poly_data)
    poly = xn_obj.polygonData[poly_data.polygon_number]
    poly_data = DataFrame(polygon_number = 1:length(poly), mapped_cell = cell_data.cell, cluster=cell_data.cluster)
    paired_obj.polygonData = poly
    meta = SpaMetaObj(cell_data, molecule_data, poly_data)
    paired_obj.spmetaData = meta
    paired_obj.pairedData.xnObj = subset_object(paired_obj.pairedData.xnObj; cells=cell_kept)
    @info("All done!")
    return paired_obj
end
