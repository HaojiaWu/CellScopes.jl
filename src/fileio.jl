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
        genes = DataFrame(load(File(format"TSV", gene_file); header_exists=false))
        cells = DataFrame(load(File(format"TSV", cell_file); header_exists=false))
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
    counts = counts[gene_kept, cell_kept]
    gene_kept, gene_removed = check_duplicates(genes)
    gene_removed = collect(values(gene_removed))
    counts = counts[Not(gene_removed), :]
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
    counts = counts[gene_kept, cell_kept]
    gene_kept, gene_removed = check_duplicates(genes)
    gene_removed = collect(values(gene_removed))
    counts = counts[Not(gene_removed), :]
    rawcount = RawCountObject(counts, cells, gene_kept)
    # prepare spatial info
    positions = DataFrame(CSV.File(position_file, header=false))
    rename!(positions, ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"])
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
    gene_rm = [blank_gene; neg_gene]
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
    count_cells = filter(:cell => ∈(Set(clustering.Barcode)), count_cells)
    count_cells.cluster = clustering.Cluster
    spObj = XeniumObject(count_molecules, count_cells, raw_count, poly, umap_obj;
            prefix = prefix, min_gene = min_gene, min_cell = min_gene)
    clustering.cell = prefix .*  "_" .* string.(clustering.Barcode)
    all_cells = spObj.rawCount.cell_name
    clustering = filter(:cell=> ∈(Set(all_cells)), clustering)
    spObj.spmetaData.cell.cluster = string.(clustering.Cluster)
    spObj.metaData.cluster = string.(clustering.Cluster)
    spObj.spmetaData.molecule = map_values(spObj.spmetaData.molecule, :cell, :cluster, clustering.cell, clustering.Cluster)
    spObj.spmetaData.molecule.cluster = string.(spObj.spmetaData.molecule.cluster)
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
    indices = h5read(h5file_path, "matrix/indices") .+= 1
    indptr = h5read(h5file_path, "matrix/indptr") .+= 1
    shape = h5read(h5file_path, "matrix/shape")
    genes = h5read(h5file_path, "matrix/features/name")
    cells = h5read(h5file_path, "matrix/barcodes")
    count_mtx = SparseMatrixCSC(shape[1], shape[2], indptr, indices, data)
    raw_ct = RawCountObject(count_mtx, cells, genes)
    return raw_ct
end

