function run_tf_idf(mtx::AbstractMatrix{<:Real}; scale_factor = 10000)
    ### compute TF
    npeaks = colSum(mtx)
    tf = *(mtx, Diagonal(1 ./ vec(npeaks))')

    ### compute IDF
    rsums = rowSum(mtx)
    idf = size(mtx)[2] ./ rsums

    ### normalize count using TF-IDF (equivalent to Signac RunTFIDF with parameter method=1)
    norm_data = *(Diagonal(vec(idf)), tf)
    norm_data = log1p.(norm_data .* scale_factor)
end

function run_tf_idf(ct_obj::RawCountObject; scale_factor = 10000)
    norm_count = run_tf_idf(ct_obj.count_mtx, scale_factor=scale_factor)
    norm_method = "TF-IDF"
    pseudocount = 1
    norm_obj = NormCountObject(norm_count, ct_obj.cell_name, ct_obj.gene_name, scale_factor, norm_method, pseudocount)
end

function run_tf_idf(atac_obj::scATACObject; scale_factor = 10000)
    norm_obj = run_tf_idf(atac_obj.rawCount; scale_factor = scale_factor)
    atac_obj.normCount = norm_obj
    return atac_obj
end

function find_top_features(ct_mtx::NormCountObject; min_cutoff="q5")
    counts = ct_mtx.count_mtx
    featurecounts = rowSum(counts)
    e_dist = StatsBase.ecdf(vec(featurecounts))
    peak_names = ct_mtx.gene_name
    hvf_info = DataFrame(
        peak_name = peak_names,
        count = vec(featurecounts),
        percentile = vec(e_dist.(featurecounts))
    )
    sort!(hvf_info, :count, rev=true)
    perc = parse(Int64, filter(x->'0'<=x<='9', min_cutoff))/100
    hvf_info = filter(:percentile => >=(perc), hvf_info)
    return hvf_info
end

function find_top_features(atac_obj::scATACObject; min_cutoff="q5")
    hvf_info = find_top_features(atac_obj.normCount; min_cutoff=min_cutoff)
    Features = hvf_info.peak_name
    var_obj = VariableGeneObject(Features, hvf_info)
    atac_obj.varPeak = var_obj
    return atac_obj
end

function run_svd(atac_obj::scATACObject;  method=:svd, pratio = 1, maxoutdim = 10)
    features = atac_obj.varPeak.var_gene
    hvf_info = atac_obj.varPeak.vst_data
    norm_data = atac_obj.normCount
    peak_names = norm_data.gene_name
    norm_data = norm_data.count_mtx
    norm_data = norm_data[indexin(hvf_info.peak_name, peak_names),:]
    pca_mat = Matrix(norm_data')
    pca_mat = convert(Matrix{Float64}, pca_mat)
    M = MultivariateStats.fit(PCA, pca_mat; method=method, pratio=pratio, maxoutdim=maxoutdim)
    proj = MultivariateStats.projection(M)
    percent_var = principalvars(M) ./ tvar(M) * 100
    key = "PC"
    pca_obj = PCAObject(proj, M, percent_var, key , method, pratio, maxoutdim)
    reduction_obj = ReductionObject(pca_obj, nothing, nothing)
    atac_obj.dimReduction = reduction_obj
    return atac_obj
end

function compute_gene_activity(atac_obj; normalize = true)
    counts = atac_obj.rawCount.count_mtx
    cells = atac_obj.rawCount.cell_name
    peak_names = atac_obj.rawCount.gene_name
    peak_anno = atac_obj.peakData
    peak_anno.peak_names = string.(peak_anno.chrom) .* "_" .* string.(peak_anno.start) .* "_" .* string.(peak_anno.stop)
    counts = DataFrame(Matrix(counts), :auto)
    rename!(counts, cells)
    peak_names = atac_obj.rawCount.gene_name
    counts.peaks = peak_names
    counts = DataFrames.stack(counts, 1:(size(counts)[2]-1))
    counts = map_values(counts, :peaks, :gene, peak_anno.peak_names, peak_anno.gene)
    counts = counts[!, 2:end]
    gdf = DataFrames.combine(DataFrames.groupby(counts, [:variable, :gene]), :value => sum)
    gdf2 = unstack(gdf, :gene, :variable, :value_sum)
    count_mat = Matrix(gdf2[!, 2:end])
    count_mat = convert(Matrix{Int64}, count_mat)
    if normalize
        count_mat = normalize_object(count_mat)
    end
    genes = gdf2.gene
    cells = names(gdf2)[2:end]
    activity_obj = GeneActivityObject(peak_anno, count_mat, cells, genes)
    atac_obj.activityData = activity_obj
    return atac_obj
end
