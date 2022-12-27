function normalize_count(mtx::Matrix{Float64}; scale_factor = 1000000, norm_method = "logarithm", pseudocount = 1)
    norm_count = log.((mtx ./ sum(mtx, dims=1)) .* scale_factor .+ pseudocount)
    return norm_count
end

function NormalizeCount(ct_obj::RawCountObject; scale_factor = 1000000, norm_method = "logarithm", pseudocount = 1)
    norm_count = normalize_count(ct_obj.count_mtx; scale_factor=scale_factor, norm_method=norm_method, pseudocount=pseudocount)
    norm_obj = NormCountObject(norm_count, ct_obj.cell_name, ct_obj.gene_name, scale_factor, norm_method, pseudocount)
    return norm_obj
end

function FindVariableGenes(ct_mtx::RawCountObject; nFeatures::Int64 = 2000, span::Float64 = 0.3)
    mean_val = mean(ct_mtx.count_mtx, dims=2)
    var_val = var(ct_mtx.count_mtx, dims=2)
    vst_data = [mean_val var_val zeros(length(mean_val)) zeros(length(mean_val))]
    vst_data = DataFrame(vst_data, :auto)
    rename!(vst_data, ["mean", "variance", "variance_expected","variance_standardized"])
    vst_data = filter(:variance => x -> x > 0.0, vst_data)
    fit = loess(log10.(vst_data.mean), log10.(vst_data.variance), span=span)
    vst_data.variance_expected = 10 .^ predict(fit, log10.(vst_data.mean))
    mat = deepcopy(ct_mtx.count_mtx);
    mat2 = (mat .- vst_data.mean) ./ sqrt.(vst_data.variance_expected)
    sd_val = var(mat2, dims=2)
    vst_data.variance_standardized = vec(sd_val)
    vst_data.gene = ct_mtx.gene_name;
    vst_data = sort(vst_data, :variance_standardized, rev=true)
    Features = vst_data.gene[1:nFeatures]
    return vst_data, Features
end
