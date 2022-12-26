function normalize_count(mtx::Matrix{Float64}; scale_factor = 1000000, norm_method = "logarithm", pseudocount = 1)
    norm_count = log.((mtx ./ sum(mtx, dims=1)) .* scale_factor .+ pseudocount)
    return norm_count
end

