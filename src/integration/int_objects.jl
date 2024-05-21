abstract type AbstractHarmony <: AbstractCellScope end

mutable struct HarmonyObject <: AbstractHarmony
    Z_corr::Array{Float64, 2}
    Z_orig::Array{Float64, 2}
    Z_cos::Array{Float64, 2}
    Phi::Array{Float64, 2}
    Phi_moe::Array{Float64, 2}
    N::Int
    Pr_b::Vector{Float64}
    B::Int
    d::Int
    window_size::Int
    epsilon_kmeans::Float64
    epsilon_harmony::Float64
    lamb::Array{Float64, 2}
    sigma::Vector{Float64}
    sigma_prior::Vector{Float64}
    block_size::Float64
    K::Int
    max_iter_harmony::Int
    max_iter_kmeans::Int
    verbose::Bool
    theta::Vector{Float64}
    objective_harmony::Vector{Float64}
    objective_kmeans::Vector{Float64}
    objective_kmeans_dist::Vector{Float64}
    objective_kmeans_entropy::Vector{Float64}
    objective_kmeans_cross::Vector{Float64}
    kmeans_rounds::Vector{Int}
    scale_dist::Array{Float64, 2}
    dist_mat::Array{Float64, 2}
    O::Array{Float64, 2}
    E::Array{Float64, 2}
    W::Array{Float64, 2}
    Phi_Rk::Array{Float64, 2}
    Y::Array{Float64, 2}
    R::Array{Float64, 2}

    function HarmonyObject(
        data_mat::Array{Float64, 2},
        meta_data::DataFrame,
        vars_use;
        theta::Union{Nothing, Float64, Int, Array{Float64}} = nothing,
        lamb::Union{Nothing, Float64, Int, Array{Float64}} = nothing,
        sigma::Union{Nothing, Float64} = 0.1,
        nclust::Union{Nothing, Int} = nothing,
        tau::Float64 = 0.0,
        block_size::Float64 = 0.05,
        max_iter_harmony::Int = 10,
        max_iter_kmeans::Int = 20,
        epsilon_cluster::Float64 = 1e-5,
        epsilon_harmony::Float64 = 1e-4,
        plot_convergence::Bool = false,
        verbose::Bool = true,
        reference_values::Nothing = nothing,
        cluster_prior::Nothing = nothing,
        random_state::Int = 0
    )

    N = size(meta_data, 1)
    if size(data_mat, 2) != N
        data_mat = transpose(data_mat)
    end
    @assert size(data_mat, 2) == N "The number of cells in the metadata must match the number of cells in the count matrix."
    if nclust === nothing
        nclust = min(round(Int, N / 30.0), 100)
    end
    K = nclust
    if typeof(sigma) == Float64 && nclust > 1
        sigma = fill(sigma, nclust)
    end       
    phi = get_dummies(metadata, [vars_use])
    phi = phi[:, 2:end]
    phi = transpose(Matrix(phi))
    phi_n = length(unique(meta_data[!, vars_use]))
    if theta === nothing
        theta = fill(1.0, phi_n)
    end
    @assert length(theta) == phi_n 
    if lamb === nothing
        lamb = fill(1.0, phi_n)
    end
    @assert length(lamb) == phi_n
    N_b = sum(phi, dims=2)
    Pr_b = vec(N_b / N)
    if tau > 0
        theta = theta .* (1 .- exp.(-(N_b / (nclust * tau)) .^ 2))
    end
    lamb_mat = diagm(0 => vcat(0.0, lamb))
    phi_moe = vcat(ones(1, N), phi)
    Random.seed!(random_state)                                        
    Z_corr = deepcopy(data_mat)
    Z_orig = deepcopy(data_mat)
    Z_cos = Z_orig' ./ maximum(Z_orig', dims=2)
    Z_cos = Matrix(Z_cos')
    Z_cos = Z_cos ./ mapslices(x -> norm(x, 2), Z_cos, dims=1)   
    N = size(Z_corr, 2)
    B = size(phi, 1)
    d = size(Z_corr, 1)
    window_size = 3
    objective_harmony = Float64[]
    objective_kmeans = Float64[]
    objective_kmeans_dist = Float64[]
    objective_kmeans_entropy = Float64[]
    objective_kmeans_cross = Float64[]
    kmeans_rounds = Int[]
    scale_dist = zeros(K, N)
    dist_mat = zeros(K, N)
    O = zeros(K, B)
    E = zeros(K, B)
    W = zeros(B + 1, d)
    Phi_Rk = zeros(B + 1, N)
    harmony_obj = new(Z_corr, Z_orig, Z_cos, phi, phi_moe, N, Pr_b, B, d, window_size,
               epsilon_cluster, epsilon_harmony, lamb_mat, sigma, sigma, block_size, K,
               max_iter_harmony, max_iter_kmeans, verbose, theta,
               objective_harmony, objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy,
               objective_kmeans_cross, kmeans_rounds, scale_dist, dist_mat, O, E, W, Phi_Rk,
               zeros(Float64, 0, 0), zeros(Float64, 0, 0))
    init_cluster!(harmony_obj, random_state)
    harmonize!(harmony_obj, max_iter_harmony, verbose)
    return harmony_obj
    end
end
