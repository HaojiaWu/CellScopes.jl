function cluster_kmeans(data, K, random_state)
    Random.seed!(random_state)
    @info "Computing centroids with kmeans..."
    result = kmeans(data, K; maxiter=25, tol=1e-4)
    km_centroids = result.centers
    km_labels = result.assignments
    @info "kmeans initialization complete!"
    return km_centroids
end

function init_cluster!(harmony_obj::HarmonyObject, random_state)
    harmony_obj.Y = cluster_kmeans(harmony_obj.Z_cos, harmony_obj.K, random_state)
    column_norms = [norm(harmony_obj.Y[:, k], 2) for k in 1:size(harmony_obj.Y, 2)]
    harmony_obj.Y = harmony_obj.Y ./ column_norms'
    harmony_obj.dist_mat = 2 * (1 .- (harmony_obj.Y' * harmony_obj.Z_cos))
    harmony_obj.R = -harmony_obj.dist_mat
    harmony_obj.R = harmony_obj.R ./ harmony_obj.sigma
    harmony_obj.R .-= maximum(harmony_obj.R, dims=1)
    harmony_obj.R = exp.(harmony_obj.R)
    harmony_obj.R = harmony_obj.R ./ sum(harmony_obj.R, dims=1)
    harmony_obj.E = harmony_obj.Pr_b' .* sum(harmony_obj.R, dims=2)
    harmony_obj.O = harmony_obj.R * transpose(harmony_obj.Phi)
    compute_objective!(harmony_obj)
    push!(harmony_obj.objective_harmony, harmony_obj.objective_kmeans[end])
end

function compute_objective!(harmony_obj::HarmonyObject)
    kmeans_error = sum(harmony_obj.R .* harmony_obj.dist_mat)
    entropy = sum(safe_entropy(harmony_obj.R) .* reshape(harmony_obj.sigma, :, 1))
    x = harmony_obj.R .* reshape(harmony_obj.sigma, :, 1)
    y = transpose(repeat(harmony_obj.theta, 1, harmony_obj.K))
    z = log.((harmony_obj.O .+ 1) ./ (harmony_obj.E .+ 1))
    w = y .* z * harmony_obj.Phi
    cross_entropy = sum(x .* w)
    push!(harmony_obj.objective_kmeans, kmeans_error + entropy + cross_entropy)
    push!(harmony_obj.objective_kmeans_dist, kmeans_error)
    push!(harmony_obj.objective_kmeans_entropy, entropy)
    push!(harmony_obj.objective_kmeans_cross, cross_entropy)
end

function safe_entropy(x::Array{Float64, 2})
    y = x .* log.(x)
    y .= ifelse.(isfinite.(y), y, 0.0)
    return y
end

function harmonize!(harmony_obj::HarmonyObject, iter_harmony::Int, verbose::Bool)
    is_converged = true
    for i in 1:iter_harmony
        if verbose
            @info "Iteration $i of $iter_harmony"
        end
        cluster!(harmony_obj)
        harmony_obj.Z_cos, harmony_obj.Z_corr, harmony_obj.W, harmony_obj.Phi_Rk = moe_correct_ridge(
            harmony_obj.Z_orig, harmony_obj.Z_cos, harmony_obj.Z_corr, harmony_obj.R, harmony_obj.W, harmony_obj.K,
            harmony_obj.Phi_Rk, harmony_obj.Phi_moe, harmony_obj.lamb
        )
        converged = check_convergence(harmony_obj, 1)
        if converged
            if verbose
                @info "Converged after $i iteration$(i > 1 ? "s" : "")"
            end
            break
        end
        is_converged = converged
    end
    if verbose && !is_converged
        @info "Stopped before convergence"
    end
end

function cluster!(harmony_obj::HarmonyObject)
    harmony_obj.dist_mat = 2 * (1 .- transpose(harmony_obj.Y) * harmony_obj.Z_cos)
    j = 0
    for i in 1:harmony_obj.max_iter_kmeans
        harmony_obj.Y = harmony_obj.Z_cos * transpose(harmony_obj.R)
        column_norms = [norm(harmony_obj.Y[:, k], 2) for k in 1:size(harmony_obj.Y, 2)]
        harmony_obj.Y = harmony_obj.Y ./ column_norms'
        harmony_obj.dist_mat = 2 * (1 .- transpose(harmony_obj.Y) * harmony_obj.Z_cos)
        update_R!(harmony_obj)
        compute_objective!(harmony_obj)
        if i > harmony_obj.window_size
            converged = check_convergence(harmony_obj, 0)
            if converged
                break
            end 
        end
        j = i
    end
    push!(harmony_obj.kmeans_rounds, j)
    push!(harmony_obj.objective_harmony, harmony_obj.objective_kmeans[end])
end

function update_R!(harmony_obj::HarmonyObject)
    harmony_obj.scale_dist = -harmony_obj.dist_mat
    harmony_obj.scale_dist = harmony_obj.scale_dist ./ reshape(harmony_obj.sigma, :, 1)
    harmony_obj.scale_dist .-= maximum(harmony_obj.scale_dist, dims=1)
    harmony_obj.scale_dist = exp.(harmony_obj.scale_dist)
    update_order = randperm(harmony_obj.N)
    n_blocks = ceil(Int, 1 / harmony_obj.block_size)
    blocks = [update_order[1 + floor(Int, harmony_obj.N / n_blocks) * (i - 1):min(floor(Int, harmony_obj.N / n_blocks) * i, harmony_obj.N)] for i in 1:n_blocks]
    for b in blocks
        harmony_obj.E .-= sum(harmony_obj.R[:, b], dims=2) * harmony_obj.Pr_b'
        harmony_obj.O .-= harmony_obj.R[:, b] * transpose(harmony_obj.Phi[:, b])
        harmony_obj.R[:, b] = harmony_obj.scale_dist[:, b]
        theta_expanded = reshape(harmony_obj.theta, 1, :)
        E_div_O = (harmony_obj.E .+ 1) ./ (harmony_obj.O .+ 1)
        power_result = E_div_O .^ theta_expanded
        dot_result = power_result * harmony_obj.Phi[:, b]
        harmony_obj.R[:, b] = harmony_obj.R[:, b] .* dot_result
        column_norms = [norm(harmony_obj.R[:, b][:, i], 1) for i in 1:size(harmony_obj.R[:, b], 2)]
        harmony_obj.R[:, b] = harmony_obj.R[:, b] ./ column_norms'
        harmony_obj.E .+= sum(harmony_obj.R[:, b], dims=2) * harmony_obj.Pr_b'
        harmony_obj.O .+= harmony_obj.R[:, b] * transpose(harmony_obj.Phi[:, b])
    end
end

function check_convergence(harmony_obj::HarmonyObject, i_type::Int)
    obj_old = 0.0
    obj_new = 0.0
    if i_type == 0
        okl = length(harmony_obj.objective_kmeans)
        for i in 1:(harmony_obj.window_size-1)
            obj_old += harmony_obj.objective_kmeans[okl - 2 - i]
            obj_new += harmony_obj.objective_kmeans[okl - 1 - i]
        end
        if abs(obj_old - obj_new) / abs(obj_old) < harmony_obj.epsilon_kmeans
            return true
        else
            return false
        end
    end
    if i_type == 1
        obj_old = harmony_obj.objective_harmony[end - 1]
        obj_new = harmony_obj.objective_harmony[end]
        if (obj_old - obj_new) / abs(obj_old) < harmony_obj.epsilon_harmony
            return true
        else
            return false
        end
    end
end

function moe_correct_ridge(Z_orig, Z_cos, Z_corr, R, W, K, Phi_Rk, Phi_moe, lamb)
    Z_corr = deepcopy(Z_orig)
    for i in 1:K
        Phi_Rk = Phi_moe .* R[i, :]'
        x = Phi_Rk * transpose(Phi_moe) + lamb
        W = inv(x) * Phi_Rk * transpose(Z_orig)
        W[1, :] .= 0.0
        Z_corr .-= transpose(W) * Phi_Rk
    end
    Z_cos = Z_corr ./ [norm(Z_cos[:, i], 2) for i in 1:size(Z_cos, 1)]
    return Z_cos, Z_corr, W, Phi_Rk
end

function get_dummies(df::DataFrame, cols::Vector{Symbol})
    result_df = select(df, Not(cols))
    for col in cols
        cat_col = categorical(df[:, col], ordered=false)
        for level in levels(cat_col)
            result_df[!, Symbol("$(col)_", level)] = Int.(cat_col .== level)
        end
    end
    return result_df
end

