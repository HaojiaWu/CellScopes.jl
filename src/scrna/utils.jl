colSum(mtx::AbstractMatrix{<:Real}) = sum(mtx, dims=1)
rowSum(mtx::AbstractMatrix{<:Real}) = sum(mtx, dims=2)
rownames(sc_obj::get_object_group("All")) = sc_obj.rawCount.gene_name
rownames(sc_obj::scATACObject) = sc_obj.rawCount.peak_name
colnames(sc_obj::get_object_group("All")) = sc_obj.rawCount.cell_name
rownames(ct_mat::AbstractCount) = ct_mat.gene_name
colnames(ct_mat::AbstractCount) = ct_mat.cell_name

function convert_geneid(ensembl_id::String)
    base_url = "https://rest.ensembl.org"
    endpoint = "/lookup/id/"
    format = "?content-type=application/json"
    try
        r = HTTP.get(base_url * endpoint * ensembl_id * format)
        data = JSON.parse(String(r.body))

        if haskey(data, "display_name")
            return data["display_name"]
        else
            return "no_gene"
        end
    catch e
        if isa(e, HTTP.ExceptionRequest.StatusError)
            return "no_gene"
        else
            rethrow(e)
        end
    end
end

function check_duplicates(arr)
    seen_values = Dict{Any, Bool}()
    removed_positions = Dict{Any, Int64}()
    uniq_arr = []
    for (i, value) in enumerate(arr)
        if !haskey(seen_values, value)
            push!(uniq_arr, value)
            seen_values[value] = true
        else
            if !haskey(removed_positions, value)
                removed_positions[value] = i
            end
        end
    end
    return uniq_arr, removed_positions
end

function is_duplicates(list_vec)
    seen = Set()
    for x in list_vec
        if x in seen
            return true
        end
        push!(seen, x)
    end
    return false
end

function subset_count(ct_obj::T; 
    genes::Union{Vector{String}, Nothing} = nothing, 
    cells::Union{Vector{String}, Nothing} = nothing) where T <: AbstractCount
    all_genes = ct_obj.gene_name
    all_cells = ct_obj.cell_name
    if isa(genes, Nothing)
        genes = all_genes
    end
    if isa(cells, Nothing)
        cells = all_cells
    end
    cells_set = Set(cells)
    genes_set = Set(genes)
    check_cell = [i in cells_set for i in all_cells]
    check_gene = [i in genes_set for i in all_genes]
    gene_name = all_genes[check_gene]
    cell_name = all_cells[check_cell]
    new_count = ct_obj.count_mtx[check_gene, check_cell]
    if isa(ct_obj, RawCountObject)
        new_obj = RawCountObject(new_count, cell_name, gene_name)
    elseif isa(ct_obj, NormCountObject)
        new_obj = NormCountObject(new_count, cell_name, gene_name,  ct_obj.scale_factor, ct_obj.norm_method, ct_obj.pseudocount)
    elseif isa(ct_obj, ScaleCountObject)
        new_obj = ScaleCountObject(new_count, cell_name, gene_name,  ct_obj.do_scale, ct_obj.do_center, ct_obj.scale_max)
    else
        new_obj = SpaCountObj(new_count, cell_name, gene_name)
    end
    return new_obj
end

function extract_cluster_count(sc_obj::get_object_group("All"), cl; count_type = "norm", anno = Union{String, Symbol}="cluster")
    df = sc_obj.clustData.clustering
    if isa(anno, String)
        anno = Symbol(anno)
    end
    cl_data = filter(anno => ==(cl), df)
    cl_cell = cl_data.cell_id
    if count_type === "raw"
        ct_obj = sc_obj.rawCount
    elseif count_type === "norm"
        ct_obj = sc_obj.normCount
    elseif count_type === "scale"
        ct_obj = sc_obj.scaleCount
    else 
        println("count_type can only be \"raw\", \"norm\" or \"scale\"!")
    end
    cl_ct = subset_count(ct_obj; cells = cl_cell)
end

function get_dim_data(dim_obj::AbstractDimReduction, dim_type::String)
    if dim_type === "pca"
       dim_data = dim_obj.pca.cell_embedding
       dim_data = DataFrame(dim_data, :auto)
       ndims = dim_obj.pca.maxoutdim
       key = dim_obj.pca.key
       rename!(dim_data, key .* "_" .* string.(collect((1:ndims))))
       x_col = key * "_1"
       y_col =  key * "_2"
       return dim_data, x_col, y_col
    elseif dim_type === "tsne"
       dim_data = dim_obj.tsne.cell_embedding
       dim_data = DataFrame(dim_data, :auto)
       ndims = dim_obj.tsne.ndim
       key = dim_obj.tsne.key
       rename!(dim_data, key .* "_" .* string.(collect((1:ndims))))
       x_col = key * "_1"
       y_col =  key * "_2"
       return dim_data, x_col, y_col
    elseif dim_type === "umap"
       dim_data = dim_obj.umap.cell_embedding
       dim_data = DataFrame(dim_data, :auto)
       ndims = dim_obj.umap.n_components
       key = dim_obj.umap.key
       rename!(dim_data, key .* "_" .* string.(collect((1:ndims))))
       x_col = key * "_1"
       y_col =  key * "_2"
       return dim_data, x_col, y_col
    else
        error("dim_type can only be \"pca\", \"tsne\", or \"umap\"!")
    end    
end

function jitter(x)
    z = abs(-(extrema(skipmissing(x))...))
    a = z/50
    if a == 0
        x = x .+ rand(size(x,1))
        return reduce(vcat, x)
    else
        x = x .+ rand(Uniform(-a, a),size(x,1))
        return reduce(vcat, x)
    end
end

function variable_genes(sc_obj::get_object_group("All"))
    vargenes = pbmc.varGene.var_gene
    return vargenes
end

function update_object(sp_obj::get_object_group("All"))
    cells = colnames(sp_obj)
    genes = rownames(sp_obj)
    all_cells = sp_obj.metaData.Cell_id
    cells_set = Set(cells)
    check_cell = [i in cells_set for i in all_cells]
    println("Updating RNA counts...")
    sp_obj = update_count(sp_obj, :normCount) 
    sp_obj = update_count(sp_obj, :scaleCount)
    cell_set = Set(cells)
    gene_set = Set(genes)
    sp_obj.metaData = filter(:Cell_id => ∈(cell_set), sp_obj.metaData)
    if isdefined(sp_obj, :dimReduction)
        println("Updating dimension reduction...")
        if isdefined(sp_obj.dimReduction, :pca)
            if getfield(sp_obj.dimReduction, :pca) !== nothing
                 sp_obj.dimReduction.pca.cell_embedding = sp_obj.dimReduction.pca.cell_embedding[check_cell, :]
            end
        end
        if isdefined(sp_obj.dimReduction, :tsne)
            if getfield(sp_obj.dimReduction, :tsne) !== nothing
                 sp_obj.dimReduction.tsne.cell_embedding = sp_obj.dimReduction.tsne.cell_embedding[check_cell, :]
            end
        end
        if isdefined(sp_obj.dimReduction, :umap)
            if getfield(sp_obj.dimReduction, :umap) !== nothing
                 sp_obj.dimReduction.umap.cell_embedding = sp_obj.dimReduction.umap.cell_embedding[check_cell, :]
            end
        end
        if isdefined(sp_obj,:clusterData)
            sp_obj.clustData.clustering = filter(:cell_id => ∈(cell_set), sp_obj.clustData.clustering)
            sp_obj.clustData.adj_mat = sp_obj.clustData.adj_mat[check_cell, check_cell]
        end
    end
    if isa(sp_obj, get_object_group("Imaging"))
        println("Updating spatial data...")
        sp_obj.spmetaData.cell = filter(:cell => ∈(cell_set), sp_obj.spmetaData.cell)
        sp_obj.spmetaData.molecule = filter(:cell => ∈(cell_set), sp_obj.spmetaData.molecule)
        sp_obj.spmetaData.molecule.gene = string.(sp_obj.spmetaData.molecule.gene)
        sp_obj.spmetaData.molecule = filter(:gene => ∈(gene_set), sp_obj.spmetaData.molecule)
        println("Updating polygons data...")
        prefix = Base.split(sp_obj.spmetaData.cell.cell[1],"_")
        if length(prefix) > 1
            poly_all_cell = prefix[1] .* "_" .* string.(sp_obj.spmetaData.polygon.polygon_number)
        else
            poly_all_cell = string.(sp_obj.spmetaData.polygon.polygon_number)
        end
        if isdefined(sp_obj.spmetaData, :polygon)
            if sp_obj.spmetaData.polygon !== nothing
                sp_obj.spmetaData.polygon = filter(:mapped_cell => ∈(cell_set), sp_obj.spmetaData.polygon)
            end
        end
        if length(prefix) > 1
            poly_cell = prefix[1] .* "_" .* string.(sp_obj.spmetaData.polygon.polygon_number)
        else
            poly_cell = string.(sp_obj.spmetaData.polygon.polygon_number)
        end
        poly_set = Set(poly_cell)
        check_poly_cell = [i in poly_set for i in poly_all_cell]
        if isdefined(sp_obj, :polyCount)
            sp_obj.polyCount = subset_count(sp_obj.polyCount; genes = genes, cells = poly_cell)
        end
        if isdefined(sp_obj, :polynormCount)
            sp_obj.polynormCount = subset_count(sp_obj.polynormCount; genes = genes, cells = poly_cell)
        end
        if isdefined(sp_obj, :polygonData)
            sp_obj.polygonData = sp_obj.polygonData[sp_obj.spmetaData.polygon.polygon_number]
            sp_obj.spmetaData.polygon.polygon_number = collect(1:length(sp_obj.spmetaData.polygon.polygon_number))
        end
        if isdefined(sp_obj, :imputeData)
            if isdefined(sp_obj.imputeData, :tgCount)
                if getfield(sp_obj.imputeData, :tgCount) !== nothing
                  sp_obj.imputeData.tgCount = subset_count(sp_obj.imputeData.tgCount; cells = cells)
                end
            end
            if isdefined(sp_obj.imputeData, :spageCount)
                if getfield(sp_obj.imputeData, :spageCount) !== nothing
                  sp_obj.imputeData.spageCount = subset_count(sp_obj.imputeData.spageCount; cells = cells)
                end
            end
            if isdefined(sp_obj.imputeData, :gimviCount)
                if getfield(sp_obj.imputeData, :gimviCount) !== nothing
                    sp_obj.imputeData.gimviCount = subset_count(sp_obj.imputeData.gimviCount; cells = cells)
                end
            end
        end
    end
    if isa(sp_obj, VisiumObject)
        println("Updating spatial data...")
        sp_obj.spmetaData.barcode = string.(sp_obj.spmetaData.barcode)
        sp_obj.spmetaData = filter(:barcode => ∈(cell_set), sp_obj.spmetaData)
    end
    return sp_obj
end

function subset_object(sp_obj::get_object_group("All"); cells = nothing, genes = nothing)
    sp_obj.rawCount = subset_count(sp_obj.rawCount; genes = genes, cells = cells)
    sp_obj = update_object(sp_obj)
    return sp_obj
end

function check_dim(sp_obj::get_object_group("All"), field_name::Union{Symbol, String})
    if isa(field_name, String)
       field_name = Symbol(field_name)
    end
    count_x = getfield(sp_obj, field_name)
    check_length = length(count_x.gene_name) !== length(sp_obj.rawCount.gene_name) || length(count_x.cell_name) !== length(sp_obj.rawCount.cell_name)
    return check_length
   end
   
function update_count(sp_obj::get_object_group("All"), ct_name::Union{Symbol, String})
    cell_id = colnames(sp_obj)
    gene_id = rownames(sp_obj)
    if isa(ct_name, String)
        ct_name = Symbol(ct_name)
    end
    if isdefined(sp_obj, ct_name)
        if check_dim(sp_obj, ct_name)
            ct_obj1 = getfield(sp_obj, ct_name)
            ct_obj2 = subset_count(ct_obj1; genes= gene_id, cells = cell_id)
            replacefield!(sp_obj, ct_name, ct_obj1, ct_obj2)
        end
    end
    return sp_obj
end

#= This function is deprecated because it is too slow
function check_vec(vec1, vec2)
    diff_elm = setdiff(vec2, vec1)
   if length(diff_elm) == 0
       gene_keep = repeat([true], length(vec2))
   else
       gene_keep = [x ∈ diff_elm ? false : true for x in vec2]
   end
end
=#

function check_vec(vec1, vec2)
    vec1_set = Set(vec1)
    return [x in vec1_set for x in vec2]
end

function find_unique_indices(input_vec)
    T = eltype(input_vec)
    counts = Dict{T, Int}()
    first_occurrence_indices = Int[]
    for (i, x) in enumerate(input_vec)
        if !haskey(counts, x)
            push!(first_occurrence_indices, i)
        end
        counts[x] = get(counts, x, 0) + 1
    end
    return first_occurrence_indices
end

function subset_matrix(count_mat, gene_name, cell_name, min_gene, min_cell)
    row_sum = sum(count_mat, dims=2)
    col_sum = sum(count_mat, dims=1)
    gene_kept = findall(x -> x > min_cell, row_sum[:])
    cell_kept = findall(x -> x > min_gene, col_sum[:])
    gene_name = gene_name[gene_kept]
    cell_name = cell_name[cell_kept]
    count_mat = count_mat[gene_kept, cell_kept]
    if is_duplicates(gene_name)
        gene_kept2 = find_unique_indices(gene_name)
        gene_name = gene_name[gene_kept2]
        count_mat = count_mat[gene_kept2, :]
    end
    return count_mat, gene_name, cell_name
end

function sparse_r_to_jl(counts::RObject{S4Sxp})
    i = rcopy(counts[:i]) .+ 1
    p = rcopy(counts[:p])
    x = rcopy(counts[:x])
    m = rcopy(counts[:Dim])[1]
    n = rcopy(counts[:Dim])[2]
    cols = vcat([fill(col, p[col + 1] - p[col]) for col in 1:n]...)
    julia_sparse_matrix = sparse(i, cols, x, m, n)
    return julia_sparse_matrix
end

function annotate_cells(sp::get_object_group("All"), cell_map::Dict; old_id_name::Union{String, Symbol}="cluster", new_id_name::Union{String, Symbol}="celltype")
    sp.metaData = map_values(sp.metaData, old_id_name , new_id_name, collect(keys(cell_map)),  collect(values(cell_map)))
    if isa(sp, Union{CartanaObject, MerfishObject, XeniumObject, seqFishObject, STARmapObject})
        sp.spmetaData.cell[!, old_id_name] = sp.metaData.cell[!, old_id_name]
        sp.spmetaData.cell = map_values(sp.spmetaData.cell,old_id_name , new_id_name, collect(keys(cell_map)),  collect(values(cell_map)))
        sp.spmetaData.molecule = map_values(sp.spmetaData.molecule, "cell", new_id_name, sp.spmetaData.cell[!, "cell"], sp.spmetaData.cell[new_id_name])
        if !isa(sp.spmetaData.polygon, Nothing)
            sp.spmetaData.polygon = map_values(sp.spmetaData.polygon, old_id_name, new_id_name, collect(keys(cell_map)),  collect(values(cell_map)))
        end
    end
    if isa(sp, Union{VisiumObject, SlideseqObject})
        sp.spmetaData[!, old_id_name] = sp.metaData[!, old_id_name]
        sp.spmetaData = map_values(sp.spmetaData, old_id_name, new_id_name, collect(keys(cell_map)),  collect(values(cell_map)))
    end
    return sp
end

function add_alpha_color_dict(color_dict::Dict{String, String}, alpha::Real)
    new_dict = Dict{String, Tuple{String, Real}}()
    for key in keys(color_dict)
        original_value = color_dict[key]
        new_value = (original_value, alpha)
        new_dict[key] = new_value
    end
    return new_dict
end

function add_alpha_color(color_vec::Dict{String, String}, alpha::Real)
    color_vec2 = map(x -> (x, alpha), color_vec)
    return color_vec2
end

## Below are the Harmony related functions. Please read the original paper for the algorithm details: 
## https://www.nature.com/articles/s41592-019-0619-0
## The Julia codes are based on a python implementation of Harmony (harmonypy):
## https://github.com/slowkow/harmonypy

function cluster_kmeans(data, K, random_state)
    Random.seed!(random_state)
    @info "Computing centroids with kmeans..."
    result = kmeans(data, K; maxiter=100, tol=1e-4)
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
#    if verbose && !is_converged
#        @info "Stopped before convergence"
#    end
end

function cluster!(harmony_obj::HarmonyObject)
    harmony_obj.dist_mat = 2 * (1 .- (harmony_obj.Y' * harmony_obj.Z_cos))
    j = 0
    for i in 1:harmony_obj.max_iter_kmeans
        harmony_obj.Y = harmony_obj.Z_cos * transpose(harmony_obj.R)
        column_norms = [norm(harmony_obj.Y[:, k], 2) for k in 1:size(harmony_obj.Y, 2)]
        harmony_obj.Y = harmony_obj.Y ./ column_norms'
        harmony_obj.dist_mat = 2 * (1 .- (harmony_obj.Y' * harmony_obj.Z_cos))
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