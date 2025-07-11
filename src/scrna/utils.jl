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
    cl_data = deepcopy(sc_obj.clustData.clustering)
    if isa(anno, String)
        anno = Symbol(anno)
    end
    cl_data=filter(anno => ==(cl), cl_data)
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
    sp_obj.metaData=filter(:Cell_id => ∈(cell_set), sp_obj.metaData)
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
            sp_obj.clustData.clustering=filter(:cell_id => ∈(cell_set), sp_obj.clustData.clustering)
            sp_obj.clustData.adj_mat = sp_obj.clustData.adj_mat[check_cell, check_cell]
        end
    end
    if isa(sp_obj, get_object_group("Imaging"))
        println("Updating spatial data...")
        sp_obj.spmetaData.cell=filter(:cell => ∈(cell_set), sp_obj.spmetaData.cell)
        sp_obj.spmetaData.molecule=filter(:cell => ∈(cell_set), sp_obj.spmetaData.molecule)
        sp_obj.spmetaData.molecule.gene = string.(sp_obj.spmetaData.molecule.gene)
        sp_obj.spmetaData.molecule=filter(:gene => ∈(gene_set), sp_obj.spmetaData.molecule)
        println("Updating polygons data...")
        prefix = Base.split(sp_obj.spmetaData.cell.cell[1],"_")
        if length(prefix) > 1
            poly_all_cell = prefix[1] .* "_" .* string.(sp_obj.spmetaData.polygon.polygon_number)
        else
            poly_all_cell = string.(sp_obj.spmetaData.polygon.polygon_number)
        end
        if isdefined(sp_obj.spmetaData, :polygon)
            if sp_obj.spmetaData.polygon !== nothing
                sp_obj.spmetaData.polygon=filter(:mapped_cell => ∈(cell_set), sp_obj.spmetaData.polygon)
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
        sp_obj.spmetaData=filter(:barcode => ∈(cell_set), sp_obj.spmetaData)
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

function match_patterns(item, patterns)
    for pattern in patterns
        if occursin(pattern, item)
            return true
        end
    end
    return false
end

function is_match(all_genes::Vector{String}, patterns::Vector{Regex})
    for gene in all_genes
        for pattern in patterns
            if occursin(pattern, gene)
                return true
            end
        end
    end
    return false
end

function fraction_expr_per_cell(expr::SparseMatrixCSC{<:Real}, gene_names::Vector{String})
    col_sums = vec(sum(expr, dims=1))
    col_sums[col_sums .== 0] .= 1e-12

    result = copy(expr)
    for col in 1:expr.n
        range = expr.colptr[col]:(expr.colptr[col+1]-1)
        result.nzval[range] ./= col_sums[col]
    end
    row_means = vec(sum(result, dims=2)) ./ size(result, 2)
    ordered_indices = sortperm(row_means; rev=true)
    ordered_genes = DataFrame(gene_name = gene_names[ordered_indices],
                              original_index = ordered_indices)
    return result, ordered_genes
end

function top_gene_fraction_df(fraction_matrix::SparseMatrixCSC{<:Real}, 
                              ordered_genes::DataFrame; 
                              top_n::Int = 18)
    top_indices = ordered_genes.original_index[1:top_n]
    top_gene_names = ordered_genes.gene_name[1:top_n]
    sub_mat = fraction_matrix[top_indices, :]
    dense_mat = Matrix(sub_mat)'
    df = DataFrame(dense_mat, Symbol.(top_gene_names))
    return df
end