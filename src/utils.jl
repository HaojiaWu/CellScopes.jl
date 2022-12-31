colSum(mtx::AbstractMatrix{<:Real}) = sum(mtx, dims=1)
rowSum(mtx::AbstractMatrix{<:Real}) = sum(mtx, dims=2)
rownames(sc_obj::scRNAObject) = sc_obj.rawCount.gene_name
colnames(sc_obj::scRNAObject) = sc_obj.rawCount.cell_name
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

function SubsetCount(ct_obj::T; 
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
    check_gene = [i in(genes) for i in all_genes]
    check_cell = [i in(cells) for i in all_cells]
    gene_name = all_genes[check_gene]
    cell_name = all_cells[check_cell]
    new_count = ct_obj.count_mtx[check_gene, check_cell]
    if isa(ct_obj, RawCountObject)
        new_obj = RawCountObject(new_count, cell_name, gene_name)
    elseif isa(ct_obj, NormCountObject)
        new_obj = NormCountObject(new_count, cell_name, gene_name,  ct_obj.scale_factor, ct_obj.norm_method, ct_obj.pseudocount)
    else
        new_obj = ScaleCountObject(new_count, cell_name, gene_name,  ct_obj.do_scale, ct_obj.do_center, ct_obj.scale_max)
    end
    return new_obj
end

function ExtractClusterCount(sc_obj::scRNAObject, cl; count_type = "norm")
    df = sc_obj.clustData.clustering
    cl_data = filter(:cluster => x -> x == cl, df)
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
    cl_ct = SubsetCount(ct_obj; cells = cl_cell)
end


