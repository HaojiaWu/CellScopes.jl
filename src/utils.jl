function colSums(df::DataFrame, row_range, col_range)
    new_df = sum.(eachcol(df[row_range, col_range]))
end

function rowSums(df::DataFrame, row_range, col_range)
    new_df = sum.(eachrow(df[row_range, col_range]))
end

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
