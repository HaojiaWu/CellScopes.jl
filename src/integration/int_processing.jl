function run_harmony(sc_obj::get_object_group("All"), batch::Union{String, Symbol}; kwargs...)        
    @info "The run_harmony function is a Julia implementation of the Harmony for data integration. Please read the original paper for the algorithm details: https://www.nature.com/articles/s41592-019-0619-0. The Julia codes are based on a python implementation of Harmony (harmonypy): https://github.com/slowkow/harmonypy"
    pca_mat = sc_obj.dimReduction.pca.cell_embedding
    metadata = sc_obj.metaData
    if isa(batch, String)
        batch = Symbol(batch)
    end
    ho = HarmonyObject(pca_mat, metadata, batch; kwargs...)
    harmony_matrix = Matrix{Float64}(ho.Z_corr')
    sc_obj.dimReduction.pca.cell_embedding = harmony_matrix
    return sc_obj
end