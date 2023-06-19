function add_genecode(atac_obj, gtf_path)
    h = ["seqnames", "source", "feature", "start", "end", "score", "strand", "phase", "info"]
    gencode = CSV.read(gtf_path, DataFrame; delim = "\t", comment = "#", header = h)
    GeneticsMakie.parsegtf!(gencode)
    select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_id, :gene_name, :gene_type, :transcript_id)
    atac_obj.genecodeData = gencode
    return atac_obj
end