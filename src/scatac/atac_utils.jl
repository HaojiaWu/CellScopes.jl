function get_gene_location(atac_obj::scATACObject, gene::String)
    if isa(atac_obj.fragmentData.genecode, Nothing)
        error("Gene annotation file is missing. Please input the gft file with the add_genecode function!")
    end
    genecode = atac_obj.fragmentData.genecode
    chr, start, stop = GeneticsMakie.findgene(gene, genecode)
    gene_loc = [chr, start, stop]
    return gene_loc
end

function add_genecode(atac_obj::scATACObject, gtf_path)
    h = ["seqnames", "source", "feature", "start", "end", "score", "strand", "phase", "info"]
    gencode = CSV.read(gtf_path, DataFrame; delim = "\t", comment = "#", header = h)
    GeneticsMakie.parsegtf!(gencode)
    select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_id, :gene_name, :gene_type, :transcript_id)
    atac_obj.fragmentData.genecode = gencode
    return atac_obj
end