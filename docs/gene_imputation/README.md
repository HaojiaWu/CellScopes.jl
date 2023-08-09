## Gene imputation for imaging-based spatial technology with SpaGE, Tangram and gimVI
Imaging based spatial technologies (such as Xenium, MERFISH, Cartana, or seqFISH+) differ from spot-based spatial transcriptomics in that they are not constrained by spatial resolution. However, they are often limited in gene throughput, usually restricted to a selection of a few hundred pre-chosen genes. To leverage this limitation, tools like SpaGE, Tangram and gimVI were developed to align spatial profiling measurements with common sc/snRNA-seq profiles and extends the spatial gene measurement to whole transcriptomes. Here we provide serveral functions in CellScopes to conveniently accomplish the gene imputation with those tools.
```julia
import CellScopes as cs
using CSV, DataFrames, GZip, CSVFiles, SparseArrays
```
### 1. Prepare scRNA-seq data
First, we need a single cell RNA-seq data for imputation. Here we download a kidney snRNA-seq data published from our lab for demo purpose. 
```julia
;wget -c "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139107&format=file&file=GSE139107%5FMouseIRI%5Fcontrol%2Edge%2Etxt%2Egz" -O control_count.dge.txt.gz
;wget -c "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139107&format=file&file=GSE139107%5FMouseIRI%2Emetadata%2Etxt%2Egz" -O meta_data.csv.gz
```
Then we follow the CellScopes [scRNA analysis tutorial](https://github.com/HaojiaWu/CellScopes.jl/tree/main/docs/scRNA_tutorial) to create a scRNAseq Object in CellScopes.

```julia
meta = CSV.File("meta_data.txt.gz", delim='\t') |> DataFrame;
rename!(meta, ["cell","group","replicate","celltype"]);
##### Here we use the sham sample from Yuhei's PNAS paper
meta_sham3 = filter([:replicate, :group] => (x, y) -> x == "1_2" && y == "Control", meta);
counts = CSV.read("control_count.dge.txt.gz", DataFrame);
rename!(counts, [["gene"] ; names(counts)[1:end-1]]);
counts = counts[!, [["gene"] ; String.(meta_sham3.cell)]];
### create scRNA object
genes = counts.gene
cells = names(counts)[2:end]
counts = counts[:, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(counts));
raw_count = cs.RawCountObject(count_df, cells, genes)
sham_sc = cs.scRNAObject(raw_count; meta_data = meta_sham3)
### clustering
sham_sc = cs.normalize_object(sham_sc)
sham_sc = cs.scale_object(sham_sc)
sham_sc = cs.find_variable_genes(sham_sc; nFeatures=2000)
sham_sc = cs.run_pca(sham_sc;  method=:svd, pratio = 0.9, maxoutdim = 20)
sham_sc = cs.run_umap(sham_sc; min_dist=0.4, n_neighbors=20)
sham_sc = cs.run_clustering(sham_sc; res=0.005)
```
We can quickly explore the clustering result.
<br>
:point_right: Cell clustering and umap:
<br>
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/kidney_sc.png" width="600"> 
<br>
:point_right: Gene expression
<br>
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/kidney_gene.png" width="800"> <br>

### 2. Prepare Cartana data


