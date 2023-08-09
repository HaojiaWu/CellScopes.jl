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
The files for the Cartana data can be downloaded from the GEO (GSE227044). To match the scRNA data, we use a sham kidney from the project for gene imputation.
```julia
molecules =  DataFrame(CSV.File("/mnt/sdc/cartana_weekend/Sham_Jia/segmentation.csv"))
count_df =  DataFrame(CSV.File("/mnt/sdc/cartana_weekend/Sham_Jia/segmentation_counts.tsv"))
cells =  DataFrame(CSV.File("/mnt/sdc/cartana_weekend/Sham_Jia/segmentation_cell_stats.csv"))

gene_name = count_df.gene
cell_name = string.(names(count_df)[2:end])
count_df = count_df[!, 2:end]
count_df = convert(SparseMatrixCSC{Int64, Int64},Matrix(count_df))
raw_count = cs.RawCountObject(count_df, cell_name, gene_name)
molecules.cell = string.(molecules.cell)
cells.cell = string.(cells.cell)
sham_sp = cs.CartanaObject(molecules, cells, raw_count;
    prefix = "sham", min_gene = 0, min_cell = 3)
```
### 3. Run SpaGE gene imputation
Before running this step, the SpaGE package need to be downloaded from the original github respoitory:
```julia
sham_sp = cs.run_spaGE(sham_sp, sham_sc, "/mnt/sdd/NC_revision/compare_seurat/SpaGE/"; 
                gene_list=["Slc4a4", "Wdr17"])
```

```julia
cs.sp_feature_plot(sham_sp, ["Wdr17"]; 
    use_imputed=true, 
    color_keys = ["gray94", "lemonchiffon", "red"])
```


