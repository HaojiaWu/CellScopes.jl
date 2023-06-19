# Tutorial for scATAC-seq analysis with CellScopes.jl

The following tutorial illustrates a standard analysis for scATAC-seq data. The current version only supports the data output produced by 10x cellranger-atac. As it continues to evolve, ```CellScopes.jl``` will support more single-cell chromatin modalities. The scATAC analysis below was inspired heavily by the [Signac package](https://stuartlab.org/signac/) from Stuart's lab.

## 1 Tutorial: Mouse kidney 12K cells
This tutorial uses a mouse kidney dataset for demo purposes. To test the scATAC-seq functionalies in ```CellScopes.jl``` using your own working environment, please download an example data from the [10x website](https://www.10xgenomics.com/resources/datasets?query=&page=1&configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000). 

### 1.1 Load data
We provided an easy-to-use function ```read_atac``` to directly read the cellrange-atac output into Julia and construct a ```scATACObject```. ```scATACObject``` is a novel data structure in ```CellScopes.jl``` designed for storing the original and processed data to facilitate the downstream analysis. The only parameter required to construct a ```scATACObject``` is the path to the cellranger-atac output. Similar to the single cell RNA-seq analysis, you can also set the min_peak and min_cell parameters to filter cell cells and peaks, respectively.

```julia
import CellScopes as cs
atac_path = "/mnt/sdd/multiomics/atac_raw/m12hr_run/outs/"
@time atac_obj = cs.read_atac(atac_path; min_peak=2500)
```
This should read the peak count, peak annotation, and fragment file into Julia and construct a complete ```scATACObject```. Note that it will take about a while to complete this step since some files are big in size. The messages below shows that a ```scATACObject``` has been successfully constracted.
```julia
This step reads all information directly from cellranger-atac output for downstream analysis. It may take 10 - 15 mins to complete as certain files (e.g. the fragment file) can be large in size.
1/3 Reading peak count data...
1/3 Peak count was loaded!
2/3 Reading the peak annotation file...
2/3 Peak annotation was loaded!
3/3 Reading the fragment file...
3/3 Fragments were loaded!
scATACObject was successfully constructed!
542.737890 seconds (1.27 G allocations: 78.840 GiB, 2.45% gc time, 10.05% compilation time: 19% of which was recompilation)
scATACObject in CellScopes.jl
Peaks x Cells = 128626 x 12620
Available data:
- Raw count
- Fragment data
- Peak data
All fields:
- rawCount
- normCount
- scaleCount
- metaData
- varPeak
- dimReduction
- clustData
- peakAnno
- fragmentData
- activityData
- undefinedData
```
### 1.2 Preprocess the scATACObject
This step performs data normalization on the raw peak count. Shadowed by the normalization methods developed by Signac, the raw read count was normalized using the term frequency inverse document frequency approach (TF-IDF). Here we computes the TF-IDF matrix as log(TF x IDF) following the **Method 1** in the RunTFIDF function of Signac.
```julia
atac_obj = cs.run_tf_idf(atac_obj)
```

### 1.3 Find top features
We use a similar approach implemented in the ```FindTopFeatures``` function in Signac to identify the top peaks. Users can set the min_cutoff parameters to choose the n% features to be selected for top features. For example, "q35" means 35% of the top peaks to be selected for dimensional reduction analysis. 

```julia
atac_obj = cs.find_top_features(atac_obj; min_cutoff="q35")
```

### 1.4 Dimensional reduction
Next, we perform linear demansional reduction using latent semantic indexing (LSI) approach in Signac (TF-IDF + SVD are knowns as LSI). We then embed the cells into low dimensional space using the UMAP methods as described in scRNA-seq and used the same graph-based approach to partition the cell distnace matrix into cell clusters.
```julia
atac_obj = cs.run_svd(atac_obj; method=:svd, pratio=0.99, maxoutdim=20)
atac_obj = cs.run_umap(atac_obj; dims_use=2:20, min_dist=0.2)
atac_obj = cs.run_clustering(atac_obj; res=0.0015)
```

## 2 Basic data visualization
### 2.1 Cell visualization
 We leverage existing visualization methods designed for scRNA-seq data to effectively visualize scATAC-seq data. For example, the ```dim_plot``` function can be employed to visualize cell clusters in a similar manner as it is used for single-cell RNA-seq data.
```julia
cs.dim_plot(atac_obj; dim_type = "umap", marker_size = 4)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/atac_clusters.png" width="600"> <br>

### 2.2 Peak visualization
Similar to gene expression visualization in scRNA-seq, TF-IDF normaized peak count can be visualized by the ```feature_plot``` function.
```julia
cs.feature_plot(atac_obj, ["chr9_2999952_3000894"]; 
    order=false, marker_size = 4, 
    count_type ="norm", color_keys=("black","yellow","red"))
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/location_plot.png" width="600"> <br>

## 3 Gene activity
We calculate the gene activity by summing the fragments intersecting the gene body and promoter region. To create a gene activity matrix, we designed a function called ```compute_gene_activity``` to extract gene coordinates, and count the number of fragments for each cell that map to gene body and promoter region of each gene. We then constructed a ```GeneActivityObject``` to store the activity data for plotting.
```julia
atac_obj = cs.compute_gene_activity(atac_obj)
```
After the ```GeneActivityObject``` is constructed, we can visualize the gene activity very easily using the ```gene_activity_plot``` function.
```julia
cs.gene_activity_plot(atac_obj, ["Nphs2", "Lrp2","Umod", "Slc12a3","Aqp2", "Pecam1"]; 
    color_keys = ["azure2","coral1","darkred"], order = false)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/gene_activity.png" width="800"> <br>

## 4 Coverage plot
We reimplemented the CoveragePlot function in Signac using pure Julia codes. In brief, fragments along a genomic region (e.g. gene body) was averaged within a group and normalized by the cell number and read counts in each group. For visiualize, counts were downsampled by setting the ```downsample_rate``` parameter and smoothed by a rolling function with a user defined window size (default is 100).
```julia
gtf_file = "/home/haojiawu/tenx_dir/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz"
atac_obj = cs.add_genecode(atac_obj, gtf_file)
cs.coverage_plot(atac_obj, "Nphs2"; downsample_rate=0.1)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/coverage.png" width="600"> <br>

