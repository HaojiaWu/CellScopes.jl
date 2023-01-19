# SpaData.jl
## A Julia package to process the imaging-based spatial transcriptomics data
Computational tools are lacking for processing the spatial transcriptomics data leading to high quality data intepretation. The ```SpaData.jl``` package is designed for processing, analyzing, and visualizing the FISH-based spatial data. The tutorial below used the spatial data from CARTANA for demo but basically it can process any FISH-based methods (such as MERFISH) after slight data formatting. If you find ```SpaData.jl``` useful to your research, please cite our paper or this github page. <br>

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/SpaData.png" width="600"> <br>

### 1. How to install SpaData.jl in Julia 1.7.3?
```julia
julia> using Pkg
julia> Pkg.add("https://github.com/HaojiaWu/SpaData.jl.git")
```

### 2. Example data for a test?
If you don't have sample data at hands, you can download our data from GEO (GSE XXXX) to test the functionality of SpaData.jl. If you use your own data, please change the column names for the transcript and cell coordinates to "x" and "y". Please also make sure your cell and molecule files contain a column to annotate the genes (gene ID) and cells (cell ID).

### 3. What is SpaObj object?
SpaObj is a new Julia object/data type we constructed to store the original and processed data including the 2D spatial coordinates, gene and cell annotations, raw/normalized count matrices and all results from analysis. We defined a broad range of computational methods to process this SpaObj object, resulting in a wide range of high quality plots for data visualization.

### 4. How to construct a SpaObj object?
**4a. Read the spatial data into Julia.** Three files as DataFrame format are required. Molecule file is the file that contains the spatial coordinates for the detected transcripts. Cell file is the spatial coordinates after cell segmentation. Count file is the gene-by-cell count matrix from cell segmentation. <br>
```julia
using Pkg, DataFrames, CSV, StatsBase, VegaLite
import SpaData as spd
molecules =  DataFrame(CSV.File("/mnt/sdc/cartana_weekend/Day2_Jia/segmentation.csv"));
count_df =  DataFrame(CSV.File("/mnt/sdc/cartana_weekend/Day2_Jia/segmentation_counts.tsv"));
cells =  DataFrame(CSV.File("/mnt/sdc/cartana_weekend/Day2_Jia/segmentation_cell_stats.csv"));
```
**4b. Construct a SpaObj object** using the count data and the cell and molecules coordinates data. You can normalize the data in this step or do it later. If you want to normalize the data in this step, set ```normalize=true```. <br>
```julia
kidney = spd.SpaObj(count_molecules,count_cells,count_df; cell_prefix="Day2", normalize=true);
```
**4c. (Optional). Normalize the raw count data** if ```normalize=false``` in 3b. <br>
```julia
kidney = spd.normalizeData(kidney);
```

### 5. SpaObj file I/O
We used the [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) package to save and load the ```SpaObj``` object. Please read the original tutorial to learn the methods from JLD2 for file input/output. Here are some example lines:
```julia
using JLD2
### save SpaObj object to disk.
save("kidney_SpaObj.jld2","kidney", kidney) 
### read SpaObj object from disk.
kidney = load("kidney_SpaObj.jld2")
kidney = kidney["kidney"]
```

### 6. Data processing and analysis
#### 6a. Cell clustering
We wrapped the codes from the popular sc tools to cluster the cells based on the transcript count matrix from cell segmentation. Currently, we provided [SpaGCN](https://github.com/jianhuupenn/SpaGCN), [Seurat](https://github.com/satijalab/seurat), and [Scanpy](https://github.com/scverse/scanpy) for cell clustering. You can skip this step if your data have completed cell clustering by other tools (e.g. Baysor). <br/>

i. **SpaGCN**: Please refer to the original tutorial to select the paramenters. Below are some example codes (Clustering results will be stored in the cell metadata of the SpaObj):
```julia
CSV.write("count_data.csv",spd.count)
kidney = spd.run_SpaGCN(kidney, "count_data.csv", "/home/users/haojiawu/anaconda3/bin/python")
```
ii. Here is how to run **Scanpy** in SpaData:
```julia
adata = spd.run_scanpy(kidney.count)
```
iii. To run **Seurat**:
```julia
seu_obj = spd.run_seurat(kidney.count)
```

#### 6b. Cell polygons and mapping
We used Baysor to create a polygon object and store it in SpaData for cell drawing. Baysor is required to install before running this step. Please refer to the original repo here: https://github.com/kharchenkolab/Baysor
```julia
import Baysor as B
scale=25;
min_pixels_per_cell = 15;
grid_step = scale / min_pixels_per_cell;
bandwidth= scale / 10;
polygons = B.boundary_polygons(kidney.molecules, kidney.molecules.cell, grid_step=grid_step, bandwidth=bandwidth);
kidney.polygons = polygons
```
Then the gene expression can be mapped to the cell polygons by approximation based on the Elucidean distance between the original cell coordinate and the center of the ploygons. This can be done with the ```polygons_cell_mapping``` function.
```julia
kidney = spd.polygons_cell_mapping(kidney)
```
Based on the cell mapping results, we can obtain a polygons count matrix. this count matrix can be used for ploting gene expression directly on cell polygons (see the Visualization section).
```julia
kidney = spd.generate_polygon_counts(kidney)
```

#### 6c. Calculate cell-cell distance
To reveal the phycical cell-cell contact, we provided two functions to calculate the distance of any given cell populations depending on their cell distribution patterns.

**i.** When the cells are confined to some specific regions (such as the glomerular cell types), we take a cell-centric approach to calculate the cell-cell distance. We take each cell from the cell population of interest, and compute the distance between this cell and the cells from other cell types. This process will be iteratively repeated until all cells from the cell type of interest are done. For example, if we want to measure which EC subtype is close to the podocyte within a given search area (radius), we can run the step below:
```julia
cell_dist = spd.compare_cell_distances(kidney, :celltype, "Podo", "gEC", "vEC", 50)
```
**ii.** When the cell distribution is very diffusive (such as the immune cells or fibroblasts), we use a cell-enrichment approach to compute the patial proximity of pairs of cell types as reported by [Lu et al](https://www.nature.com/articles/s41421-021-00266-1). We calculate the probability of cell type pairs in a neighborhood within a given searching area (radius =50). We then compute the enrichment of cell type pairs in spatial proximity after normalized to the control probability based on random pairing. 
```julia
cell_dist = spd.run_cell_pairing(kidney.cells, :cell2, :celltype, 50)
```
#### 6d. Convert the xy coordinates to kidney coordinates
In ```SpaData.jl```, we created a new coordinate system, namely **kidney coordinate system**, to precisely depict the position of every single cell in the kidney. In this system, the position of a cell is defined by the kidney depth, and the kidney angle. To transform the xy coordinate system to kidney coordinate system, we first define the origin of the coordinate by finding the center point in the papilla. For each cell, we compute the kidney depth by calculating the distance of the cell to the kidney boundary, and divided by the distance of the kidney boundary to the origin of the coordinate. We can define the kidney angle of the cells by measuring the angle of the slope and the new x coordinate (in tangent value). The schematic below explains the coordinate transformation.

<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/kidney_coordinate.png" width="300"> <br>

This kidney coordinate system can help define the kidney compartment where the cell resides, how the cell type and transcript distribution changes from outer cortex to papilla, and how the gene expression changes in different condiitons. Here are some steps to complete the transformation.
```julia
spd.plot_fov(kidney, 40,40; group_label="celltype", cell_highlight="CD-PC", shield=true) ### use grid to find the papilla area
cell_sub=spd.subset_fov(kidney, [661,671,941,951], 40,40); ### select the papilla area
cell_sub=filter(:celltype => x -> x =="CD-PC", cell_sub); ### select the PC cells in the papilla
center=[mean(cell_sub.x),mean(cell_sub.y)]; ### find the center point as origin
kidney = spd.compute_kidney_coordinates(kidney, center);
```
#### 6e. scRNA-seq integration
```SpaData.jl``` uses [Tangram](https://github.com/broadinstitute/Tangram) to integrate the scRNA-seq data. Please refer to the tangram [tutorial](https://tangram-sc.readthedocs.io/en/latest/) for more detail. We wrapped the python codes into a julia function as below. ```SpaData.jl``` also provides function to plot the results (See the Visualization session).
```julia
using SpaData, CSV, DataFrames, JLD2, PyCall,Pkg
ENV["PYTHON"]="/home/users/haojiawu/anaconda3/envs/tangram-env/bin/python"
Pkg.build("PyCall")
kidney=spd.run_tangram(kidney, "/home/users/haojiawu/anaconda3/envs/tangram-env/bin/python",
    "IRI_sham/sc_matrix_w6.mtx",
    "IRI_sham/sc_genes_w6.tsv",
    "IRI_sham/sc_barcodes_w6.tsv",
    "IRI_sham/sc_metadat_w6.csv",
    "markers_mapping.csv",
    "name"
)
```

### 7. Visualization
We provided a number of functions to visualize the results from the above analysis. 
#### 7a. Plot gene expression on segmented cells. 
Gene expression can be easily visualized by the ```featurePlot``` function. Here is an example for ploting the whole kidney.
```julia
alpha_trans=1
anno2 = Dict("Podo" => ("magenta1",alpha_trans), "HealthyPT"=>("green3",alpha_trans), "InjPT"=>("#f92874",alpha_trans),"TAL"=>("lightslateblue",alpha_trans),"DCT"=>("blue",alpha_trans),"CD-PC"=>("turquoise1",alpha_trans),"CD-IC"=>("#924cfa",alpha_trans),"vEC"=>("firebrick",alpha_trans),"gEC"=>("dodgerblue",alpha_trans),"Fib"=>("#edff4d",alpha_trans),"JGA"=>("sienna2",alpha_trans),"Immune"=>("darkgreen",alpha_trans),"Uro"=>("black",alpha_trans));
p3=spd.plot_cell_polygons(kidney, "celltype"; 
    anno_color=anno2,x_lims=(0,35000), 
    y_lims=(0,40000),canvas_size=(5000,6000),
    stroke_color="gray80")
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/whole.jpg" width="300"> <br>

Usually it's hard to see the delicate tructure when ploting gene on the whole kidney. Therefore, we provided three ways to plot gene expression in a selected field of view. <br/>

```julia
##i.plot gene on cells as data points
p1=spd.featurePlot(kidney, ["Podxl"]; 
    layer="molecules",molecule_colors=["red"], x_lims=(21600,24400), y_lims=(4200,6200),order=true,
    pt_size=20,fig_width=600, fig_height=500)
##ii. plot gene on cells as segmented polygons
cmap=ColorSchemes.ColorScheme([colorant"gray98",colorant"red", colorant"red4"])
p2=spd.plot_gene_polygons(kidney, "Podxl",cmap; 
    x_lims=(21600,24400), y_lims=(4200,6200),
    canvas_size=(600,500))
##iii. plot transcripts on cells as segmented polygons
alpha_trans=0.5
anno2 = Dict("Podo" => ("fuchsia",alpha_trans), "HealthyPT"=>("green",alpha_trans),"InjPT"=>("lime",alpha_trans),"TAL"=>("cyan4",alpha_trans),"DCT"=>("yellow",alpha_trans),"CD-PC"=>("gray95",alpha_trans),
            "CD-IC"=>("gray95",alpha_trans),"aEC"=>("red",alpha_trans),"gEC"=>("blue",alpha_trans),"Fib"=>("gray95",alpha_trans),"MC"=>("coral",alpha_trans),"Immune"=>("gray95",alpha_trans),"Uro"=>("gray95",alpha_trans))
@time spd.plot_transcript_polygons(kidney; 
    genes=["Podxl","Ehd3","Ren1"], colors=["fuchsia","blue","coral"],canvas_size=(600,600),
    markersize=2,annotation=:celltype, ann_colors=anno2, is_noise=:is_noise,noise_kwargs=(markersize=0, color="transparent"),
    show_legend=false,bg_color="transparent",x_lims=(17600,20000), y_lims=(7700,10000),segline_size=1, transparency=0.3
)
```
<p float="left">
  <img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/datapoints.png" width=32% height=250>
  <img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/polygons.png" width=32% height=250> 
  <img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/glom.png" width=32% height=250>
</p>

#### 7b. Plot cell annotation
Here are some examples to show how to visualize the kidney structure in different region.
```julia
### tubule
alpha_trans=1
anno2 = Dict("Podo" => ("magenta1",alpha_trans), "HealthyPT"=>("green3",alpha_trans),"InjPT"=>("#f92874",alpha_trans),"TAL"=>("lightslateblue",alpha_trans),"DCT"=>("blue",alpha_trans),"CD-PC"=>("turquoise1",alpha_trans),"CD-IC"=>("#924cfa",alpha_trans),"vEC"=>("firebrick",alpha_trans),"gEC"=>("dodgerblue",alpha_trans),"Fib"=>("#edff4d",0.5),"JGA"=>("sienna2",alpha_trans),"Immune"=>("darkgreen",alpha_trans),"Uro"=>("black",alpha_trans));
spd.plot_cell_polygons(kidney, "celltype"; 
    anno_color=anno2, x_lims=(18300,19800), y_lims=(10700,14000), 
    canvas_size=(300,600),stroke_color="gray80")

### renal artery
alpha_trans=1
anno2 = Dict("Podo" => ("white",alpha_trans), "HealthyPT"=>("white",alpha_trans),"InjPT"=>("white",alpha_trans),"TAL"=>("white",alpha_trans),"DCT"=>("white",alpha_trans),"CD-PC"=>("white",alpha_trans),"CD-IC"=>("white",alpha_trans),"vEC"=>("firebrick",alpha_trans),"gEC"=>("white",alpha_trans),"Fib"=>("#edff4d",alpha_trans),"JGA"=>("white",alpha_trans),"Immune"=>("white",alpha_trans),"Uro"=>("white",alpha_trans));
spd.plot_cell_polygons(kidney, "celltype"; 
    anno_color=anno2, x_lims=(21900,24500), y_lims=(11500,19500),
    canvas_size=(300,900),stroke_color="gray80")

### renal cortex
alpha_trans=1
anno2 = Dict("Podo" => ("magenta1",alpha_trans), "HealthyPT"=>("green3",alpha_trans),"InjPT"=>("#f92874",alpha_trans),"TAL"=>("lightslateblue",alpha_trans),"DCT"=>("blue",alpha_trans),"CD-PC"=>("turquoise1",alpha_trans),"CD-IC"=>("#924cfa",alpha_trans),"vEC"=>("firebrick",alpha_trans),"gEC"=>("dodgerblue",alpha_trans),"Fib"=>("#edff4d",0.5),"JGA"=>("sienna2",alpha_trans),"Immune"=>("darkgreen",alpha_trans),"Uro"=>("black",alpha_trans));
spd.plot_cell_polygons(kidney, "celltype"; 
    anno_color=anno2, x_lims=(23000,24800), y_lims=(7400,9000), 
    canvas_size=(450,400),stroke_color="gray80")
```
<p float="left">
  <img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/tubule.png" height=400>
  <img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/artery.png" height=400> 
  <img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/cortex.png" height=400> 
</p>

#### 7c. Plot gene expression across clusters.
We provided a gene rank function to identify the marker genes for each cluster. If cell types were annotated, the ```dotPlot``` function can visualize the gene expression across cell types.

**i.** plot gene rank
```julia
spd.plot_marker_rank(kidney, "celltype","vEC"; num_gene=20)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/gene_rank.png" height="200"> <br>

**ii.** Plot gene expression in dotplot format
```julia
cell_order=["Podo", "HealthyPT", "InjPT","TAL","DCT","CD-PC","CD-IC","Uro","gEC","aEC","MC","Fib","Immune"];
genes=["Podxl","Slc26a4","Havcr1","Slc5a2","Krt19","Aqp2","Slc12a3","Eln","Ehd3","Acta2","Col1a1"]
spd.dotPlot(kidney, genes, :celltype; cell_order=cell_order, expr_cutoff=0.1,fontsize=16,fig_height=500, fig_width=300)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/dotplot2.png" width="400"> <br>

#### 7d. Plot cell fraction.
We also provided a way to plot cell fraction across multiple conditions.
```julia
using VegaLite
sham_cell=spd.make_cell_proportion_df(sham.cells; nfeatures=0)
sham_cell.time.="Sham";
hour4_cell=spd.make_cell_proportion_df(hour4.cells; nfeatures=0)
hour4_cell.time.="Hour4";
hour12_cell=spd.make_cell_proportion_df(hour12.cells; nfeatures=0)
hour12_cell.time.="Hour12";
day2_cell=spd.make_cell_proportion_df(day2.cells; nfeatures=0)
day2_cell.time.="Day2";
week6_cell=spd.make_cell_proportion_df(week6.cells; nfeatures=0)
week6_cell.time.="Week6";
all_time=[sham_cell; hour4_cell; hour12_cell; day2_cell; week6_cell];
p=all_time |> @vlplot()+@vlplot(mark={:area, opacity=0.6}, x={"index", axis={grid=false} }, y={:fraction, stack=:zero, axis={grid=false}}, color={"celltype2:n",  scale={
            domain=cell_order2,
            range=cell_color
        }},width=200,height=200) +
@vlplot(mark={:bar, width=1, opacity=1}, x={"index", title="Time point"}, y={:fraction, stack=:zero, title="Cell proportion"}, color={"celltype2:n",  scale={
            domain=cell_order2,
            range=cell_color
        }},
    width=200,height=200)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/cellfrac.png" width="400"> <br>

#### 7e. Select and plot field of view (fov).
```SpaData.jl``` allows you to select the field of view for further analysis. First, we provided a function to draw grid on the spatial graph. Then the fov of interest can be selected using the ```subset_fov``` function.
```julia
spd.plot_fov(kidney, 10,10; group_label="celltype", cell_highlight="CD-PC", shield=true)
cell_sub=spd.subset_fov(kidney, [47,48,57,58], 10,10);
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/grid.jpg" width="400"> <br>

#### 7f. Plot cell type/transcript distribution from cortex to papilla.
After we transform the xy coordinates to the kidney coordinates (See the **Data processing and analysis** section), we can plot the cell type and transcript distribution from outer cortex to papilla.
```julia
markers=["Slc5a2","Slc12a3","Podxl","Slc26a4","Ehd3","Havcr1","Cd74","Ren1","Col1a1",
"Eln","Umod","Krt19","Aqp2"]
markers=reverse(markers);
celltypes=["HealthyPT","DCT","Podo","CD-IC","gEC","InjPT","Immune","gEC","Fib",
"vEC","TAL","Uro","CD-PC"]
celltypes=reverse(celltypes);
spd.plot_depth(kidney, celltypes=celltypes, markers=markers)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/kidney_depth.png" height="300"> <br>
We can make this into animation too.
```julia
spd.plot_depth_animation(kidney, celltypes=celltypes, markers=markers)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/cartana_tutorial/img/animations.gif" height="300"> <br>

