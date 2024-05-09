# Analysis and visualization of Visium HD data with CellScopes.jl

Visium HD is an enhanced version of Visium that identifies transcripts arrayed across millions of 2 x 2 µm barcoded squares with no gaps between them. Space Ranger outputs data in multiple bin sizes: 2 µm, 8 µm, and 16 µm. CellScopes allows reading and aligning data from all these bin sizes with histology images for data analysis and visualization. With Visium HD, we incorporated layers in CellScopes to manage different bin sizes independently, allowing users to analyze and visualize data at each bin size without interference from other bins.

### 1. Example dataset
This tutorial uses a mouse brain dataset provided by 10x Genomics to illustrate how to use CellScopes to read, analyze, and visualize Visium HD data. You can download the example dataset from 10x website:
https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he

### 2. Importing data
The ```read_visiumHD``` function in CellScopes is designed to import the output from Space Ranger, including spatial coordinates, histology images, and scale factors for various bin sizes. Data from each bin size is stored in separate layers, which allow users to select data at different resolutions to analyze.

```Julia
import CellScopes as cs
hd_dir = "/mnt/sdb/visiumHD/hd_output/outs/"
hd = cs.read_visiumHD(hd_dir)
```
This will create a VisiumHDObject that hold all the information for downstream analysis.

```
1. loading 2um binned data...
Formatting cell polygons...
Progress: 100%[==================================================] Time: 0:00:5039m
1. 2um binned data loaded!
2. loading 8um binned data...
Formatting cell polygons...
Progress: 100%[==================================================] Time: 0:00:08
2. 8um binned data loaded!
3. loading 16um binned data...
Formatting cell polygons...
3. 16um binned data loaded!
VisiumHDObject in CellScopes.jl
Genes x Cells = 18988 x 397258
Available data:
- layerData
- rawCount
- normCount
- scaleCount
- metaData
- spmetaData
- varGene
- dimReduction
- clustData
- imageData
- alterImgData
- polygonData
- defaultData
All fields:
- layerData
- rawCount
- normCount
- scaleCount
- metaData
- spmetaData
- varGene
- dimReduction
- clustData
- imageData
- alterImgData
- polygonData
- defaultData
```

### 3. Layer selection
After the Visium HD object is successfully constructed, you can select the data layer to analyze. By default, the data from 8 µm was selected. You can switch to other bin sizes if you want. For example, to analyze the 2 µm and 16 µm bin sizes, simply specify these options in the "layer_slot" parameter of the ```set_default_layer``` function.
```Julia
hd = cs.set_default_layer(hd; layer_slot = "2_um")
## or
hd = cs.set_default_layer(hd; layer_slot = "16_um")
```
When you switch layers, all data from the current layer will be automatically updated and saved. You can use the ```default_layer``` function to check which layer CellScopes is currently using.

```Julia
cs.default_layer(hd)
```

### 4. (Optional) Clustering and dimensional reduction
This step is optional, as Space Ranger already provides clustering and UMAP results for the 8 µm and 16 µm bin sizes. If you wish to recluster the cells using CellScopes' internal functions, please refer to the cell clustering workflow detailed in the Visium SD tutorial available at:
https://github.com/HaojiaWu/CellScopes.jl/tree/main/docs/visium_tutorial <br >
Note that reclustering cells from the 8 µm bin size can be time-consuming, as the number of cells often exceeds 300,000.

### 5. Data visualization
The clustering results can be directly visualized using the standard plotting functions provided by CellScopes.
#### 5.1 Visualize the clusters
```julia
hd = cs.set_default_layer(hd; layer_slot = "8_um")
cs.sp_dim_plot(hd, "cluster"; width=1300, height=1000, 
        do_legend=true, legend_size=40, legend_fontsize=40)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd1.png" width="1000"> <br>

You may notice that the tissue orientation in the graph looks different from what is displayed in the Loupe Browser. This difference can be attributed to the use of CairoMakie as the backend plotting tool for most visualizations in CellScopes, which handles the coordinates differently. To align the orientation with the original tissue orientation, you can use the ```convert_image_data``` function. Here’s how you can adjust the orientation:
```julia
@time hd = cs.convert_image_data(hd)
```
After converting the image data, all graphs will use the adjusted coordinates for visualization. You can verify the results using the same command line.
```julia
hd = cs.set_default_layer(hd; layer_slot = "8_um")
cs.sp_dim_plot(hd, "cluster"; width=1300, height=1000, 
        do_legend=true, legend_size=40, legend_fontsize=40)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd2.png" width="1000"> <br>

We can also highlight the dentate gyrus (DG) of the hippocampus on the same graph (cluster 13 in the 16 µm or cluster 11 in the 8 µm clustering analysis.).
```julia
hd = cs.set_default_layer(hd; layer_slot = "16_um")
hd = cs.convert_image_data(hd; layer_slot = "16_um")
cs.sp_dim_plot(hd, "cluster"; width=1000, height=1200, 
          anno_color = Dict("13" => "green2") ,do_legend=true, img_res = "high", stroke_width=0.2, 
            cell_highlight = "13",legend_size=30, adjust_contrast =1, adjust_brightness=0.0, alpha=0.4)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd3.png" width="1000"> <br>

#### 5.2 Select a region of interest for detailed visualization
To crop a specific region of interest in your spatial data, you can initially use the ```plot_fov``` function to locate the coordinates of the four corners.
```julia
cs.plot_fov(hd, 20, 20; marker_size = 0, width=2000, height=2400)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd4.png" width="1000"> <br>

 Here's how you can subset the field of view (FOV) using specific square numbers for the corners:
```julia
df = cs.subset_fov(hd, [193,197,333,337], 20,20);
xlim, ylim = cs.get_subset_cood(df)
#((9518.974171506527, 17960.0047340106), (14384.796062010153, 20199.649061752563))
```
This initial cropping may not precisely align with the desired region. Based on the coordinates obtained, you can adjust the limits slightly to better fit the area you intend to focus on:
```julia
xlim = (9418, 17960)
ylim=(15500, 19799)
```
This method allows you to fine-tune the coordinates manually for an optimal close-up view of your region of interest.
```julia
cs.sp_dim_plot(hd, "cluster"; width=1800, height=1000, y_lims=ylim , x_lims=xlim, 
          do_legend=true, img_res = "high", stroke_width=0.4, 
            anno_color = Dict("11" => "slateblue1", "2"=>"green2","13"=>"yellow1"),
            cell_highlight =["11", "2","13"],legend_size=40, 
        adjust_contrast =1, adjust_brightness=0.0, alpha=1)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd5.png" width="1000"> <br>

#### 5.3 Compare the spatial resolution of different bin sizes
You can easily compare the spatial distributions of cell types across different bin sizes using CellScopes. Here's how you can analyze the distributions in 8 µm and 16 µm resolutions:
```julia
xlim = (11818, 16100)
ylim=(15500, 17599)
hd = cs.set_default_layer(hd; layer_slot = "8_um")
cs.sp_dim_plot(hd, "cluster"; width=1500, height=900, y_lims=ylim , x_lims=xlim, 
          do_legend=true, img_res = "high", stroke_width=0.4, 
            anno_color = Dict("11" => "slateblue1", "2"=>"green2"),
            cell_highlight =["11", "2"],legend_size=40, 
        adjust_contrast =1, adjust_brightness=0.0, alpha=0.7)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd6.png" width="1000"> <br>

```julia
hd = cs.set_default_layer(hd; layer_slot = "16_um")
cs.sp_dim_plot(hd, "cluster"; width=1500, height=900, y_lims=ylim , x_lims=xlim, 
          do_legend=true, img_res = "high", stroke_width=0.4, 
            anno_color = Dict("13" => "slateblue1", "16"=>"green2"),
            cell_highlight =["13", "16"],legend_size=40, 
        adjust_contrast =1, adjust_brightness=0.0, alpha=0.7)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd7.png" width="1000"> <br>

Since no clustering analysis is available for the 2 µm data, it cannot be directly visualized in the same manner. However, you can still compare the 2 µm, 8 µm, and 16 µm resolutions in terms of spatial gene expression distribution by adjusting the layer selection and plotting parameters accordingly in CellScopes. Here is how to use the ```sp_feature_plot``` to visualize the gene expression at various resolutions.

```julia
xlim = (11818, 16100)
ylim=(15500, 17599)
hd = cs.set_default_layer(hd; layer_slot="2_um")
hd = cs.convert_image_data(hd; layer_slot = "2_um")
cs.sp_feature_plot(hd, ["Pcp4"]; color_keys=["gray94", "cyan", "blue", "darkblue"],  
    width=800, height=500, x_lims=xlim , y_lims=ylim, img_res = "high", alpha=1)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd8.png" width="1000"> <br>

```julia
hd = cs.set_default_layer(hd; layer_slot="8_um")
cs.sp_feature_plot(hd, ["Pcp4"]; color_keys=["gray94", "cyan", "blue", "darkblue"],  
    width=800, height=500, x_lims=xlim , y_lims=ylim, img_res = "high", alpha=1)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd9.png" width="1000"> <br>

```julia
hd = cs.set_default_layer(hd; layer_slot="16_um")
cs.sp_feature_plot(hd, ["Pcp4"]; color_keys=["gray94", "cyan", "blue", "darkblue"],  
    width=800, height=500, x_lims=xlim , y_lims=ylim, img_res = "high", alpha=1)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/hd10.png" width="1000"> <br>

### 6. Exporting data
```julia
cs.save(hd; filename="visiumHD.jld2")
```
