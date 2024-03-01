# Use a high-resolution histology/DAPI image for Xesium data visualization
```CellScopes.jl``` provides functions that allow for the integration of an external histology or immunofluorescence image for visualizing Xenium data. Given that Xenium data provide three layers of spatial information (transcripts, genes, and cells), this tutorial will walk you through how to visualize this information on top of a high-resolution histology or DAPI image.

### 1. Input image: a registered H&E or DAPI image
This tutorial requires a pre-registered H&E/DAPI image that has been accurately scaled and aligned with the Xenium transcript coordinates. The process of image registration can be quite challenging. For H&E or DAPI images derived from the identical section that has undergone Xenium processing, linear transformation (affine) is generally sufficient. However, when attempting to register images from consecutive tissue sections, more advanced algorithms, including those based on deep learning, are often required to achieve accurate alignment. This tutorial will not cover image registration. Users are encouraged to explore the existing methods to achieve a well-registered image before beginning this tutorial. For demo purpose, we will use the pre-registered H&E and DAPI image provided by the Xenium dataset from 10x website (https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast).

### 2. Xenium dataset
For this tutorial, we will analyze a breast cancer Xenium dataset from 10x Genomics. The dataset can be downloaded from here: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast. 
```julia
import CellScopes as cs
breast_cancer = read_xenium("/mnt/sdb/breast_cancer_xenium/")
```
This dataset come with a H&E image taken on that same Xenium slide. Xenium analyzer also output a DAPI tiff image that has been aligned to the xenium transcript coordinates. We will use these two images to showcase the ploting functions provided by CellScopes. 

```julia
he_file = "he.tif"
dapi_file = "if.tif"
```

### 3. Add H&E image to the xenium object.
We provide a function that can add any pre-aligned image taken by confocal microscope to the existing xenium object. Here is an example to add the H&E image.
```julia
breast_cancer = cs.add_xenium_img(breast_cancer;img_path=he_file)
```
#### 3a. Visualize cell type annotation on the H&E image.
Once the H&E image is added to the Xenium object, we can overlay the cell type annotations to the histological image using the ```sp_dim_plot``` function in ```CellScopes```. By default, the ```custom_img``` parameter is set to ```false```. We need to change the value to ```true``` in order to display the H&E image. Another important parameter is ```adjust_coord_to_img``` which adjusts the coordinates based on the image dimension. In most case, setting it to ```auto``` will be fine. The ```adjust_contrast``` and ```adjust_brightness``` parameters allow for the adjustment of the contrast and brightness of the H&E image. ```adjust_contrast=1.0``` and ```adjust_brightness=0.0``` indicates that there will be no adjustments made..
```julia
cs.sp_dim_plot(breast_cancer, "cluster"; 
    do_label = false,marker_size = 5, 
    width = 2500, height=2000, do_legend=false, custom_img=true, adjust_coord_to_img="auto", 
    adjust_contrast = 1.0, adjust_brightness = 0.0
    )
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/merge.jpg" width="800"> <br>

#### 3b. Visualize cell type annotation on a cropped field of view.
To focus on a specific tissue structure and examine the distribution of cell types within it, simply specify the x and y coordinate limits of the desired region to be cropped in the ```sp_dim_plot``` function. 
```julia
x1 = 5000
x2 = 5700
y1 = 800
y2 = 1600
cs.sp_dim_plot(breast_cancer, "cluster"; 
    do_label = false, marker_size = 8, x_lims = (x1, x2),
    y_lims = (y1, y2), width=600, height=650, do_legend=false, alpha=0.5,
    custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 1.0, adjust_brightness = 0.0
    )
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/cell.png" width="600"> <br>
<br>
The image above has no H&E image. To add the H&E image to the background, we can do:
```julia
### Cell annotation on H&E
cs.sp_dim_plot(breast_cancer, "cluster"; 
    do_label = false, marker_size = 8, x_lims = (x1, x2),
    y_lims = (y1, y2), width=600, height=650, do_legend=false, alpha=0.5,
    custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 1.0, adjust_brightness = 0.0
    )
### H&E only
cs.sp_dim_plot(breast_cancer, "cluster"; 
    do_label = false, marker_size = 0, x_lims = (x1, x2),
    y_lims = (y1, y2), width=600, height=650, do_legend=false, alpha=0.5,
    custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 1.0, adjust_brightness = 0.0
    )
```
<p float="left">
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/cell_he.png" width=400>
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/he.png" width=400>
</p>

#### 3c. How to find the limits of x, y coordinates?
In most cases, determining the precise x and y coordinate limits for the field of view (FOV) to be cropped can be quite challenging. ```CellScopes``` provides functions designed to simplify this process. The approach involves creating a numbered grid system over the graph, with each square of the grid assigned a unique number. The following image illustrates how this grid system works.
```julia
cs.plot_fov(breast_cancer, 10, 10)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/grid.png" width="600"> <br>
Then just select the numbers, or the numbers at the 4 corners (in numerical order) of a square or rectangle to crop your fov of interest, and provide these numbers to the ```subset_fov``` function. Here is an example how we get the x_lims and y_lims to visualize the structure provided in this tutorial.
```julia
fov1 = cs.subset_fov(breast_cancer, [62, 63, 72, 73],10,10)
x1 = minimum(fov1.x)
x2 = maximum(fov1.x)
y1 = minimum(fov1.y)
y2 = maximum(fov1.y)
```
### 4. Add the DAPI image to the xenium object.
Similarly, we can add the DAPI image to the xenium object and project the cells or gene expression on the DAPI image.
```julia
breast_cancer = cs.add_xenium_img(breast_cancer;img_path=dapi_file)
```
The subsequent sections of this tutorial will use this DAPI image to demonstrate how to overlay various layers of spatial information, including cell annotations, gene expression, and transcripts, onto the image. Note that the codes provided below are applicable to any type of pre-aligned images that have been added to the Xenium object.

#### 4a. Visualize cell type annotation on the DAPI image.
Below are the example codes to show cell annotations within the same cropped area  highlighted prevsiously.
```julia
### Cell annotations on DAPI
cs.sp_dim_plot(breast_cancer, "cluster"; 
    do_label = false, marker_size = 8, x_lims = (x1, x2),
    y_lims = (y1, y2), width=600, height=650, do_legend=false, alpha=0.5,
    custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 4.0, adjust_brightness = 0.3
    )
### DAPI only
cs.sp_dim_plot(breast_cancer, "cluster"; 
    do_label = false, marker_size = 8, x_lims = (x1, x2),
    y_lims = (y1, y2), width=600, height=650, do_legend=false, alpha=0.5,
    custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 4.0, adjust_brightness = 0.3
    )
```
<p float="left">
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/cell_dapi.png" width=400>
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/dapi.png" width=400>
</p>

#### 4b. Visualize gene expression on the DAPI image.
With the function ```sp_feature_plot```, we can overlay gene expression on the DAPI image. Here is the way to do this.
```julia
### DAPI off
cs.sp_feature_plot(breast_cancer, "TACSTD2"; x_lims = (x1, x2),y_lims = (y1, y2),
    color_keys=["gray94", "thistle1", "magenta"], 
    height=650, width=600, marker_size = 8)
### DAPI on
cs.sp_feature_plot(breast_cancer, "TACSTD2"; x_lims = (x1, x2),y_lims = (y1, y2),
    color_keys=["gray94", "thistle1", "magenta"], 
    height=650, width=600, marker_size = 8, clip=4,
    custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 4.0, adjust_brightness = 0.3)
```
<p float="left">
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/gene.png" width=400>
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/gene_dapi.png" width=400>
</p>

Visualization of multiple genes simultaneously is also possible. Simply provide a list of genes, rather than a single gene, to the ```sp_feature_plot``` function.
```julia
### DAPI off
cs.sp_feature_plot(breast_cancer, ["ACTA2","TACSTD2"]; x_lims = (x1, x2),y_lims = (y1, y2),
    color_keys=["gray94", "thistle1", "magenta"], 
    height=650, width=600, marker_size = 8)
### DAPI on
cs.sp_feature_plot(breast_cancer, ["ACTA2","TACSTD2"]; x_lims = (x1, x2),y_lims = (y1, y2),
    color_keys=["gray94", "thistle1", "magenta"], 
    height=650, width=600, marker_size = 8, clip=4,
    custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 4.0, adjust_brightness = 0.3)
```
<p float="top">
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/two%20genes.png" height=400>
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/twogenes_dapi.png" height=400>
</p>

#### 4c. Visualize transcripts on the DAPI image.
Finally, we can visualize the transcripts directly on the DAPI image. Given that transcripts are individual spots, we have the option to overlay multiple transcripts on a single graph or to distribute the transcripts across multiple graphs. To combine multiple transcripts within a single graph, please follow the codes below.
```julia
### DAPI off
cs.sp_feature_plot(breast_cancer, ["ACTA2","TACSTD2"]; layer="transcripts",
    x_lims = (x1, x2),y_lims = (y1, y2), alpha=0.5,
    gene_colors=["cyan", "magenta"], 
    height=650, width=600, marker_size = 2)
### DAPI on
cs.sp_feature_plot(breast_cancer, ["ACTA2","TACSTD2"]; layer="transcripts",
    x_lims = (x1, x2),y_lims = (y1, y2), alpha=0.5,
    gene_colors=["cyan", "magenta"], 
    height=650, width=600, marker_size = 2,
custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 4.0, adjust_brightness = 0.3)
```
<p float="left">
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/transcripts.png" width=500>
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/transcript_dapi.png" width=500>
</p>

To split them, please use the following codes.
```julia
### DAPI off
cs.sp_feature_plot(breast_cancer, ["ACTA2","TACSTD2"]; layer="transcripts",
    x_lims = (x1, x2),y_lims = (y1, y2), alpha=0.5, combine=false,
    gene_colors=["gray94", "chartreuse"], height=650, width=600, marker_size = 2)
### DAPI on
cs.sp_feature_plot(breast_cancer, ["ACTA2","TACSTD2"]; layer="transcripts",
    x_lims = (x1, x2),y_lims = (y1, y2), alpha=0.5, combine=false,
    height=650, width=600, marker_size = 2, gene_colors=["gray94", "chartreuse"],
    custom_img=true, adjust_coord_to_img="auto", adjust_contrast = 4.0, adjust_brightness = 0.3
)
```
<p float="top">
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/transcript_split2.png" height=400>
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/xenium_he/transcript_split_dapi.png" height=400>
</p>
