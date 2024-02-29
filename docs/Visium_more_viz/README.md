# Use a high-resolution histology image for Visium data visualization
SpaceRanger generates output files that include both high and low-resolution H&E images. These images are very helpful for data intepretation because they can be used to overlay the cell annotations or gene expression on the Visium spots. For overall views of the whole tissue, the high-resolution images usually offer enough detail. However, when attempting to zoom in to closely inspect more intricate structures, the default images often lack the necessary resolution. This tutorial guides you through the steps of adding a high-resolution H&E image captured with a confocal microscope with the Visium spots using ```CellScopes.jl```. 

### 1. Input file: a full-resolution H&E image
Since SpaceRanger outputs a high-resolution tiff image, to avoid confusion on the names, we refer to the H&E image captured by the confocal microscope as a "full-resolution" image. This full-resolution H&E image must be aligned with the fiducial marker locations on the Visium slide. This can be done in Loupe Browser. Please refer to the guides [Manual Fiducial Alignment for Visium](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/image-fiducial-alignment) provided by 10x Genomics. 

### 2. Visium dataset
For this tutorial, we will be analyzing a breast cancer Visium dataset from 10x Genomics. https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast.
We start by reading the data into Julia, creating a CellScopes VisiumObject, and further processing the object for data normalization, dimensional reduction and clustering. For more details, please refer to https://github.com/HaojiaWu/CellScopes.jl/tree/main/docs/visium_tutorial
```julia
visium = cs.read_visium("/mnt/sdb/breast_cancer_visium/")
visium = cs.normalize_object(visium; scale_factor = 10000)
visium = cs.scale_object(visium)
visium = cs.find_variable_genes(visium; nFeatures = 1000)
visium = cs.run_pca(visium;  method=:svd, pratio = 1, maxoutdim = 10)
visium = cs.run_umap(visium; min_dist=0.2, n_neighbors=20)
visium = cs.run_clustering(visium; res=0.01, n_neighbors=20)
```

### 3. Add the aligned full-resolution image into a Visium CellScopes object
The full resolution H&E tiff file can be downloaded from the 10x website. 
```sh
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_tissue_image.tif
```
This tiff image is in 16-bit format. It can be converted into 8-bit format with Fiji to reduce the file size. Since the image has been aligned to the visium fiducials, we can directly add it to the VisiumObject using the CellScopes function ```add_visium_img```.
```julia
visium = cs.add_visium_img(visium, 
        "/mnt/sdb/breast_cancer_visium/CytAssist_FFPE_Human_Breast_Cancer_tissue_image_8bit.tif")
```

### 4. Visualize gene expression on the full-res H&E image
Now the Visium object in CellScopes should have three H&E images: low res, high res and full res. They are stored in the ```imageData``` slot. The low and high resolution images were from the Space Ranger output. For visualizing genes in the entire organ, we recommend to use the low resolution image (by setting ```img_res = "low"```) since this can save time and there is minimal noticeable difference to the naked eye between low and full resolution when showing genes in a large view.
```julia
cs.sp_feature_plot(visium, ["CEACAM6"]; 
    marker_size = 7, color_keys=["azure1", "lightsteelblue1" ,"blue"], 
    adjust_contrast=1, adjust_brightness = 0.0, scale=true, alpha=[0,0.6],clip=0.5,
    height=500, width=800,img_res="low"
)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/visium_he/gene_whole.png" width="600"> <br>

When we focus on a specific area of the field, the image at full resolution provides clearer histological details of the tissue. Here, we have placed low, high, and full resolution images side by side to highlight the differences.
#### 4a. Gene expression on low-res image
The low-res image is from the Space Ranger output and has been incorporated in the VisiumObject by default (in the ```read_visium``` step).
```julia
cs.sp_feature_plot(visium, ["KRT17","TACSTD2"]; 
    marker_size = 50, color_keys=["azure1", "lightsteelblue1" ,"blue"], 
    adjust_contrast=1, adjust_brightness = 0.0, scale=true, alpha=[0,0.6],clip=0.6,
    height=600, width=400,img_res="low", x_lims = (7119, 8246), y_lims=(9539, 11828)
)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/visium_he/gene_low.png" width="600"> <br>
#### 4b. Gene expression on high-res image
Similar to the low-res image, the high-res image is also included in the VisiumObject by default in the ```read_visium``` step. The figure below demonstrates that even the high-resolution image offered by Visium does not adequately reveal the tissue structure. I
```julia
cs.sp_feature_plot(visium, ["KRT17","TACSTD2"]; 
    marker_size = 50, color_keys=["azure1", "lightsteelblue1" ,"blue"], 
    adjust_contrast=1, adjust_brightness = 0.0, scale=true, alpha=[0,0.6],clip=0.6,
    height=600, width=400,img_res="high", x_lims = (7119, 8246), y_lims=(9539, 11828)
)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/visium_he/gene_high.png" width="600"> <br>
#### 4c. Gene expression on full-res image
In contrast, the full-res image captured by confocal delivers a significantly clearer resolution.
```julia
cs.sp_feature_plot(visium, ["KRT17","TACSTD2"]; 
    marker_size = 50, color_keys=["azure1", "lightsteelblue1" ,"blue"], 
    adjust_contrast=1, adjust_brightness = 0.0, scale=true, alpha=[0,0.6], clip=0.6,
    height=600, width=400,img_res="full", x_lims = (7119, 8246), y_lims=(9539, 11828)
)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/visium_he/gene.png" width="600"> <br>
In the ```sp_feature_plot``` function provided above, there are several parameters that could be changed to create a figure that is easier to interpret. The ```img_res``` parameter specifies the resolution of the images to be used. The ```adjust_contrast``` and ```adjust_brightness``` parameters allow for the adjustment of the contrast and brightness of the H&E image, respectively. The ```alpha``` parameter is used to adjust the transparency of the spots to allow the underlying tissue structure to be seen. The ```clip``` parameter (0,1) is designed to hide spots below a certain threshold, and show only those spots whose expression is above the clip threshold. Users can play around those parameters to get a better output image.

### 5. Visualize cell annotation on the full-res H&E image
Cell type annotation can also be plotted on the full-res H&E image too. Here is the way to do it:
```julia
cs.sp_dim_plot(visium, "cluster"; 
    marker_size = 50, height=600, width=400, img_res="full", alpha=0.5,
    do_label=false, adjust_contrast=1, adjust_brightness = 0, 
    x_lims = (7119, 8246), y_lims=(9539, 11828))
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/visium_he/cell_he.png" width="400"> <br>
If you want to highlight only a particular cell type, just feed the cell type name into the ```cell_highlight``` parameter. For example:
```julia
cs.sp_dim_plot(visium, "cluster"; 
    marker_size = 50, height=600, width=400, img_res="full", cell_highlight="10",
    do_label=false, adjust_contrast=1, adjust_brightness = 0, alpha=0.5,
    x_lims = (7119, 8246), y_lims=(9539, 11828))
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/visium_he/cell_he_highlight.png" width="400"> <br>

### 6. Select the region of interest to visualize
Now you might be curious about how to easily obtain the ```x_lims``` and ```y_lims``` parameters to zoom into your area of interest for a closer lookup. The ```plot_fov``` and ```subset_fov``` functions we developped are designed to assist you in achieving this.
```julia
cs.plot_fov(visium, 10, 10)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/visium_he/fov_he.jpg" width="600"> <br>
You can also highlight a particular cell type in the background to help you better locate the field of view (fov). Again, just feed the cell type name to the ```cell_highlight``` parameter.
```julia
cs.plot_fov(visium, 10, 10; group_label="cluster", cell_highlight="10", alpha=0.7)
```
<img src="https://github.com/HaojiaWu/CellScopes.jl/blob/main/data/visium_he/fov_highlight.jpg" width="600"> <br>
Then simply select the numbers, or the numbers at the 4 corners (in numerical order) of a square or rectangle to crop your fov of interest. Here is an example how we get the x_lims and y_lims to visualize the structure provided in this tutorial.
```julia
df = cs.subset_fov(visium, [26,27], 10, 10)
x1 = minimum(df.x)
x2 = maximum(df.x)
y1 = minimum(df.y)
y2 = maximum(df.y)
x_lims = (x1, x2)
y_lims = (y1, y2)
```


