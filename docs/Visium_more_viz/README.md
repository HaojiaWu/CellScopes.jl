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
Now the Visium object in CellScopes should have three H&E images: low res, high res and full res. They are stored in the ```imageData``` slot. The low and high resolution images were from the Space Ranger output. For visualizing genes in the entire organ, we recommend to use the low resolution image (by setting ```img_res = "low"```) since this can save time and there is minimal noticeable difference to the naked eye between low and full resolution when showing genes in large view.
```julia
cs.sp_feature_plot(visium, ["CEACAM6"]; 
    marker_size = 7, color_keys=["azure1", "lightsteelblue1" ,"blue"], 
    adjust_contrast=1, adjust_brightness = 0.0, scale=true, alpha=[0,0.6],clip=0.5,
    height=500, width=800,img_res="low"
)
```



