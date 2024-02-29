# Use a high-resolution histology/IF image for Xesium data visualization
```CellScopes.jl``` provides functions that allow for the integration of an external histology or immunofluorescence image for visualizing Xenium data. Given that Xenium data provide three layers of spatial information (transcripts, genes, and cells), this tutorial will walk you through how to visualize this information on top of a high-resolution histology or immunofluorescence (IF) image.
### 1. Input image: a registered H&E or IF image
This tutorial requires a pre-registered H&E/IF image that has been accurately scaled and aligned with the Xenium transcript coordinates. The process of image registration can be quite challenging. For H&E or IF images derived from the identical section that has undergone Xenium processing, linear transformation (affine) is generally sufficient. However, when attempting to register images from consecutive tissue sections, more advanced algorithms, including those based on deep learning, are often required to achieve accurate alignment. This tutorial will not cover image registration. Users are encouraged to explore the existing methods to achieve a well-registered image before beginning this tutorial. For demo purpose, we will use the pre-registered H&E and IF image provided by the Xenium dataset from 10x website.
### 2. Xenium dataset
For this tutorial, we will analyze a breast cancer Xenium dataset from 10x Genomics. The dataset can be downloaded from here: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast. 
```julia
import CellScopes as cs
breast_cancer = read_xenium("/mnt/sdb/breast_cancer_xenium/")
```
This dataset come with a H&E image taken on that same Xenium slide. Xenium analyzer also output a DAPI tiff image that has been aligned to the xenium transcript coordinates. We will use these two images to showcase the ploting functions provided by CellScopes.
