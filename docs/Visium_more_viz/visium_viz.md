# Use a high-resolution histology image for Visium data visualization
SpaceRanger generates output files that include both high and low-resolution H&E images. These images are very useful for data intepretation because they can be used to overlay the cell annotations or gene expression on the Visium spots. For overall views of the whole tissue, the high-resolution images usually offer enough detail. However, when attempting to zoom in to closely inspect more intricate structures, the default images often lack the necessary resolution. This tutorial guides you through the steps of adding a high-resolution H&E image captured with a confocal microscope with the Visium spots using ```CellScopes.jl```. 

### 1. Input file: a full-resolution H&E image
Since SpaceRanger output a high-resolution tiff image, to avoid confusion on the names, we refer to the H&E image captured by the confocal microscope as a "full-resolution" image. This full-resolution H&E image must be aligned with the fiducial marker locations on the Visium slide. This can be done in Loupe Browser. Please refer to the guides [Manual Fiducial Alignment for Visium](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/image-fiducial-alignment) provided by 10x Genomics. 

### 2. Add the aligned full-resolution image into a Visium CellScopes object



