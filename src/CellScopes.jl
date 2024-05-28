__precompile__()

module CellScopes

function __init__()
    println("Welcome to use CellScopes.jl!")
end

using Pkg
using TOML
using MatrixMarket
using DataFrames
using CSV
using HTTP
using JSON
using ProgressMeter
using Statistics
using StatsModels
using MultivariateStats
using Plots
using UMAP
using NearestNeighborDescent
using Distances
using SparseArrays
using Random
using Leiden
using StatsBase
using VegaLite
using Loess
using JLD2
using TSne
using HypothesisTests
using MultipleTesting
using Match
using Colors
using ColorSchemes
import CairoMakie as MK
using StatsPlots
using Folds
using PyCall
using RCall
import PlotlyJS as plyjs
using Images
using Polynomials
using FileIO
using CSVFiles
using GZip
using Augmentor
using Grep
using LinearAlgebra
using Glob
using RollingFunctions
using DataFramesMeta
using GeneticsMakie
using Arrow
using HDF5
using PyCallUtils
import FFTW
import ImageMorphology
import Graphs
using NearestNeighbors: KDTree, knn
using SimpleWeightedGraphs
using StatsBase: countmap
using KernelDensity
using Graphs: src, dst
using DuckDB
using MLJMultivariateStatsInterface: PCA as MLJ_PCA
using MLJ
using Clustering: kmeans
using Logging
using CategoricalArrays

include("scrna/objects.jl")
include("spatial/sp_objects.jl")
include("scatac/atac_objects.jl")
include("spatial_highres_full/sp_hf_objects.jl")
include("integration/int_objects.jl")
include("properties.jl")
include("spatial_highres_full/sp_hf_utils.jl")
include("scrna/processing.jl")
include("scrna/utils.jl")
include("scrna/visualization.jl")
include("fileio.jl")
include("spatial/sp_plots.jl")
include("spatial/sp_utils.jl")
include("spatial/sp_analysis.jl")
include("spatial/sp_processing.jl")
include("scatac/atac_processing.jl")
include("scatac/atac_plots.jl")
include("scatac/atac_utils.jl")
include("spatial/baysor_boundary.jl")
include("spatial_highres_full/sp_hf_plots.jl")
include("spatial_highres_full/sp_hf_processing.jl")
include("integration/int_utils.jl")
include("integration/int_processing.jl")

end
