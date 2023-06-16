__precompile__()

module CellScopes

function __init__()
    println("Welcome to use CellScopes.jl!")
end

using Pkg
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

include("scrna/objects.jl")
include("spatial/sp_objects.jl")
include("scatac/atac_objects.jl")
include("properties.jl")
include("scrna/processing.jl")
include("scrna/utils.jl")
include("scrna/visualization.jl")
include("scrna/fileio.jl")
include("spatial/sp_plots.jl")
include("spatial/sp_utils.jl")
include("spatial/sp_analysis.jl")
include("spatial/sp_processing.jl")
include("scatac/atac_processing.jl")
include("scatac/atac_plots.jl")


end
