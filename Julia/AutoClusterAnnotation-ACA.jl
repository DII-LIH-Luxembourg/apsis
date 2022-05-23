## ----------------------------------------------------------------------------
##
##                  ACA - Automated Cluster Annotation
## Annotate each cluster based on their expression profile
## defined by PBMC_Cell_Matrix.csv and marker_thresholds.csv
##
## ----------------------------------------------------------------------------
using Pkg;
Pkg.activate(".")
# Pkg.instantiate()
using CSV
using DataFrames
using Plots
using StatsBase
using XLSX
include("ACAhelper.jl")
using .ACAhelper

# -----------------------------------------------------------------------------
# Main settings
# -----------------------------------------------------------------------------
maindir = "" # set the path
datadir = ""
outdir = datadir

clusterSize = 1024

exprFile = "median_expr_$clusterSize.csv"
freqFile = "cluster_freq_$clusterSize.csv"
cMatrixFile = "PBMC_Cell_Matrix.csv"
thFile = "marker_thresholds.csv"

annotatedDF = ACAhelper.doAnnotation(maindir, datadir, outdir, exprFile, freqFile , thFile, cMatrixFile)

CSV.write(outdir*"/annotated_df_$clusterSize.csv", annotatedDF)

# create the mean table 
tmp = annotatedDF[:, Not(:cluster)]
groupby(tmp, :cell)
X = combine(groupby(tmp, :cell), 1:35 .=> mean, :size => sum)

CSV.write(outdir*"/mean_marker_expr_df.csv", X)
# output the cluster merging file_id
CSV.write(outdir*"/cluster_merging.csv" ,annotatedDF[!, [:cell, :cluster]] |> DataFrame)
CSV.write(outdir*"/som_df_norm_$clusterSize.csv", annotatedDF)











