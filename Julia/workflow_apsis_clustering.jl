using GigaSOM, Distributed
import NearestNeighbors
using Random
using Clustering
import Distances
using DataFrames
using CSV
using StatsBase
using Statistics
using DataFrames

# add the path
datapath = ""
cd(datapath)

# som grid size
gridSize = 32
# training iterations
nEpochs = 40

# adjust the number of workers for parallel processing
addprocs(20)
@everywhere using GigaSOM
using XLSX
md = GigaSOM.DataFrame(XLSX.readtable("metadata.xlsx", "Sheet1", infer_eltypes=true)...)
panel = GigaSOM.DataFrame(XLSX.readtable("panel.xlsx", "Sheet1", infer_eltypes=true)...)

_, fcsParams = loadFCSHeader(md[1, :file_name])

_, fcsAntigens = getMarkerNames(getMetaData(fcsParams))

antigens = panel[panel[:,:Lineage].==1, :Antigen]

# cleanNames functions to check for unwanted characters in column names
cleanNames!(antigens)
cleanNames!(fcsAntigens)
lineage_markers = panel.Antigen[panel.Lineage.==1]
cleanNames!(lineage_markers)

di = loadFCSSet(:fcsData, md[:,:file_name]; postLoad=selectFCSColumns(lineage_markers))

cols = Vector(1:length(lineage_markers)) # shortcut for "all rows"
dtransform_asinh(di, cols, 5)

Random.seed!(1)

som = initGigaSOM(di, gridSize, gridSize, seed=13) # set a seed value
som = trainGigaSOM(som, di, epochs = nEpochs)

e = embedGigaSOM(som, di)
e = distributed_collect(e) 

somClusters = mapToGigaSOM(som, di)
mapping = distributed_collect(mapToGigaSOM(som, di), free=true)

CSV.write("mapping.csv", DataFrame(hcat(mapping, 1:length(mapping)), :auto))

sCodes = DataFrame(som.codes)
rename!(sCodes, Symbol.(lineage_markers))
CSV.write("somCodes.csv", sCodes)

StatsBase.counts(mapping)

# avoid any meta-clustering by seeting the value to its actual grid size
mc = 1024 

metaClusters =
cutree(k=mc,
        hclust(linkage=:average,
            GigaSOM.distMatrix(Distances.Euclidean())(som.codes)))

distributed_transform(mapping, (m)->metaClusters[m])
clusterFreq = dcount(mc, mapping)
df = DataFrame(column_name = clusterFreq)
CSV.write("cluster_freq_$mc.csv", df)

files=distributeFCSFileVector(:fileIDs, md[:,:file_name])

# Get the count table per fileID
count_tbl = dcount_buckets(mc, mapping, size(md,1), files)
ct = convert(DataFrame, count_tbl)
rename!(ct, md.sample_id)
# export the count talbe
CSV.write("cluster_counts_$mc.csv", ct)

# Get the median expression per cluster
expr_tbl = dmedian_buckets(di, mc, mapping, cols)
# expr_tbl = dstat_buckets
et = convert(DataFrame, expr_tbl)
rename!(et, lineage_markers)

# export median marker expression
CSV.write("median_expr_$mc.csv", et)

# now try to select embed coordinates for specific condition
distributed_collect(files)[distributed_collect(files) .== 1]
# attach the file index to every x y value
eCoord = DataFrame(e)
eCoord[:ind] = distributed_collect(files)
eCoord[:meta] = distributed_collect(mapping)

# write embed to csv
CSV.write("embed_XY_$mc.csv", eCoord)
