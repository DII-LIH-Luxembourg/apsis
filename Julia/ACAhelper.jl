module gigaSomHelper
using CSV
using DataFrames
using Plots
using StatsBase
using XLSX


struct CellType
	
	posMarker::Array{String,1}
	negMarker::Array{String,1}
	name::String
	
end


"""
    plotClusterMarkerHist(datadir, outdir, filename)

Plot histogram of marker expression and save it to .png file.
"""
function plotClusterMarkerHist(datadir, outdir, filename)

    exprDF = CSV.File(datadir * "/" * filename) |> DataFrame

    lineage_markers = names(exprDF)
    ds = DataFrame()
    # plot density per condition
    for (i, m) in enumerate(lineage_markers)
        # get the marker to display
            # ds[!, :marker] = scaleNorm(exprDF[:, i])
            # ds[!, :marker] = exprDF[:, i]
        a = histogram(scaleNorm(exprDF[:, i]), bins=200, xticks=0:0.05:10)
        savefig(a, outdir * "/density_plot_$m.png")    
        
    end

end

"""
    getNormDF(datadir, filename)

Scale each marker expression bewteen 0 and 1.
"""
function getNormDF(datadir, filename)

    exprDF = CSV.File(datadir * "/" * filename) |> DataFrame
    m = exprDF |> Matrix
    dt = fit(UnitRangeTransform, m, dims=1)
    m = StatsBase.transform(dt, m)

    normDF = m |> DataFrame
    rename!(normDF, names(exprDF)) 

    return normDF
end

"""
    getCellDefinitions(maindir, cMatrixFile)

Create a vector of cell type specific objects with their negative and positive
Expression markers. 
"""
function getCellDefinitions(maindir, cMatrixFile)

    cMatrix = CSV.File(maindir * "/" * cMatrixFile) |> DataFrame
    allCells = Vector{gigaSomHelper.CellType}()

    for i in 1:size(cMatrix, 1)
        posExpr = names(cMatrix)[Array(cMatrix[i,:]) .== 1]
        negExpr = names(cMatrix)[Array(cMatrix[i,:]) .== -1]
        newCell = CellType(posExpr, negExpr, cMatrix[i, :cell])
        push!(allCells, newCell)
    end
    
    return allCells
end

function annotateRow(df, rownumber, cell, cellName)

    if df[rownumber, cell] == "notDefined"
        df[rownumber, cell] = cellName
    else
        df[rownumber, cell] = df[rownumber, cell] * ',' * cellName
    end
end

"""
    doAnnotation(maindir, datadir, outdir, exprFile, freqFile, thFile, cMatrixFile)

Annotates each cluster by the matching expression values defined in cell definitions.
"""
function doAnnotation(maindir, datadir, outdir, exprFile, freqFile, thFile, cMatrixFile)

    # Read the cluster expression table
    normDF = gigaSomHelper.getNormDF(datadir, exprFile)
    thresholds = CSV.File(maindir * "/" * thFile) |> DataFrame
    cluster_freq = CSV.File(datadir * "/" * freqFile) |> DataFrame

    clustering_prop = round.((cluster_freq.column_name ./ sum(cluster_freq.column_name) * 100), digits=3)

    normDF[!, :size] = clustering_prop
    normDF[!, :cell] .= "notDefined"
    normDF[!, :cluster] .= 0

    allCells = gigaSomHelper.getCellDefinitions(maindir, cMatrixFile)

    for ct in allCells
        for i in 1:size(normDF, 1)
            mNeg = normDF[i, ct.negMarker]
            mPos = normDF[i, ct.posMarker]

            # if we have pos and neg marker defined
            if !isempty(mNeg) && !isempty(mPos)
                thNeg = thresholds[∈(ct.negMarker).(thresholds.cell), :value]
                thPos = thresholds[∈(ct.posMarker).(thresholds.cell), :value]
        
                if all(Array(mNeg) .<= thNeg) && all(Array(mPos) .> thPos)
                    gigaSomHelper.annotateRow(normDF, i, :cell, ct.name)
                end
            # only pos markers defined
            elseif isempty(mNeg)
                thPos = thresholds[∈(ct.posMarker).(thresholds.cell), :value]
        
                if all(Array(mPos) .> thPos)
                    gigaSomHelper.annotateRow(normDF, i, :cell, ct.name)
                end
            # only neg markers defined
            elseif isempty(mPos)
                thNeg = thresholds[∈(ct.negMarker).(thresholds.cell), :value]
        
                if all(Array(mNeg) .<= thNeg)
                    gigaSomHelper.annotateRow(normDF, i, :cell, ct.name)
                end
            end
        end
    end    

    return normDF
end

end

