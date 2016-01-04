
require(ezRun)
WORK_DIR = "~/R/checkYourReps-data"

setwd(WORK_DIR)



setwdNew("check-p1314-EBO6-Muscle-20151125-run1")
datasetFile = "/srv/gstore/projects/p1314/Count_FeatureCounts_8878_2015-11-25--07-44-26/muscle-dataset.tsv"
countDs = EzDataset(file=datasetFile)
replicateGroupingVariable = "Genotype"

param = list()
param[['refBuild']] = countDs$getColumn("refBuild")[1]
param$dataRoot = "/srv/gstore/projects"
param$sigThresh = 50
param[['refFeatureFile']] = countDs$getColumn("refFeatureFile")[1]
param[['featureLevel']] = 'gene'
param[['normMethod']] = 'logMean'
param[['expressionName']] = ''
param = ezParam(param)

# countRd = loadCountDataset(countDs, param = param) 
# save(countRd, file="countRd.RData")
load("countRd.RData")
countRd$signal = ezNorm(countRd$counts, method=param$normMethod, presentFlag = countRd$presentFlag)
seqAnno = countRd$seqAnno
geneCats = inverseMapping(goStringsToList(seqAnno$`GO BP`, seqAnno$gene_name))
geneCats = geneCats[sapply(geneCats, length) >= 5]
geneCatIdxs = lapply(geneCats, function(gn){match(gn, seqAnno$gene_name)})


replicateGroups = split(countDs$getNames(), countDs$getColumn(replicateGroupingVariable))

require(edgeR)

repSamples = replicateGroups[[1]] ## repeat for each set of replicates
countDataUse = selectSamples(countRd, repSamples)

selectedSample = repSamples[1] ## repeat for each sample


###### analysis using mroast
repFactor = repSamples %in% selectedSample
design = model.matrix(~ repFactor)

mroastResult = mroast(countDataUse$counts, index=geneCatIdxs, design = design)
head(mroastResult[order(mroastResult$PValue.Mixed), ])


#### analysis using camera
cameraResult = camera(countDataUse$counts, index=geneCatIdxs, design = design)
head(cameraResult[order(cameraResult$PValue), ])


### analysis using GSEA







#goana.default


