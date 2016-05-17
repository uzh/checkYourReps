
require(ezRun)
WORK_DIR <<- "~/R/checkYourReps-data"
ENRICHR_GENELIST_DIRECTORY <<- "/srv/GT/databases/enrichr"
MIN_GENESET_SIZE <<- 5

setwd(WORK_DIR)


setwdNew("check-p1314-EBO6-Muscle-20151125-run1")
datasetFile = "/srv/gstore/projects/p1314/Count_FeatureCounts_8878_2015-11-25--07-44-26/muscle-dataset.tsv"
countDs = EzDataset(file=datasetFile)
replicateGroupingVariable = "Genotype"

param = list()
param[['refBuild']] = countDs$getColumn("refBuild")[1]
param$dataRoot = "/srv/gstore/projects"
param$sigThresh = 50
param$bgExpression = 25
param[['refFeatureFile']] = countDs$getColumn("refFeatureFile")[1]
param[['featureLevel']] = 'gene'
param[['normMethod']] = 'logMean'
param[['expressionName']] = ''
param = ezParam(param)

# countRd = loadCountDataset(countDs, param = param) 
# save(countRd, file="countRd.RData")
load("countRd.RData")
countRd$signal = ezNorm(countRd$counts, method=param$normMethod, presentFlag = countRd$presentFlag)
geneSymbols = countRd$seqAnno$hgnc_symbol


# ########## get gene lists
# geneCats = inverseMapping(goStringsToList(seqAnno$`GO BP`, seqAnno$gene_name))
# geneCats = geneCats[sapply(geneCats, length) >= 5]
# geneCatIdxs = lapply(geneCats, function(gn){match(gn, geneSymbols)})

geneListFiles = list.files(ENRICHR_GENELIST_DIRECTORY, pattern=".gmt$", full.names = TRUE)
names(geneListFiles) = sub(".gmt$", "", basename(geneListFiles))

dataset = countRd$dataset


screenAllSamples(exprCounts, exprConditions, geneListFiles)


## prepare the data ------------
replicateGroups = split(countDs$getNames(), countDs$getColumn(replicateGroupingVariable))
repSamples = replicateGroups[[1]] ## repeat for each set of replicates
countDataUse = selectSamples(countRd, repSamples)
selectedSample = repSamples[1] ## repeat for each sample
controlSamples = setdiff(repSamples, selectedSample)

exprCounts = countDataUse$counts

## get gene set information ready ----
# repeat for each gene list file
gmtFile = geneListFiles[1]
geneSets = readGmt(gmtFile)
geneSets = geneSets[sapply(geneSets, length) >= MIN_GENESET_SIZE]
geneSetIndices = lapply(geneSets, function(gn){na.omit(match(gn, geneSymbols))})


###### analysis using mroast
require(edgeR)
repFactor = repSamples %in% selectedSample
design = model.matrix(~ repFactor)
mroastResult = mroast(log2(exprCounts + param$bgExpression), index=geneSetIndices, design = design)
head(mroastResult[order(mroastResult$PValue.Mixed), ])


#### analysis using camera
require(edgeR)
cameraResult = camera(log2(exprCounts + param$bgExpression), index=geneSetIndices, design = design)
head(cameraResult[order(cameraResult$PValue), ])

## analysis using clusterProfiler
require(clusterProfiler)
## we need to generate a ranked gene list. We rank according to the log2 ratio
x = log2(countRd$signal + param$bgExpression)
scores = x[ , selectedSample] - rowMeans(x[ , controlSamples])
names(scores) = geneSymbols
sortedScores = sort(scores)
term2Gene = ezFrame(term=rep(names(geneSets), sapply(geneSets, length)), gene=unlist(geneSets))

## Fisher's Exact Test:
enricherResult = enricher(names(sortedScores)[1:50], minGSSize = MIN_GENESET_SIZE, TERM2GENE=term2Gene, pvalueCutoff = 1 )
head(summary(enricherResult))

## threshold-free analysis
gseaResult = GSEA(geneList = sortedScores, minGSSize = 10, TERM2GENE=term2Gene[term2Gene$gene %in% names(sortedScores), ], 
                  verbose = TRUE, pvalueCutoff = 1, nPerm = 1000  )
## WARNING this function will use all cores on your system
## you can change this only be overwriting DOSE::gsea
## with nperm > 30 and minGSSize = 3 this tends to fail with the error: Error in sign(ES) : non-numeric argument to mathematical function
## and warnings saying:  In mean.default(x[x >= 0]) : argument is not numeric or logical: returning NA
## apparently a problem of two small gene set size which can cause problems in some permutations
head(summary(gseaResult))


## other approaches:
## - gage
## - segGSEA
## - Category::gseattperm

## different scores as input for GSEA!
## normalization with different gene universe sets!





