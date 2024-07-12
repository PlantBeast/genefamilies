## tutorial
# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
## phew

library("topGO")
library("Rgraphviz")

#do first with OGs unique to Iochroma from Orthogroup analysis

geneID2GO <- readMappings(file = "~/Downloads/TopGO Analyses/IC_OG_w_GOs.tsv")   
names(geneID2GO)->geneUniverse
genesOfInterest <- read.table("~/Downloads/ICUniqueOGs.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="IC unique OGs", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO) ## should repeat for BP and CC
resultFisher<-runTest(myGOdata, algorithm="classic", statistic="fisher")
allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

#redo w/expanded and contracted set from CAFE

geneID2GO <- readMappings(file = "~/Downloads/IC_OG_w_GOs.tsv")   
names(geneID2GO)->geneUniverse
genesOfInterest <- read.table("~/Downloads/TopGO Analyses/expanded_in_IC.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="Expanded IC OGs", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO) ## should repeat for BP and CC
resultFisher<-runTest(myGOdata, algorithm="classic", statistic="fisher")
allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 100)
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)



