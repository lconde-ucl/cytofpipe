library(citrus)
library(hash)
library(ini)


options(stringsAsFactors = F)
rm(list = ls())

#------------------------------------------------------------------
#- Parse arguments
#------------------------------------------------------------------


jobid <- as.character(Sys.getenv("RAND_ID"))
input <- paste0(jobid, ".txt")

args<-read.ini(input)

dataDirectory=args$paramscitrus$INPUTFILE
conditionsFile=args$paramscitrus$CONDITIONS
outputDirectory=args$paramscitrus$OUTPUTFILE
markersFile=args$paramscitrus$MARKERSFILE
asinh_cofactor=args$paramscitrus$ASINH
mergeMethod = args$paramscitrus$MERGE
fileSampleSize = args$paramscitrus$DOWNSAMPLE
medians = args$paramscitrus$MEDIANS

if(asinh_cofactor == '-'){
	asinh_cofactor=5;
}


#——————————————————————--------
#- READ FCS AND EXTRACT MARKERS
#——————————————————————--------

files <- list.files(dataDirectory,pattern='.fcs$', full=TRUE)
usermarkers <- as.character(read.table(markersFile, header = FALSE)[,1])
 
fcs1<-read.FCS(files[1])

allMarkerNames<-as.vector(pData(parameters(fcs1))$name)
allMarkerDesc<-as.vector(pData(parameters(fcs1))$desc)
 
UserName2Desc <-hash()
for(i in 1:length(allMarkerDesc)){
	if(is.na(allMarkerDesc[i])){
 		UserName2Desc[[ allMarkerNames[i] ]] <- allMarkerNames[i]
	}else{
		UserName2Desc[[ allMarkerDesc[i] ]] <- allMarkerDesc[i]
	}
}
if (sum(has.key( usermarkers, UserName2Desc )) == 0) {
  	clear(UserName2Desc)
 	clear(Desc2UserName)
  	for(i in 1:length(allMarkerDesc)){
  		if(!is.na(allMarkerDesc[i])){
 			id <- gsub( "^[^_]+_", "", allMarkerDesc[i])
  			UserName2Desc[[ id ]] <- allMarkerDesc[i]
  			Desc2UserName[[ allMarkerNames[i] ]] <- id
		}
  	}
}
Desc2Name <-hash()
for(i in 1:length(allMarkerDesc)){
 	if(is.na(allMarkerDesc[i])){
 		Desc2Name[[ allMarkerNames[i] ]] <- allMarkerNames[i]
	}else{
		Desc2Name[[ allMarkerDesc[i] ]] <- allMarkerNames[i]
	}
}
  
markersDesc <-vector()
for(i in 1:length(usermarkers)){
	markersDesc[i]<-values(UserName2Desc,  keys=usermarkers[i])
}
markersName <-vector()
for(i in 1:length(markersDesc)){
	markersName[i]<-values(Desc2Name, keys=markersDesc[i])
}


#———————————----———
#- PARAMETERS
#——————————————---


## @knitr fixedparameters

family = "classification"
minimumClusterSizePercent = 0.05
nFolds = 1
featureType = ""

if(mergeMethod == '-'){
	if(fileSampleSize == '-'){
		fileSampleSize = 10000
	}else{
		fileSampleSize = as.numeric(fileSampleSize)
	}	
}else{
	fileSampleSize="NULL";
}

## @knitr medians

mediansmarkersName<-vector()
if(basename(medians) == '-'){
	featureType=c("abundances")
}else{
	featureType=c("medians")

	mediansmarkers <- as.character(read.table(medians, header = FALSE)[,1])
 
	mediansmarkersDesc <-vector()
	for(i in 1:length(mediansmarkers)){
		mediansmarkersDesc[i]<-values(UserName2Desc,  keys=mediansmarkers[i])
	}
	mediansmarkersName <-vector()
	for(i in 1:length(mediansmarkersDesc)){
		mediansmarkersName[i]<-values(Desc2Name, keys=mediansmarkersDesc[i])
	}

}

## @knitr parameters

clusteringColumns<-markersName
transformColumns<-allMarkerNames[-c(grep("Time|Event|Cell_length|viability", allMarkerNames, ignore.case = TRUE))]
scaleColumns = transformColumns
transformCofactor <- as.numeric(asinh_cofactor) 
medianColumns<-mediansmarkersName

data_conditions<-read.table(conditionsFile, sep="\t",header=F)
fileList = data.frame(defaultCondition=data_conditions[,1])
labels = as.factor(data_conditions[,2])

## @knitr modelTypes

modelTypes<-vector()
if(length(levels(labels)) > 2){
	modelTypes = c("pamr","sam")
}else{
	modelTypes = c("pamr","glmnet","sam")
}



#———————————----———
#- CITRUS
#——————————————---

## @knitr citrus

# Read Data
citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList,fileSampleSize,transformColumns,transformCofactor)

# Cluster all the data
citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels,nFolds)

# Make vector of conditions for analysis. If comparing two conditions, should be 
# two elements - first element is baseline condition and second is comparison condition.
conditions = colnames(fileList)[1]

# Build cluster features
citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering,citrus.combinedFCSSet,
                                                         featureType=featureType,
                                                         minimumClusterSizePercent=minimumClusterSizePercent,
                                                         conditions=conditions,
                                                         medianColumns=medianColumns
                                                         )

# Build regression models for each model type
citrus.regressionResults = mclapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family)

# Plot Results for each model
lapply(citrus.regressionResults,plot,outputDirectory=outputDirectory,citrus.foldClustering=citrus.foldClustering,citrus.foldFeatureSet=citrus.foldFeatureSet,citrus.combinedFCSSet=citrus.combinedFCSSet,
	theme="white")


#-- export largeenoughclusters
allLargeEnoughClusters<-citrus.foldFeatureSet$allLargeEnoughClusters 
dir.create(file.path(outputDirectory, "exportedClusters"))
for (i in 1:length(allLargeEnoughClusters)){
	dir.create(file.path(paste0(outputDirectory,"/exportedClusters"), allLargeEnoughClusters[i]))
	citrus.exportCluster(allLargeEnoughClusters[i],citrus.clustering=citrus.foldClustering$allClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,outputDirectory=paste0(outputDirectory,"/exportedClusters/",allLargeEnoughClusters[i]))
}

#-- Save object
save(citrus.combinedFCSSet,citrus.foldFeatureSet, file = paste0(outputDirectory, "/cytofpipe_citrusClustering.rData"))


paste0("files: ", files)
paste0("data_conditions: ", data_conditions)
paste0("markersName: ", markersName)
if(mergeMethod == '-'){
	paste0("fileSampleSize: ", as.numeric(fileSampleSize))
}else{
	paste0("mergeMethod: ", mergeMethod)
}
paste0("asinh_cofactor: ", as.numeric(asinh_cofactor))
paste0("family: ", family)
paste0("minimumClusterSizePercent: ", minimumClusterSizePercent)
paste0("nFolds: ", nFolds)
paste0("featureType: ", featureType)
paste0("modelTypes: ", modelTypes)


sessionInfo()

