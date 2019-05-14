## @knitr libraries

library(cytofkit) 
library(flowCore)
library(ini)
library(hash)
library(openCyto)
library(mvtnorm)

require(reshape2)
require(VGAM)
require(colourpicker)
require(gplots)

#------------------------------------------------------------------
#- Parse parameters
#------------------------------------------------------------------

## @knitr parameters

jobid <- as.character(Sys.getenv("JOB"))
input <- paste0(jobid, ".txt")

args<-read.ini(input)

inputfiles=args$paramsclustering$INPUTFILE
outputdir=args$paramsclustering$OUTPUTFILE
markersFile=args$paramsclustering$MARKERSFILE
configFile =args$paramsclustering$CONFIGFILE
template=args$paramsclustering$GATINGFILE
transformMethod = args$paramsclustering$TRANSFORM
mergeMethod = args$paramsclustering$MERGE
fixedNum = args$paramsclustering$DOWNSAMPLE
displayAll = args$paramsclustering$DISPLAY_ALL
groupfile = args$paramsclustering$GROUPS
randomSampleSeed = args$paramsclustering$RANDOM_SAMPLE_SEED
randomTsneSeed = args$paramsclustering$RANDOM_TSNE_SEED
randomFlowSeed = args$paramsclustering$RANDOM_FLOW_SEED
array = args$paramsclustering$ARRAY


#---------------------------------------------------------------------------------------------
#- Functions
#---------------------------------------------------------------------------------------------

## @knitr functions_opencyto


#- A custom gating function for a DNA/DNA gate on CyTOF data.
#- Finds the intersection between a quantile of a multivariate normal fit
#- of a population and a boundary along y = -x+b 
#- author: jfreling@fhcrc.org
boundry <-  function(xs) {
    # find the boundry events that are above a quantile and below a line 

    cxs <- scale(xs) # scale data so that it can be compaired to the results from qnorm
    f <- qnorm(0.95) # set a boundry level
    pd <- dmvnorm(c(f, f))[1] # and find the p(x) for that level

    pxs <- dmvnorm(x=cxs)  
    idxs <- (pxs > pd) # find those points who are above the boundy level

    idxs2 <- ((-1*cxs[,1]) + 1.96) > cxs[,2] # find points that are below the line with y=-1*x+b 
    pos_xs <- xs[idxs&idxs2,] # intersection of points below line and above threshold level

    hpts <- chull(pos_xs) # find the boundry points of the intersection of cells
    return(pos_xs[hpts,])
}

.dnaGate <- function(fr, pp_res, channels = NA, filterId="", ...){
   xs <- exprs(fr[,channels])
  pnts <- boundry(xs)
  return(polygonGate(.gate=pnts, filterId=filterId))
}

registerPlugins(fun=.dnaGate,methodName='dnaGate', dep='mvtnorm','gating')


## @knitr functions_cytofkit

#- Function to plot all level plots for all markers (https://github.com/JinmiaoChenLab/cytofkit/blob/master/inst/shiny/global.R)
cytof_wrap_colorPlot <- function(data, xlim=NULL, ylim=NULL, xlab, ylab, markers, scaleMarker = FALSE,
                             colorPalette = c("bluered", "spectral1", "spectral2", "heat"), 
                             pointSize=1, min=min, max=max,
                             removeOutlier = TRUE){
     
     remove_outliers <- function(x, na.rm = TRUE, ...) {
         qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
         H <- 1.5 * IQR(x, na.rm = na.rm)
         y <- x
         y[x < (qnt[1] - H)] <- qnt[1] - H
         y[x > (qnt[2] + H)] <- qnt[2] + H
         y
     }
     
     data <- as.data.frame(data)
     title <- "Marker Expression Level Plot"
     data <- data[,c(xlab, ylab, markers)]
     
     if(removeOutlier){
         for(m in markers){
             data[[m]] <- remove_outliers(data[ ,m])
         }
     }
     
     if(scaleMarker){
         data[ ,markers] <- scale(data[ ,markers], center = TRUE, scale = TRUE)
         ev <- "ScaledExpression"
         data <- melt(data, id.vars = c(xlab, ylab), 
                      measure.vars = markers,
                      variable.name = "markers", 
                      value.name = ev)
     }else{
         ev <- "Expression"
         data <- melt(data, id.vars = c(xlab, ylab), 
                      measure.vars = markers,
                      variable.name = "markers", 
                      value.name = ev)
     }
     
 
     colorPalette <- match.arg(colorPalette)
     switch(colorPalette,
            bluered = {
                myPalette <- colorRampPalette(c("blue", "white", "red"))
            },
            spectral1 = {
                myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                                "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                                "#F46D43", "#D53E4F", "#9E0142"))
            },
            spectral2 = {
                myPalette <- colorRampPalette(rev(c("#7F0000","red","#FF7F00","yellow","white", 
                                                    "cyan", "#007FFF", "blue","#00007F")))
            },
            heat = {
                myPalette <- colorRampPalette(heat.colors(50))
            }
     )
     zlength <- nrow(data)
     grid_row_num <- round(sqrt(length(markers)))

     if(!is.null(xlim)){
	     gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = ev)) + 
	         facet_wrap(~markers, nrow = grid_row_num, scales = "fixed") +
	         scale_colour_gradientn(name = ev, colours = myPalette(zlength), limits=c(min, max)) +
	         geom_point(size = pointSize) + theme_bw() + coord_fixed() +
                 xlim(xlim) +
                 ylim(ylim) +
	         theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
	         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	         theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
     }else{
	     gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = ev)) + 
	         facet_wrap(~markers, nrow = grid_row_num, scales = "fixed") +
	         scale_colour_gradientn(name = ev, colours = myPalette(zlength)) +
	         geom_point(size = pointSize) + theme_bw() + coord_fixed() +
	         theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
	         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	         theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
     }     
     return(gp)
}


## Heat Map (modified from cytof_heatmap to show all percentage heat maps in groups with same scale

cytof_heatmap_LC <- function (data, baseName = "Cluster", scaleMethod = "none", dendrogram = c("both", 
    "row", "column", "none"), colPalette = c("bluered", "greenred", 
    "spectral1", "spectral2"), cex_row_label = NULL, cex_col_label = NULL, 
    key.par = list(mgp = c(1.5, 0.5, 0), mar = c(3, 2.5, 3.5, 
        1)), keysize = 1.4, margins = c(6, 6), max=100, min=0) 
{
    data <- as.matrix(data)
    dendrogram <- match.arg(dendrogram)
    colPalette <- match.arg(colPalette)
    if (is.null(cex_row_label)) {
        cex_row_label <- (11 - ceiling(nrow(data)/10))/10
    }
    if (is.null(cex_col_label)) {
        cex_col_label <- (11 - ceiling(ncol(data)/10))/10
    }
    if (dendrogram == "row") {
        dendrogramRowv <- TRUE
        dendrogramColv <- FALSE
    }
    else if (dendrogram == "column") {
        dendrogramRowv <- FALSE
        dendrogramColv <- TRUE
    }
    else if (dendrogram == "none") {
        dendrogramRowv <- FALSE
        dendrogramColv <- FALSE
    }
    else {
        dendrogramRowv <- TRUE
        dendrogramColv <- TRUE
    }
  
    
    par(cex.main=0.8)
    heatmap.2(x = data, Rowv = dendrogramRowv, Colv = dendrogramColv, 
        dendrogram = dendrogram, col = colPalette, trace = "none", 
        symbreaks = FALSE, scale = scaleMethod, cexRow = cex_row_label, 
        cexCol = cex_col_label, srtCol = 30, symkey = FALSE, 
        key.par = key.par, margins = margins, keysize = keysize, 
        main = paste(baseName, "Heat Map"), breaks=seq(min, max, length.out=51))
}


#- A function to normalize expression values to a 0-1 range
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

#- Leave quietly in --array mode
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}



#-----------------
#- Get input data
#-----------------

## @knitr fcs

files <- list.files(inputfiles,pattern='.fcs$', full=TRUE)
files_short <- list.files(inputfiles,pattern='.fcs$', full=F)
usermarkers <- as.character(read.table(markersFile, header = FALSE)[,1])


## @knitr fcs1

fcs1<-read.FCS(files[1])

#— Ideally the user should upload the marker names as provided in the “Description” column of Flowjo. However, if the user uploads the shorter version (i.e., CD38 instead of 141Pr_CD38), the below will do the trick. This is irrelevant for flow data, as the description field given by Flowjo for Flow is always in the short form
#- Also, the colnames in the expression data are on the “Name<Desc>” form, so I am making another hash to be able to substitute these colnames for the markers given by the user (i.e., description or the short version of the description)
 
allMarkerNames<-pData(parameters(fcs1))$name
allMarkerDesc<-pData(parameters(fcs1))$desc
allMarkerNameAndDesc<-paste0(allMarkerNames, "<", allMarkerDesc, ">")
 
UserName2Desc <-hash()
Desc2UserName <-hash()
for(i in 1:length(allMarkerDesc)){
	if(is.na(allMarkerDesc[i])){
 		UserName2Desc[[ allMarkerNames[i] ]] <- allMarkerNames[i]
 		Desc2UserName[[ allMarkerNames[i] ]] <- allMarkerNames[i]
	}else{
		UserName2Desc[[ allMarkerDesc[i] ]] <- allMarkerDesc[i]
 		Desc2UserName[[ allMarkerDesc[i] ]] <- allMarkerDesc[i]
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
 
Desc2NameDesc <-hash()
NameDesc2Desc <- hash()
for(i in 1:length(allMarkerDesc)){
	if(is.na(allMarkerDesc[i])){
		Desc2NameDesc[[ allMarkerNames[i] ]] <- allMarkerNameAndDesc[i]
		NameDesc2Desc[[ allMarkerNameAndDesc[i] ]] <- allMarkerNames[i]
	}else{
		Desc2NameDesc[[ allMarkerDesc[i] ]] <- allMarkerNameAndDesc[i]
		NameDesc2Desc[[ allMarkerNameAndDesc[i] ]] <- allMarkerDesc[i]
	}
}
  
NameDesc2UserName<-hash()
for(i in 1:length(NameDesc2Desc)){
	Desc=values(NameDesc2Desc, keys= allMarkerNameAndDesc[i])
	NameDesc2UserName[[ allMarkerNameAndDesc[i] ]] <- values(Desc2UserName, keys= Desc)
}  

markersDesc <-vector()
for(i in 1:length(usermarkers)){
	markersDesc[i]<-values(UserName2Desc,  keys=usermarkers[i])
}
markersNameDesc <-vector()
for(i in 1:length(markersDesc)){
	markersNameDesc[i]<-values(Desc2NameDesc, keys=markersDesc[i])
}
markersUserName <-vector()
for(i in 1:length(markersNameDesc)){
	markersUserName[i] <-values(NameDesc2UserName, keys=markersNameDesc[i])
}


#------------------------------------------------------------------
#- Parse config file
#------------------------------------------------------------------

## @knitr parseConfig

projectName = "cytofpipe"

dimReductionMethod="tsne"
clusterMethods<-vector()
visualizationMethods<-vector()
visualizationMethods<-c(visualizationMethods,"tsne")
groups<-hash()

config<-read.ini(configFile)

autogating=config$clustering$GATING
if(transformMethod == '-'){transformMethod = config$clustering$TRANSFORM}
if(mergeMethod == '-'){mergeMethod = config$clustering$MERGE}
if(fixedNum == '-'){
	if(is.null(config$clustering$DOWNSAMPLE)){
		fixedNum = 100
	}else{
		fixedNum = config$clustering$DOWNSAMPLE
	}
}
if(displayAll == '-'){
	if(is.null(config$clustering$DISPLAY_ALL)){
		displayAll = "no"
	}else{
		displayAll = config$clustering$DISPLAY_ALL
	}
}
if(basename(groupfile) != '-'){

	conditions<-read.table(groupfile, header=F, sep="\t")
	colnames(conditions)<-c("Sample","Group")
	
	for(i in 1:length(conditions$Sample)){
		groups[[ sub(".fcs$", "", conditions$Sample[i], ignore.case=T) ]] <- as.character(conditions$Group[i])
	}
}
if(randomSampleSeed == '-'){
	if(is.null(config$clustering$RANDOM_SAMPLE_SEED)){
		randomSampleSeed = "no"
	}else{
		randomSampleSeed = config$clustering$RANDOM_SAMPLE_SEED
	}
}
if(randomTsneSeed == '-'){
	if(is.null(config$clustering$RANDOM_TSNE_SEED)){
		randomTsneSeed = "no"
	}else{
		randomTsneSeed = config$clustering$RANDOM_TSNE_SEED
	}
}
if(randomFlowSeed == '-'){
	if(is.null(config$clustering$RANDOM_FLOW_SEED)){
		randomFlowSeed = "no"
	}else{
		randomFlowSeed = config$clustering$RANDOM_FLOW_SEED
	}
}
flowsom_num = 15
perplexity = 30
theta = 0.5
max_iter = 1000


if(length(config$clustering$PERPLEXITY)==1){perplexity=config$clustering$PERPLEXITY}
if(length(config$clustering$THETA)==1){theta=config$clustering$THETA}
if(length(config$clustering$MAX_ITER)==1){max_iter=config$clustering$MAX_ITER}

if(length(config$clustering$PHENOGRAPH)==1){tolower(config$clustering$PHENOGRAPH);if(config$clustering$PHENOGRAPH == "yes"){clusterMethods<-c(clusterMethods,"Rphenograph")}}
if(length(config$clustering$CLUSTERX)==1){tolower(config$clustering$CLUSTERX);if(config$clustering$CLUSTERX == "yes"){clusterMethods<-c(clusterMethods,"ClusterX")}}
if(length(config$clustering$DENSVM)==1){tolower(config$clustering$DENSVM);if(config$clustering$DENSVM == "yes"){clusterMethods<-c(clusterMethods,"DensVM")}}
if(length(config$clustering$FLOWSOM)==1){tolower(config$clustering$FLOWSOM);if(config$clustering$FLOWSOM == "yes"){clusterMethods<-c(clusterMethods,"FlowSOM");flowsom_num=config$clustering$FLOWSOM_K}}
if(length(clusterMethods) == 0){clusterMethods<-c(clusterMethods,"NULL")}

if(length(config$clustering$PCA)==1){tolower(config$clustering$PCA);if(config$clustering$PCA == "yes"){visualizationMethods<-c(visualizationMethods,"pca")}}
if(length(config$clustering$ISOMAP)==1){tolower(config$clustering$ISOMAP);if(config$clustering$ISOMAP == "yes"){visualizationMethods<-c(visualizationMethods,"isomap")}}


#------------------------------------------------------------------
#- Do automatic gating
#------------------------------------------------------------------

## @knitr gating

if(autogating == 'yes'){

	gt<-gatingTemplate(template)

	#------------------------------------------------------------------------------------------------
	#- gates are based on arcSinh transformed data, so raw files need to be transformed before gating
	#------------------------------------------------------------------------------------------------
	
	arcsinh <- arcsinhTransform("arcsinh transformation")	
	dataTransform <- transform(read.flowSet(files), 
			"arcsinh_Ce142Di"= arcsinh(Ce142Di),
			"arcsinh_Ce140Di"= arcsinh(Ce140Di),
			"arcsinh_Ir191Di"= arcsinh(Ir191Di),
			"arcsinh_Ir193Di"= arcsinh(Ir193Di),
			"arcsinh_Y89Di"= arcsinh(Y89Di),
			"arcsinh_Pt195Di"= arcsinh(Pt195Di)
	)

	gs <- GatingSet(dataTransform)
	gating(x = gt, y = gs)
	fs_live <- getData(gs,"Live")

	pdf(paste0(outputdir,"/gating_scheme.pdf"))
	plot(gs)
	dev.off()

	write.flowSet(fs_live, paste0(outputdir, "/gating_fs_live"))

	rm(files)
	rm(files_short)

	files<-list.files(paste0(outputdir, "/gating_fs_live"), patter=".fcs", full=T)
	files_short<-list.files(paste0(outputdir, "/gating_fs_live"), patter=".fcs", full=F)

	for(i in 1:length(files_short)){
		pdf(paste0(outputdir,"/gating_",files_short[i],".pdf"))
		plotGate(gs[[i]])
		dev.off()
	}

}


#------------------------------------------------------------------
#- Run cytofkit wraper
#------------------------------------------------------------------

## @knitr cytofkit


#- Set seeds
sampleSeed=123
if(randomSampleSeed == 'yes'){
	sampleSeed=sample(1:1000000,1)
}
tsneSeed=42
if(randomTsneSeed == 'yes'){
	tsneSeed=sample(1:1000000,1)
}
flowSeed=100
if(randomFlowSeed == 'yes'){
        flowSeed=sample(1:1000000,1)
}


#- cytof_exprsMerge calls cytof_exprsExtract, which excludes Time and Event channels from the expression matrix, and excludes FSC/SSC from transformation
#- By default has sampleSeed = 123
exprs_data <- cytof_exprsMerge(fcsFiles = files, comp = FALSE, verbose = FALSE, 
                                   transformMethod = transformMethod, sampleSeed = sampleSeed,
                                   mergeMethod = mergeMethod, fixedNum = as.numeric(fixedNum))

#- change the colnames here so that the plots show the markers as uploaded by the user
for(i in 1:length(colnames(exprs_data))){
	colnames(exprs_data)[i]<-values(NameDesc2UserName, keys=colnames(exprs_data)[i])

 }


#- If array mode, just get the median_data and cluster_percentage files and leave quietly
#- This avaiods going through the tsne step which takes a long time
#- This is only if clustering is not DensVM or ClusterX, as these need the tSNE data

if(array == 'yes'){

	if(!dir.exists(outputdir)){
    		dir.create(outputdir)
	}
	curwd <- getwd()
	setwd(outputdir)

	exprs <- exprs_data[, markersUserName]
	ifMultiFCS <- length(unique(sub("_[0-9]*$", "", row.names(exprs)))) > 1
	
	if("DensVM" %in% clusterMethods || "ClusterX" %in% clusterMethods){

		## dimension reduced data, a list
		#- By default has tsneSeed = 42
		alldimReductionMethods <- unique(c(visualizationMethods, dimReductionMethod))
		allDimReducedList <- lapply(alldimReductionMethods,
		                             cytof_dimReduction, data = exprs_data,
		                             markers = markersUserName,
		                             tsneSeed = tsneSeed,
		                             perplexity = as.numeric(perplexity),
		                             theta = as.numeric(theta),
		                             max_iter = as.numeric(max_iter))
		names(allDimReducedList) <- alldimReductionMethods
				
		## cluster results, a list
		#- by default has flowSeed = NULL, I was using flowSeed=100 in v1.0 as that's what they used when they changed the code to make FlowSOM reproducible
		cluster_res <- lapply(clusterMethods, cytof_cluster,
		                          ydata = allDimReducedList[[dimReductionMethod]],
		                          xdata = exprs_data[, markersUserName],
		                          FlowSOM_k = as.numeric(flowsom_num),
		                          flowSeed = flowSeed)
	}else{
		
		cluster_res <- lapply(clusterMethods, cytof_cluster, 
	                          ydata = NULL, 
	                          xdata = exprs_data[, markersUserName],
	                          FlowSOM_k = as.numeric(flowsom_num),
	                          flowSeed = flowSeed)
		
	}
	names(cluster_res) <- clusterMethods
	clusterData <- cluster_res
	
	## save clusterData
	if(!is.null(clusterData) && length(clusterData) > 0){
		for(j in 1:length(clusterData)){
			methodj <- names(clusterData)[j]
			dataj <- clusterData[[j]]
			if(!is.null(dataj)){
	
				exprs_cluster_sample <- data.frame(exprs, cluster = dataj, check.names = FALSE)
	                   
				## cluster median
				cluster_median <- cytof_clusterStat(data= exprs_cluster_sample, cluster = "cluster", statMethod = "median")
				write.csv(cluster_median, paste(projectName, methodj, "cluster_median_data.csv", sep = "_"))
	                     
				## cluster percentage
				if (ifMultiFCS) {
					cluster_percentage <- cytof_clusterStat(data= exprs_cluster_sample, cluster = "cluster", statMethod = "percentage")
					write.csv(cluster_percentage, paste(projectName, methodj, "cluster_cell_percentage.csv", sep = "_"))
				}
			}
		}
	}
	stop_quietly()
}


## dimension reduced data, a list
#- By default has tsneSeed = 42
alldimReductionMethods <- unique(c(visualizationMethods, dimReductionMethod))
allDimReducedList <- lapply(alldimReductionMethods, 
                             cytof_dimReduction, data = exprs_data, 
			     markers = markersUserName,
			     tsneSeed = tsneSeed,
			     perplexity = as.numeric(perplexity),
			     theta = as.numeric(theta),
			     max_iter = as.numeric(max_iter))
names(allDimReducedList) <- alldimReductionMethods


## cluster results, a list
#- by default has flowSeed = NULL, I was using flowSeed=100 in v1.0 as that's what they used when they changed the code to make FlowSOM reproducible
cluster_res <- lapply(clusterMethods, cytof_cluster, 
                          ydata = allDimReducedList[[dimReductionMethod]], 
                          xdata = exprs_data[, markersUserName],
                          FlowSOM_k = as.numeric(flowsom_num),
                          flowSeed = flowSeed)
names(cluster_res) <- clusterMethods



## select the markers to display in the plots
displayMarkers <- vector()
if(displayAll == 'yes'){
	exclude <- grep("FSC|SSC|viability", colnames(exprs_data), ignore.case = TRUE)
	displayMarkers <- setdiff(colnames(exprs_data), colnames(exprs_data)[exclude])
}else{
	displayMarkers = markersUserName
}


## wrap the results
message("Stashing sample names...")
names <- sub("^.+/", "", unique(sub(".fcs$", "", files)))
samples <- as.list(NULL)
for(i in seq_along(names)){
	samples[[i]] <- names[i]
}
analysis_results <- list(expressionData = exprs_data[,displayMarkers],
                             dimReductionMethod = dimReductionMethod,
                             visualizationMethods = alldimReductionMethods,
                             dimReducedRes = allDimReducedList,
                             clusterRes = cluster_res, 
                             projectName = projectName,
                             rawFCSdir = inputfiles,
                             resultDir = outputdir,
			     dimRedMarkers = markersUserName,
			     sampleNames = samples)
        
## save the results
cytof_writeResults(analysis_results = analysis_results,
                       saveToRData = TRUE,
                       saveToFCS = TRUE,
                       saveToFiles = TRUE)


#------------------------------------------------------------------
#- Get scaled and norm01 heatmaps for median and percentage
#-	and level Plots
#------------------------------------------------------------------

## @knitr scaledHeatmaps

exprs <- as.data.frame(analysis_results$expressionData)
clusterData <- analysis_results$clusterRes
dimRed<-as.data.frame(analysis_results$dimReducedRes)

numFCS <- length(unique(sub("_[0-9]*$", "", row.names(exprs))))
ifMultiFCS <- length(unique(sub("_[0-9]*$", "", row.names(exprs)))) > 1

visualizationData <- analysis_results$dimReducedRes[analysis_results$visualizationMethods]

## Level plots
data_all<-cbind(exprs, dimRed)

visualizationData <- analysis_results$dimReducedRes[analysis_results$visualizationMethods]
for(i in 1:length(visualizationData)){
	if(!is.null(visualizationData[[i]])){
		methodi <- names(visualizationData)[i]
		datai <- as.data.frame(visualizationData[[i]])
		
		## Level plots
		pdf(paste0(outputdir,"/",projectName, "_", methodi, "_level_plot.pdf"))
		gp<-cytof_wrap_colorPlot(data=data_all,xlab=paste0(methodi,".", methodi,"_1"), ylab=paste0(methodi,".", methodi, "_2"), markers=displayMarkers, colorPalette = c("spectral1"), pointSize=0.1)
		print(gp)
		dev.off()

		## Level plots per file
		#- get min and max expression values for marker level plot colour scale
		min=min(apply(exprs,2,min))
		max=max(apply(exprs,2,max))

		## if multiple files, do level plot per sample, and if groups, also per group
		if (ifMultiFCS) {
			# combine datai and dataj
			xlab <- colnames(datai)[1]
			ylab <- colnames(datai)[2]
			dataik<- datai
			dataik$sample <- sub("_[0-9]*$", "", row.names(dataik))
			sample <- "sample"
			
			X<-split(dataik, dataik$sample)				
			for (d in 1:length(X)){
				samplename=X[[d]]$sample[1]
				data_all_sample <- subset(data_all, rownames(data_all) %in% rownames(X[[d]]))
	
				#- so that all the plots have the same x and y scales
				xlab=paste0(methodi,".", methodi,"_1")
				ylab=paste0(methodi,".", methodi,"_2")
							
				range.x<-max(data_all[xlab])-min(data_all[xlab])
				range.y<-max(data_all[ylab])-min(data_all[ylab])
				xlim=c(min(data_all[xlab]), max(data_all[xlab]))
				ylim=c(min(data_all[ylab]), max(data_all[ylab]))
	
				pdf(paste0(outputdir,"/",projectName, "_", methodi,  "_", samplename,  "_sample_level_plot.pdf"))
				gp<-cytof_wrap_colorPlot(data=data_all_sample, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, markers=displayMarkers, colorPalette = c("spectral1"), pointSize=0.1, min=min, max=max)
				print(gp)
				dev.off()		
			}

			#- add a colum with the groups 
			if(basename(groupfile) != '-'){
				for(l in 1:length(dataik$sample)){
					dataik$group[l] <- values(groups, keys=dataik$sample[l])
				}
	
				## Level plots per group
				X2<-split(dataik, dataik$group)				
				for (d2 in 1:length(X2)){
					groupname=X2[[d2]]$group[1]
					data_all_group <- subset(data_all, rownames(data_all) %in% rownames(X2[[d2]]))
		
					#- so that all the plots have the same x and y scales
					xlab=paste0(methodi,".", methodi,"_1")
					ylab=paste0(methodi,".", methodi,"_2")
								
					range.x<-max(data_all[xlab])-min(data_all[xlab])
					range.y<-max(data_all[ylab])-min(data_all[ylab])
					xlim=c(min(data_all[xlab]), max(data_all[xlab]))
					ylim=c(min(data_all[ylab]), max(data_all[ylab]))
		
					pdf(paste0(outputdir,"/",projectName, "_", methodi,  "_", groupname,  "_group_level_plot.pdf"))
					gp<-cytof_wrap_colorPlot(data=data_all_group, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, markers=displayMarkers, colorPalette = c("spectral1"), pointSize=0.1, min=min, max=max)
					print(gp)
					dev.off()
				}
			}
		}
		## if multiple files, redo he cluster grid plot to correct label size
		if (ifMultiFCS) {
 			if(!is.null(clusterData) && length(clusterData) > 0){
				for(j in 1:length(clusterData)){
					if(!is.null(clusterData[[j]])){
						methodj <- names(clusterData)[j]
						dataj <- clusterData[[j]]
                        	    
						# combine datai and dataj
						xlab <- colnames(datai)[1]
						ylab <- colnames(datai)[2]
						dataij <- datai
						dataij$sample <- sub("_[0-9]*$", "", row.names(dataij))
						dataij$cluster <- factor(dataj)
						cluster <- "cluster"
						sample <- "sample"
                        	    
						## cluster grid plot if multiple files
						figName <- paste(projectName, methodi, methodj, sep=" ")
						labelsizesscaled=floor(10/numFCS)-1
						labelsizesscaled <- ifelse(labelsizesscaled > 2, labelsizesscaled , 2)
						
						pdf(paste0(outputdir,"/",projectName, "_", methodi, "_", methodj, "_cluster_grid_scatter_plot.pdf"))
						cluster_grid_plot <- cytof_clusterPlot(dataij, xlab, ylab, cluster, sample, figName, 2, point_size =0.5, labelSize= labelsizesscaled)
						print(cluster_grid_plot)
						dev.off()

					}
				}
			}
		}  
	}
}


if(!is.null(clusterData) && length(clusterData) > 0){

	## Heatmaps
	for(j in 1:length(clusterData)){
		methodj <- names(clusterData)[j]
		dataj <- clusterData[[j]]
		if(!is.null(dataj)){
                    
			exprs_cluster_sample <- data.frame(exprs[, displayMarkers], cluster = dataj, check.names = FALSE)
		
			## cluster median 
			cluster_median <- cytof_clusterStat(data=exprs_cluster_sample, cluster = "cluster", statMethod = "median")

			## Heatmap scaled
			pdf(paste0(outputdir,"/",projectName, "_",methodj, "_cluster_median_heatmap_scaled.pdf"))
			cytof_heatmap(cluster_median, scaleMethod="column", paste(projectName, methodj, "\ncluster median (scaled)", sep = " "))
			dev.off()

			## Heatmap norm01
			cluster_median_norm01<-as.data.frame( apply(cluster_median, 2, range01))
			pdf(paste0(outputdir,"/",projectName, "_",methodj, "_cluster_median_heatmap_norm01.pdf"))
			cytof_heatmap(cluster_median_norm01, paste(projectName, methodj, "\ncluster median (norm01)", sep = " "))
			dev.off()

			if(basename(groupfile) != '-'){
				## cluster percentage
				if (ifMultiFCS) {
					cluster_percentage <- cytof_clusterStat(data=exprs_cluster_sample, cluster = "cluster", statMethod = "percentage")
	
					exprs_cluster_sample$sample <- sub("_[0-9]*$", "", row.names(exprs_cluster_sample))
					sample <- "sample"
	                        	#- add a colum with the groups 
					for(l in 1:length(exprs_cluster_sample$sample)){
						exprs_cluster_sample$group[l] <- values(groups, keys= exprs_cluster_sample$sample[l])
					}		
				
					## Percentage plots per group
					X3<-split(exprs_cluster_sample, exprs_cluster_sample$group)				
					for (d3 in 1:length(X3)){
						groupname=X3[[d3]]$group[1]
						exprs_group <- subset(exprs_cluster_sample, rownames(exprs_cluster_sample) %in% rownames(X3[[d3]]))
						exprs_group$sample <- NULL
						exprs_group$group <- NULL
						cluster_percentage_group <- cytof_clusterStat(data=exprs_group, cluster = "cluster", statMethod = "percentage")

						write.table(cluster_percentage_group, file=paste0(outputdir,"/",projectName, "_",methodj, "_", groupname, "_group_cluster_percentage_data.csv"), col.names=NA, sep=",")
						if (ncol(cluster_percentage_group) > 1){
							pdf(paste0(outputdir,"/",projectName, "_",methodj, "_", groupname, "_group_cluster_percentage_heatmap.pdf"))
							cytof_heatmap_LC(cluster_percentage_group,paste(projectName, methodj, groupname,"cluster\ncell percentage", sep = " "), max=max(cluster_percentage), 
min=min(cluster_percentage))	
							dev.off()
						}
					}		
                      	    	}	

			}

		}
	}
}

paste0("files: ", files)
paste0("transformMethod: ", transformMethod)
paste0("markersUserName: ", markersUserName)
paste0("groupfile: ", groupfile)
if(basename(groupfile) != '-'){
	groups
}
paste0("displayMarkers: ", displayMarkers)
paste0("mergeMethod: ", mergeMethod)
paste0("fixedNum: ", as.numeric(fixedNum))
paste0("randomSampleSeed: ", randomSampleSeed," (sampleSeed = ", sampleSeed,")")
paste0("randomTsneSeed: ", randomTsneSeed," (tsneSeed = ", tsneSeed,")")
paste0("randomFlowSeed: ", randomFlowSeed," (flowseed = ", flowSeed,")")
paste0("perplexity: ", as.numeric(perplexity))
paste0("theta: ", as.numeric(theta))
paste0("max_iter: ", as.numeric(max_iter))


sessionInfo()

