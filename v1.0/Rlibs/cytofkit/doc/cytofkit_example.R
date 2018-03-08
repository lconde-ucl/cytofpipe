## ---- eval=FALSE-----------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("cytofkit")

## ---- message=FALSE--------------------------------------------------------
library("cytofkit") 

## ---- eval=FALSE-----------------------------------------------------------
#  ?"cytofkit-package"

## ---- eval=FALSE-----------------------------------------------------------
#  cytofkit_GUI()

## ---- eval=FALSE-----------------------------------------------------------
#  set.seed(100)
#  dir <- system.file('extdata',package='cytofkit')
#  file <- list.files(dir ,pattern='.fcs$', full=TRUE)
#  parameters <- list.files(dir, pattern='.txt$', full=TRUE)
#  res <- cytofkit(fcsFiles = file,
#                  markers = parameters,
#                  projectName = 'cytofkit_test',
#                  transformMethod = "cytofAsinh",
#                  mergeMethod = "ceil",
#                  fixedNum = 500,                                    ## set at 500 for faster run
#                  dimReductionMethod = "tsne",
#                  clusterMethods = c("Rphenograph", "ClusterX"),    ## accept multiple methods
#                  visualizationMethods = c("tsne", "pca"),          ## accept multiple methods
#                  progressionMethod = "isomap",
#                  clusterSampleSize = 500,
#                  resultDir = getwd(),
#                  saveResults = TRUE,
#                  saveObject = TRUE)

## --------------------------------------------------------------------------
## Loading the FCS data:  
dir <- system.file('extdata',package='cytofkit')
file <- list.files(dir ,pattern='.fcs$', full=TRUE)
paraFile <- list.files(dir, pattern='.txt$', full=TRUE)
parameters <- as.character(read.table(paraFile, header = TRUE)[,1])

## File name
file

## parameters
parameters

## --------------------------------------------------------------------------
## Extract the expression matrix with transformation
data_transformed <- cytof_exprsExtract(fcsFile = file, 
                                       comp = FALSE, 
                                       transformMethod = "cytofAsinh")
## If analysing flow cytometry data, you can set comp to TRUE or 
## provide a transformation matrix to apply compensation

## If you have multiple FCS files, expression can be extracted and combined
combined_data_transformed <- cytof_exprsMerge(fcsFiles = file, comp=FALSE,
                                              transformMethod = "cytofAsinh",
                                              mergeMethod = "all")
## change mergeMethod to apply different combination strategy

## Take a look at the extracted expression matrix
head(data_transformed[ ,1:3])

## ---- message=FALSE--------------------------------------------------------
## use clustering algorithm to detect cell subsets
## to speed up our test here, we only use 100 cells
data_transformed_1k <- data_transformed[1:100, ]

## run PhenoGraph
cluster_PhenoGraph <- cytof_cluster(xdata = data_transformed_1k, method = "Rphenograph")

## run ClusterX
data_transformed_1k_tsne <- cytof_dimReduction(data=data_transformed_1k, method = "tsne")
cluster_ClusterX <- cytof_cluster(ydata = data_transformed_1k_tsne,  method="ClusterX")

## ---- eval=FALSE-----------------------------------------------------------
#  ## run DensVM (takes long time, we skip here)
#  cluster_DensVM <- cytof_cluster(xdata = data_transformed_1k,
#                                  ydata = data_transformed_1k_tsne, method = "DensVM")

## ---- message=FALSE--------------------------------------------------------
## run FlowSOM with cluster number 15
cluster_FlowSOM <- cytof_cluster(xdata = data_transformed_1k, method = "FlowSOM", FlowSOM_k = 12)

## combine data
data_1k_all <- cbind(data_transformed_1k, data_transformed_1k_tsne, 
                     PhenoGraph = cluster_PhenoGraph, ClusterX=cluster_ClusterX, 
                     FlowSOM=cluster_FlowSOM)
data_1k_all <- as.data.frame(data_1k_all)

## ---- message=FALSE--------------------------------------------------------
## PhenoGraph plot on tsne
cytof_clusterPlot(data=data_1k_all, xlab="tsne_1", ylab="tsne_2", 
                  cluster="PhenoGraph", sampleLabel = FALSE)

## PhenoGraph cluster heatmap
PhenoGraph_cluster_median <- aggregate(. ~ PhenoGraph, data = data_1k_all, median)
cytof_heatmap(PhenoGraph_cluster_median[, 2:37], baseName = "PhenoGraph Cluster Median")

## --------------------------------------------------------------------------
## ClusterX plot on tsne
cytof_clusterPlot(data=data_1k_all, xlab="tsne_1", ylab="tsne_2", cluster="ClusterX", sampleLabel = FALSE)

## ClusterX cluster heatmap
ClusterX_cluster_median <- aggregate(. ~ ClusterX, data = data_1k_all, median)
cytof_heatmap(ClusterX_cluster_median[, 2:37], baseName = "ClusterX Cluster Median")

## --------------------------------------------------------------------------
## FlowSOM plot on tsne
cytof_clusterPlot(data=data_1k_all, xlab="tsne_1", ylab="tsne_2", 
                  cluster="FlowSOM", sampleLabel = FALSE)

## FlowSOM cluster heatmap
FlowSOM_cluster_median <- aggregate(. ~ FlowSOM, data = data_1k_all, median)
cytof_heatmap(FlowSOM_cluster_median[, 2:37], baseName = "FlowSOM Cluster Median")

## ---- message=FALSE--------------------------------------------------------
## Inference of PhenoGraph cluster relatedness
PhenoGraph_progression <- cytof_progression(data = data_transformed_1k, 
                                            cluster = cluster_PhenoGraph, 
                                            method="isomap", clusterSampleSize = 50, 
                                            sampleSeed = 5)
p_d <- data.frame(PhenoGraph_progression$sampleData, 
                  PhenoGraph_progression$progressionData, 
                  cluster = PhenoGraph_progression$sampleCluster, 
                  check.names = FALSE)

## cluster relatedness plot
cytof_clusterPlot(data=p_d, xlab="isomap_1", ylab="isomap_2", 
                  cluster="cluster", sampleLabel = FALSE)

## marker expression profile
markers <- c("(Sm150)Di<GranzymeB>", "(Yb173)Di<Perforin>")

cytof_colorPlot(data=p_d, xlab="isomap_1", ylab="isomap_2", zlab = markers[1], limits = range(p_d[,1:52]))
cytof_colorPlot(data=p_d, xlab="isomap_1", ylab="isomap_2", zlab = markers[2], limits = range(p_d[,1:52]))

cytof_progressionPlot(data=p_d, markers=markers, orderCol="isomap_1", clusterCol = "cluster")

## ---- message=FALSE--------------------------------------------------------
## Inference of ClusterX cluster relatedness
ClusterX_progression <- cytof_progression(data = data_transformed_1k, 
                                          cluster = cluster_ClusterX, 
                                          method="isomap", 
                                          clusterSampleSize = 30, 
                                          sampleSeed = 3)
c_d <- data.frame(ClusterX_progression$sampleData, 
                  ClusterX_progression$progressionData,
                  cluster=ClusterX_progression$sampleCluster, 
                  check.names = FALSE)

## cluster relatedness plot
cytof_clusterPlot(data=c_d, xlab="isomap_1", ylab="isomap_2", 
                  cluster="cluster", sampleLabel = FALSE)

## marker expression profile
markers <- c("(Sm150)Di<GranzymeB>", "(Yb173)Di<Perforin>")

cytof_colorPlot(data=c_d, xlab="isomap_1", ylab="isomap_2", zlab = markers[1], limits = range(c_d[,1:52]))
cytof_colorPlot(data=c_d, xlab="isomap_1", ylab="isomap_2", zlab = markers[2], limits = range(c_d[,1:52]))

cytof_progressionPlot(data=c_d, markers, orderCol="isomap_1", clusterCol = "cluster")

## ---- eval=FALSE-----------------------------------------------------------
#  ## save analysis results to FCS file
#  cytof_addToFCS(data_1k_all, rawFCSdir=dir, analyzedFCSdir="analysed_FCS",
#                 transformed_cols = c("tsne_1", "tsne_2"),
#                 cluster_cols = c("PhenoGraph", "ClusterX", "FlowSOM"))

## --------------------------------------------------------------------------
## See documentation, this function uses the output of main cytofkit cuntion as its input
?cytof_clusterMtrx

## --------------------------------------------------------------------------
cytofkitNews()

## --------------------------------------------------------------------------
sessionInfo()

