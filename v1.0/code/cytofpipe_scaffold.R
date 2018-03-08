## @knitr libraries

library(scaffold)
library(flowCore)
library(tools)
library(igraph)
library(reshape)
library(ggrepel)
library(hash)
library(ini)

options(stringsAsFactors = F)
rm(list = ls())


#------------------------------------------------------------------
#- Parse arguments
#------------------------------------------------------------------

## @knitr arguments

jobid <- as.character(Sys.getenv("RAND_ID"))
input <- paste0(jobid, ".txt")

args<-read.ini(input)

working.dir=args$paramsscaffold$INPUTFILE
ref.file=paste0(args$paramsscaffold$REF, ".clustered.txt")
outputdir=args$paramsscaffold$OUTPUTFILE
markersFile=args$paramsscaffold$MARKERSFILE
asinh_cofactor=args$paramsscaffold$ASINH
mergeMethod = args$paramsscaffold$MERGE
fixedNum = args$paramsscaffold$DOWNSAMPLE

if(asinh_cofactor == '-'){
	asinh_cofactor=5;
}

#——————————————————————
#- DOWNSAMPLE FCS FILES
#——————————————————————

#— Downsample to 10000 events to speed up clustering (unless they add the --all parameter)

files.list <- list.files(path = working.dir, pattern = "*.fcs$")
if(mergeMethod == '-'){

	if(fixedNum == '-'){fixedNum = 10000}

	downsampled.dir<-paste0(outputdir,"/downsampled_",fixedNum)
	dir.create(file.path(outputdir, paste0("downsampled_",fixedNum)))
	
	downsampleFCS <- function(x) {
		x[sample(nrow(x), size = as.numeric(fixedNum), replace = ifelse(nrow(x) < as.numeric(fixedNum), TRUE, FALSE)), , drop=FALSE]
	}
	
	for(file in files.list) {
		outfile<-basename(file_path_sans_ext(file))
		fcs<-read.FCS(paste0(working.dir, "/", file))
	
		exp<-exprs(fcs)
		downsampled<-downsampleFCS(exp)
		exprs(fcs)<-downsampled
		write.FCS(fcs, filename=paste0(downsampled.dir, "/", outfile, ".fcs"))
	}
}else{
	downsampled.dir<-working.dir
}

#————————————
#- CLUSTERING
#————————————


usermarkers <- as.character(read.table(markersFile, header = FALSE)[,1])
fcs1<-read.FCS(paste0(working.dir, "/", files.list[1]))

allMarkerNames<-pData(parameters(fcs1))$name
allMarkerDesc<-pData(parameters(fcs1))$desc
 
UserName2Desc <-hash()
for(i in 1:length(allMarkerDesc)){
	UserName2Desc[[ allMarkerDesc[i] ]] <- allMarkerDesc[i]
}
if (sum(has.key( usermarkers, UserName2Desc )) == 0) {
 	clear(UserName2Desc)
 	for(i in 1:length(allMarkerDesc)){
 		id <- gsub( "^[^_]+_", "", allMarkerDesc[i])
 		UserName2Desc[[ id ]] <- allMarkerDesc[i]
 	}
}

markersDesc<-vector()
for(i in 1:length(usermarkers)){
  	markersDesc[i]<-values(UserName2Desc, keys=usermarkers[i])
}
 

## @knitr fixedparameters

num.cores = 1
num_clusters = 200
num_samples = 50

## @knitr parameters

col.names = markersDesc
asinh.cofactor = as.numeric(asinh_cofactor)

scaffold:::cluster_fcs_files_in_dir(downsampled.dir, num.cores, col.names, num_clusters, num_samples, asinh.cofactor)


#——————————————
#- RUN SCAFFOLD
#——————————————


ew_influence<- ceiling(length(col.names) / 3)

files.list <- list.files(path = downsampled.dir, pattern = "*.clustered.txt$")
files.list <- files.list[files.list != ref.file]
print(sprintf("Markers used for SCAFFoLD: %s", paste(col.names, collapse = ", ")))

files.list <- c(ref.file, files.list)
print(paste("Using as reference", files.list[1], sep = " "))

ref.dir <- paste(working.dir, "gated/", sep = "/")
gated_data <- scaffold:::load_attractors_from_gated_data(ref.dir, asinh.cofactor)
tab.attractors <- gated_data$tab.attractors
att.labels <- gated_data$cellType_key$population
G.attractors <- NULL

ret <- list(graphs = list(), clustered.data = list())

for(fi in files.list)
{
    f <- paste0(downsampled.dir,"/",fi)
    print (paste("Processing", f, sep = " "))
    tab <- read.table(f, header = T, sep = "\t", quote = "", check.names = F, comment.char = "", stringsAsFactors = F)
    col.names.inter_cluster <- col.names
    tab <- tab[!apply(tab[, col.names], 1, function(x) {all(x == 0)}),]    
    names(tab) <- gsub("cellType", "groups", names(tab))
    names(tab) <- gsub("^X", "", names(tab))
    print(sprintf("Running with Edge weight: %f", ew_influence))
    res <- scaffold:::process_data(tab, G.attractors, tab.attractors,
         col.names = col.names, att.labels = att.labels, already.clustered = T, ew_influence = ew_influence, 
         col.names.inter_cluster = col.names.inter_cluster, inter.cluster.connections = T, overlap_method = NULL)
    G.complete <- scaffold:::get_highest_scoring_edges(res$G.complete)
    clustered.data <- tab 
    ret$graphs[basename(f)] <- list(G.complete)
    ret$clustered.data[basename(f)] <- list(clustered.data)
        
    G.attractors <- res$G.attractors
}

ret <- c(ret, list(dataset.statistics = scaffold:::get_dataset_statistics(ret)))
ret <- c(list(scaffold.col.names = col.names, landmarks.data = gated_data$downsampled.data), ret)
ref=gsub(pattern = "\\.fcs.*", "", ref.file)
scaffold:::my_save(ret, paste(outputdir, "cytofpipe.scaffold", sep = "/"))


#———————————----———
#- SCAFFOLD PLOTS
#——————————————---

f_name <- paste(outputdir, "cytofpipe.scaffold", sep = "/")

con <- file(f_name, "rb")
data <- unserialize(con)
close(con)

#- Get min and max x and y coordinates
x<-vector()
y<-vector()
for (i in 1:length(data$graphs)) {
	G <- data$graphs[[i]]
	
#	layout<-layout.auto(G)
#	layout<-layout.forceatlas2(G, iterations=1000, plotstep=500)
#	fixed <- rep(FALSE, vcount(G))
#	fixed[1:length(att.labels)] <- TRUE
#	layout<-scaffold:::layout.forceatlas2(G, fixed = fixed)
	
	layout<-cbind(V(G)$x,V(G)$y)  
	x<-c(x, layout[,1])
	y<-c(y, layout[,2])
}
range.x<-max(x)-min(x)
range.y<-max(y)-min(y)
xlim=c(min(x), max(x))
ylim=c(min(y), max(y))

for (i in 1:length(data$graphs)) {
	G <- data$graphs[[i]]
	
	name=names(data$graphs[i])
	name2=gsub(pattern = "\\.fcs.*", "", name)    

	#-get the node coordinates
#	plotcord <- data.frame(layout.auto(G))
	plotcord <- data.frame(cbind(V(G)$x,V(G)$y), row.names=V(G)$name)
	colnames(plotcord) = c("X1","X2")

	#-reverse data so that lanmark populatios are plotted last
	plotcord <- plotcord[rev(rownames(plotcord)),]
	
	#-get edges, which are pairs of node IDs
	edgelist <- get.edgelist(G)
	
	#-convert to a four column edge data frame with source and destination coordinates, and reverse
	edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
	colnames(edges) <- c("X1","Y1","X2","Y2")
	edges <- edges[rev(rownames(edges)),]

	#- labels, colours, node sizes.. and reverse
	labels<-c(V(G)$name[1:length(att.labels)], rep(NA, vcount(G)-length(att.labels)))
	labels<-rev(labels)
	colores<-c(rep(rgb(255/255,117/255,128/255, 0.9),length(att.labels)), rep(rgb(79/255, 147/255, 222/255, 0.3), vcount(G)-length(att.labels)))
	colores<-rev(colores)
#	V(G)$label.cex<-c(rep(1, length(att.labels)), rep(0.8, vcount(G)-length(att.labels)))
	popsize<-(V(G)$popsize/sum(V(G)$popsize, na.rm = T))*(range.x/8)
	popsize[is.na(V(G)$popsize)]<-quantile(popsize, .80,na.rm = TRUE)
	popsize <- rev(popsize)

	pdf(paste0(outputdir,"/scaffold_map_",name2,".pdf"))
	p<-ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey") + 
		geom_point(aes(X1, X2), size=popsize, colour=colores, data=plotcord) +
		geom_text_repel(aes(X1, X2),data=plotcord, label = labels) +
		xlim(xlim) +
		ylim(ylim) +
		theme(axis.line=element_blank(),
			axis.text.x=element_blank(),
			axis.text.y=element_blank(),
			axis.ticks=element_blank(),
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),
			legend.position="none",
			panel.background=element_blank(),
			panel.border=element_blank(),
			panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),
			plot.background=element_blank())

	print(p)
	dev.off()

}

paste0("FCS_files: ",files.list)
paste0("Ref: ",ref)
paste0("Clustering_markers: ",usermarkers)
if(mergeMethod == '-'){
	paste0("Events: ", fixedNum)
}else{
        paste0("Events: ", mergeMethod)
}
paste0("Asinh cofactor: ", asinh.cofactor)
paste0("Num. cores: ",num.cores)
paste0("Num. clusters: ",num_clusters)
paste0("Num. samples: ",num_samples)


sessionInfo()

