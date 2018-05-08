required.packages <- c("cytofkit","ini","hash","flowCore","openCyto","mvtnorm","reshape2","VGAM","colourpicker","gplots", "scaffold","tools","igraph","reshape","ggrepel","citrus")
not.installed <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]

if(length(not.installed)>0){
	paste0("ERROR: package not installed: ", not.installed)
}else{
	paste0("All required R packages installed")
}

