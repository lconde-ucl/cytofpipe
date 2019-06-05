required.packages <- c("cytofkit2","ini","hash","flowCore","reshape2")
not.installed <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]

if(length(not.installed)>0){
	paste0("ERROR: package not installed: ", not.installed)
}

