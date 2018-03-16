## ----requirements, echo=FALSE--------------------------------------------
if (!require(flowWorkspaceData)) {
  stop("Cannot build the vignettes without 'flowWorkspaceData'")
}

## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", message = FALSE, warning = FALSE)

## ----load-flowWorkspace, echo=F------------------------------------------
library(flowWorkspace)

## ----load-xml, eval=TRUE-------------------------------------------------
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
wsfile <- list.files(flowDataPath, pattern="manual.xml",full = TRUE)
wsfile

## ----openWorkspace, eval=F-----------------------------------------------
#  ws <- openWorkspace(wsfile)

## ----parseWorkspace, eval=F----------------------------------------------
#  gs <- parseWorkspace(ws, name= "T-cell", subset =1, isNcdf = TRUE)

## ----load_gs_manual, echo = FALSE----------------------------------------
gs <- load_gs(file.path(flowDataPath,"gs_manual"))

## ----plot-manual-GatingHierarchy-----------------------------------------
gh <- gs[[1]]
plot(gh)

## ----plot-manual-gates, fig.width = 9------------------------------------
plotGate(gh)

## ----gatingTemplate, eval = T--------------------------------------------
library(openCyto)
library(data.table)
gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
dtTemplate <- fread(gtFile, autostart = 1L)
dtTemplate

## ----gatingTemplate-nonDebris, eval = T----------------------------------
dtTemplate[1,]

## ----gatingTemplate-singlets, eval = T-----------------------------------
dtTemplate[2,]

## ----gatingTemplate-lympth, eval = T-------------------------------------
dtTemplate[3,]

## ----gatingTemplate-cd3, eval = T----------------------------------------
dtTemplate[4,]

## ----gatingTemplate-cd4cd8, eval = T-------------------------------------
dtTemplate[5,]

## ----gatingTemplate-expand, echo = F, results = F------------------------
expanded <- openCyto:::.preprocess_csv(dtTemplate)
rownames(expanded) <- NULL

## ----gatingTemplate-expand1, echo = F------------------------------------
expanded[5:6,]

## ----gatingTemplate-expand2, echo = F------------------------------------
expanded[7:10,]

## ----load-gt, eval = T---------------------------------------------------
gt_tcell <- gatingTemplate(gtFile, autostart = 1L)
gt_tcell

## ----plot-gt, eval = T---------------------------------------------------
plot(gt_tcell)

## ----load-fcs------------------------------------------------------------
fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
ncfs  <- read.ncdfFlowSet(fcsFiles)
fr <- ncfs[[1]]
gs <- GatingSet(ncfs)
gs

## ----compensate----------------------------------------------------------
compMat <- getCompensationMatrices(gh)
gs <- compensate(gs, compMat)

## ----compensate_plot, echo = F, fig.width = 4, fig.height = 4------------
library(ggcyto)
sub_chnl <- c("V545-A","V450-A")
fr <- fr[,sub_chnl]
fr_comp <- getData(gs[[1]])[,sub_chnl]
fs <- as(list(fr = fr, fr_comp = fr_comp), "flowSet")
#transform data to better visualize the compensation effect
fs <- transform(fs, estimateLogicle(fr,sub_chnl))
autoplot(fs, "V545", "V450")

## ----transformation, eval = T--------------------------------------------
chnls <- parameters(compMat)
trans <- estimateLogicle(gs[[1]], channels = chnls)
gs <- transform(gs, trans)

## ----transformation_plot, echo = F, fig.width = 5, fig.height = 5--------
fr_trans <- getData(gs[[1]])[,sub_chnl[1]]
fr_comp <- fr_comp[,sub_chnl[1]]
fs <- as(list(fr = fr_comp, fr_trans = fr_trans), "flowSet")
p1 <- as.ggplot(autoplot(fs[1], "V545"))
p2 <- as.ggplot(autoplot(fs[2], "V545"))
plot(gridExtra::arrangeGrob(p1,p2))

## ----gating, eval = TRUE-------------------------------------------------
gating(gt_tcell, gs)

## ----gating_par, eval = FALSE--------------------------------------------
#  gating(gt_tcell, gs, mc.cores=2, parallel_type = "multicore")

## ----plot_afterGating----------------------------------------------------
plot(gs[[1]])

## ----hideGate, results = "hide"------------------------------------------
dodesToHide <- c("cd8+", "cd4+"
				, "cd4-cd8-", "cd4+cd8+"
				, "cd4+cd8-/HLA+", "cd4+cd8-/CD38+"
				, "cd4-cd8+/HLA+", "cd4-cd8+/CD38+"
				, "CD45_neg/CCR7_gate", "cd4+cd8-/CD45_neg"
				, "cd4-cd8+/CCR7+", "cd4-cd8+/CD45RA+"
				)
lapply(dodesToHide, function(thisNode)setNode(gs, thisNode, FALSE))

## ----rename, results = "hide"--------------------------------------------
setNode(gs, "cd4+cd8-", "cd4")
setNode(gs, "cd4-cd8+", "cd8")

## ----plot_afterHiding----------------------------------------------------
plot(gs[[1]])

## ----plotGate_autoGate, fig.width = 9------------------------------------
plotGate(gs[[1]])

