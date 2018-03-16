## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----echo=FALSE, include=FALSE-------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE, warning = FALSE
)

## ----loadPackage, echo=FALSE,results='hide'------------------------------
library(flowCore)

## ----ReadFiles, echo=TRUE, results='markup'------------------------------
file.name <- system.file("extdata","0877408774.B08", package="flowCore")
x <- read.FCS(file.name, transformation=FALSE)
summary(x)

## ----SearchKeywords, echo=TRUE,results='markup'--------------------------
keyword(x,c("$P1E", "$P2E", "$P3E", "$P4E"))

## ----PrintSummary, echo=TRUE, results='markup'---------------------------
summary(read.FCS(file.name))

## ----PrintSummary2, echo=TRUE,results='markup'---------------------------
summary(read.FCS(file.name,transformation="scale")) 

## ----ReadFiles2, echo=TRUE,results='markup'------------------------------
read.FCS(file.name,alter.names=TRUE) 

## ----RedFiles3, echo=TRUE, results='markup'------------------------------
x <- read.FCS(file.name, column.pattern="-H") 
x 

## ----RedFiles4, echo=TRUE, results='markup'------------------------------
lines <- sample(100:500, 50)
y <- read.FCS(file.name, which.lines = lines) 
y 

## ----Plot2, echo=TRUE, results='hide'------------------------------------
library(ggcyto)
autoplot(x, "FL1-H", "FL2-H")

## ----plot3, echo=TRUE, results='hide'------------------------------------
autoplot(x, "FL1-H")

## ----Frames1, echo=TRUE, results='markup'--------------------------------
frames <- lapply(dir(system.file("extdata", "compdata", "data",
                           package="flowCore"), full.names=TRUE), 
                 read.FCS)
as(frames, "flowSet")

## ----Frames2, echo=TRUE,results='markup'---------------------------------
names(frames) <- sapply(frames, keyword, "SAMPLE ID")
fs <- as(frames, "flowSet")
fs

## ----metaData, echo=TRUE, results='markup'-------------------------------
phenoData(fs)$Filename <- fsApply(fs,keyword, "$FIL")
pData(phenoData(fs))

## ----ReadFlowSet, echo=TRUE,results='markup'-----------------------------
read.flowSet(path = system.file("extdata", "compdata", "data", 
             package="flowCore"))

## ----ReadFowSet2, echo=TRUE,results='markup'-----------------------------
fs <- read.flowSet(path=system.file("extdata", "compdata", "data",
                   package="flowCore"), name.keyword="SAMPLE ID",
                   phenoData=list(name="SAMPLE ID", Filename="$FIL"))
fs
pData(phenoData(fs))

## ----fsApply1, echo=TRUE, results='markup'-------------------------------
fsApply(fs, each_col, median)

## ----fsApply2, echo=TRUE,results='markup'--------------------------------
fsApply(fs,function(x) apply(x, 2, median), use.exprs=TRUE)

## ----Transfo1, echo=TRUE, results='hide'---------------------------------
autoplot(transform(fs[[1]]
                   , `FL1-H`=log(`FL1-H`)
                   , `FL2-H`=log(`FL2-H`)
                   )
         , "FL1-H","FL2-H")

## ----Transfo2, results='hide'--------------------------------------------
autoplot(transform(fs[[1]]
                   , log.FL1.H=log(`FL1-H`)
                   , log.FL2.H=log(`FL2-H`)
                   )
         , "log.FL1.H", "log.FL2.H")

## ----Transfo3, echo=TRUE-------------------------------------------------
aTrans <- truncateTransform("truncate at 1", a=1)
aTrans

## ----Transfo4, echo=TRUE,results='markup'--------------------------------
transform(fs,`FL1-H`=aTrans(`FL1-H`))

## ----Transfo4.1, echo=TRUE,results='markup'------------------------------
f1 <- function(fs,...){
  transform(fs, ...)[,'FL1-H']
}

f2 <- function(fs){
  aTrans <- truncateTransform("truncate at 1", a=1)
  f1(fs, `FL1-H` = aTrans(`FL1-H`))
}
res <- try(f2(fs), silent = TRUE)
res

## ----Transfo4.2, echo=TRUE,results='markup'------------------------------
myTrans <- transformList('FL1-H', aTrans)
transform(fs, myTrans)

## ----rectGate, echo=TRUE, results='markup'-------------------------------
rectGate <- rectangleGate(filterId="Fluorescence Region", 
                          "FL1-H"=c(0, 12), "FL2-H"=c(0, 12))

## ----echo=TRUE,results='markup'------------------------------------------
result = filter(fs[[1]],rectGate)
result

## ----Summary3, echo=TRUE, results='markup'-------------------------------
summary(result)
summary(result)$n
summary(result)$true
summary(result)$p

## ----SummarFilter, echo=TRUE, results='markup'---------------------------
summary(filter(fs[[1]], kmeansFilter("FSC-H"=c("Low", "Medium", "High"),
                                     filterId="myKMeans")))

## ----echo=TRUE,results='markup'------------------------------------------
filter(fs,rectGate)

## ----Norm2Filter, echo=TRUE, results='markup'----------------------------
morphGate <- norm2Filter("FSC-H", "SSC-H", filterId="MorphologyGate", 
                         scale=2)
smaller <- Subset(fs, morphGate)
fs[[1]]
smaller[[1]]

## ----Split, echo=TRUE, results='markup'----------------------------------
split(smaller[[1]], kmeansFilter("FSC-H"=c("Low","Medium","High"),
                                 filterId="myKMeans"))

## ----Split2, echo=TRUE, results='markup'---------------------------------
split(smaller, kmeansFilter("FSC-H"=c("Low", "Medium", "High"),
                            filterId="myKMeans"))

## ----CombineFilter, echo=TRUE, results='markup'--------------------------

rectGate & morphGate
rectGate | morphGate
!morphGate


## ----Summary5, echo=TRUE, results='markup'-------------------------------
summary(filter(smaller[[1]],rectGate %&% morphGate))

## ----Transfo5, echo=TRUE,results='markup'--------------------------------
tFilter <- transform("FL1-H"=log,"FL2-H"=log)
tFilter

## ----TectGate3, echo=TRUE, results='markup'------------------------------
rect2 <- rectangleGate(filterId="Another Rect", "FL1-H"=c(1,2), 
"FL2-H"=c(2,3)) %on% tFilter
rect2

## ----Plot6,echo=TRUE, results='hide'-------------------------------------
autoplot(tFilter %on% smaller[[1]], "FL1-H","FL2-H")

## ----loadData------------------------------------------------------------
library(flowWorkspace)
fcsfiles <- list.files(pattern = "CytoTrol"
                       , system.file("extdata", package = "flowWorkspaceData")
                       , full = TRUE)
fs <- read.flowSet(fcsfiles)


## ----createGatingSet-----------------------------------------------------
gs <- GatingSet(fs)
gs

## ----getComp, echo=FALSE-------------------------------------------------
gs_manual <- load_gs(list.files(pattern = "gs_manual"
                                , system.file("extdata", package = "flowWorkspaceData")
                                , full = TRUE))
comp <- getCompensationMatrices(gs_manual[[1]])

## ----Compensate----------------------------------------------------------
comp
gs <- compensate(gs, comp)

## ----plotComp,echo=TRUE, results='hide'----------------------------------
fs_comp <- getData(gs)
transList <- estimateLogicle(fs[[1]], c("V545-A","V450-A"))
library(gridExtra)
p1 <- autoplot(transform(fs[[1]], transList)
               , 'V545-A', 'V450-A') + ggtitle("Before")
p2 <- autoplot(transform(fs_comp[[1]], transList)
               , 'V545-A', 'V450-A') + ggtitle("After")
grid.arrange(as.ggplot(p1), as.ggplot(p2), ncol = 2)

## ----nodes---------------------------------------------------------------
getNodes(gs)

## ----addTrans------------------------------------------------------------
biexpTrans <- flowJo_biexp_trans(channelRange=4096, maxValue=262144
                          , pos=4.5,neg=0, widthBasis=-10)
chnls <- parameters(comp)
tf <- transformerList(chnls, biexpTrans)

#or use estimateLogicle directly on GatingHierarchy object to generate transformerList automatically
#tf <- estimateLogicle(gs[[1]], chnls)

gs <- transform(gs, tf)

## ----plotTrans,echo=TRUE, results='hide'---------------------------------
p1 <- autoplot(fs_comp[[1]], "B710-A") + ggtitle("raw")
p2 <- autoplot(flowData(gs)[[1]], "B710-A") + 
          ggtitle("trans") + 
          ggcyto_par_set(limits = "instrument")
grid.arrange(as.ggplot(p1), as.ggplot(p2), ncol = 2)

## ----addGate-nonDebris---------------------------------------------------
rg1 <- rectangleGate("FSC-A"=c(50000, Inf), filterId="NonDebris")
add(gs, rg1, parent = "root")
getNodes(gs)
# gate the data
recompute(gs)

## ----plotGate,echo=TRUE, results='hide'----------------------------------
autoplot(gs, "NonDebris")

## ----plotGate-density,echo=TRUE, results='hide'--------------------------
ggcyto(gs, aes(x = `FSC-A`)) + geom_density() + geom_gate("NonDebris")

## ----getStats1-----------------------------------------------------------
getTotal(gs[[1]], "NonDebris")#counts
getProp(gs[[1]], "NonDebris")#proportion

## ----addGate-singlets----------------------------------------------------
# add the second gate
mat <- matrix(c(54272,59392,259071.99382782
                ,255999.994277954,62464,43008,70656
                ,234495.997428894,169983.997344971,34816)
              , nrow = 5)
colnames(mat) <-c("FSC-A", "FSC-H")
mat
pg <- polygonGate(mat)
add(gs, pg, parent = "NonDebris", name = "singlets")

# add the third gate
rg2 <- rectangleGate("V450-A"=c(2000, Inf))
add(gs, rg2, parent = "singlets", name = "CD3")
getNodes(gs)


## ----addQuadGate---------------------------------------------------------
qg <- quadGate("B710-A" = 2000, "R780-A" = 3000)
add(gs, qg, parent="CD3", names = c("CD8", "DPT", "CD4", "DNT"))
getChildren(gs[[1]], "CD3")
# gate the data from "singlets"
recompute(gs, "singlets")

## ----plotgs, eval=FALSE--------------------------------------------------
#  plot(gs)

## ----plotwfdo, echo=FALSE, results='hide'--------------------------------
if(suppressWarnings(require(Rgraphviz))){
    plot(gs)
}else{
    plot(1,1, type="n", axes=FALSE, ann=FALSE)
    text(1,1,"Need to install Rgraphviz")
}

## ----plotGateAll, results='hide'-----------------------------------------
autoplot(gs[[1]])

## ----getData-------------------------------------------------------------
fs_nonDebris <- getData(gs, "NonDebris")
fs_nonDebris 
nrow(fs_nonDebris[[1]])
nrow(fs[[1]])

## ----getStats2-----------------------------------------------------------
getPopStats(gs)

## ----Rm------------------------------------------------------------------
Rm('CD3', gs)
getNodes(gs)
Rm('NonDebris', gs)
getNodes(gs)

## ----openCyto-nonDebris--------------------------------------------------
library(openCyto)
thisData <- getData(gs)
nonDebris_gate <- fsApply(thisData
                          , function(fr)
                            openCyto:::.mindensity(fr, channels = "FSC-A"))
add(gs, nonDebris_gate, parent = "root", name = "nonDebris")
recompute(gs)

## ----openCyto-singletGate------------------------------------------------
thisData <- getData(gs, "nonDebris") #get parent data
singlet_gate <- fsApply(thisData
                        , function(fr)
                          openCyto:::.singletGate(fr, channels =c("FSC-A", "FSC-H")))
add(gs, singlet_gate, parent = "nonDebris", name = "singlets")
recompute(gs)

## ----openCyto-CD3--------------------------------------------------------
thisData <- getData(gs, "singlets") #get parent data
CD3_gate <- fsApply(thisData
                    , function(fr)
                        openCyto:::.mindensity(fr, channels ="V450-A"))
add(gs, CD3_gate, parent = "singlets", name = "CD3")
recompute(gs)

## ----openCyto-Tsub-------------------------------------------------------
thisData <- getData(gs, "CD3") #get parent data
Tsub_gate <- fsApply(thisData
                     , function(fr)
                        openCyto::quadGate.seq(fr
                             , channels = c("B710-A", "R780-A")
                             , gFunc = 'mindensity'
                            )
                  )
add(gs, Tsub_gate, parent = "CD3", names = c("CD8", "DPT", "CD4", "DNT"))
recompute(gs)

## ----plotALL-openCyto, results='hide'------------------------------------
autoplot(gs[[1]])

