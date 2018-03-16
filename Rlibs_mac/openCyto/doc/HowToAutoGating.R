## ---- echo=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(message = FALSE, warning = FALSE, fig.height= 2, fig.width= 3, fig.show = "hold")

## ------------------------------------------------------------------------
library(openCyto)
library(ggcyto)

gs <- load_gs(system.file("extdata/gs_bcell_auto", package = "flowWorkspaceData"))

## ------------------------------------------------------------------------
fr <- getData(gs[[2]], "Live")
chnl <- "CD3"
g <- openCyto:::.mindensity(fr, channels = chnl)
autoplot(fr, chnl) + geom_gate(g)
autoplot(fr, chnl, "SSC-A") + geom_gate(g)

## ------------------------------------------------------------------------
fr <- getData(gs[[1]], "boundary")
chnl <- "FSC-A"
g <- openCyto:::.mindensity(fr, channels = chnl)
mylimits <- ggcyto_par_set(limits = "instrument")
p <- autoplot(fr, chnl) + mylimits
p + geom_gate(g)
autoplot(fr, chnl, "SSC-A") + geom_gate(g)

## ------------------------------------------------------------------------
g <- openCyto:::.mindensity(fr, channels = chnl, gate_range=c(7e4,1e5), adjust = 1.5)
p + geom_gate(g)
autoplot(fr, chnl, "SSC-A") + geom_gate(g)

## ------------------------------------------------------------------------
g <- openCyto:::.mindensity(fr, channels = chnl, min = 7e4, max = 1e5)
p + geom_gate(g)

## ------------------------------------------------------------------------
fr <- getData(gs[[1]], "lymph")
chnl <- "Live"
g <- openCyto:::.tailgate(fr, channels = chnl, tol = 0.05)
p <- autoplot(fr, chnl) + mylimits
p + geom_gate(g)
autoplot(fr, chnl, "SSC-A") + geom_gate(g)

## ------------------------------------------------------------------------
g <- openCyto:::.quantileGate(fr, channels = chnl, probs = 0.99)
p <- autoplot(fr, chnl) + mylimits
p + geom_gate(g)
autoplot(fr, chnl, "SSC-A") + geom_gate(g)

## ------------------------------------------------------------------------
fr <- getData(gs[[1]], "root")
chnl <- c("FSC-A", "SSC-A")
g <- openCyto:::.boundary(fr, channels = chnl, min = c(0, 0), max=c(2.5e5,2.5e5))
p <- autoplot(fr, x = chnl[1], y = chnl[2])
p + geom_gate(g)

## ------------------------------------------------------------------------
fr <- read.FCS(system.file("extdata/CytoTrol_CytoTrol_1.fcs", package = "flowWorkspaceData"))
chnl <- c("FSC-A", "FSC-H")
g <- openCyto:::.singletGate(fr, channels = chnl)
p <- autoplot(fr, x = chnl[1], y = chnl[2])
p + geom_gate(g)

## ------------------------------------------------------------------------
fr <- getData(gs[[1]], "nonDebris")
chnl <- c("FSC-A", "SSC-A")
g <- openCyto:::.flowClust.2d(fr, channels = chnl, K=2, target=c(1e5,5e4), quantile=0.95)
p <- autoplot(fr, x = chnl[1], y = chnl[2]) + mylimits
p + geom_gate(g)

## ------------------------------------------------------------------------
fr <- getData(gs[[1]], "CD19andCD20")
chnl <- c("CD38", "CD24")
g <- openCyto:::.flowClust.2d(fr, channels = chnl, K=6,transitional=TRUE,target=c(3.5e3,3.5e3), quantile=0.95,translation=0.15, pp_res = NULL)
p <- autoplot(fr, x = chnl[1], y = chnl[2]) + mylimits
p + geom_gate(g)

## ------------------------------------------------------------------------
gs <- load_gs(system.file("extdata/gs_DC_auto", package = "flowWorkspaceData"))
fr <- getData(gs[[2]], "HLADR+")
chnl <- c("CD11c", "CD123")
p <- autoplot(fr, chnl[1], chnl[2])
g <- openCyto:::.quadGate.tmix(fr, channels = chnl, K = 3, usePrior = "no")
p + geom_gate(g)

