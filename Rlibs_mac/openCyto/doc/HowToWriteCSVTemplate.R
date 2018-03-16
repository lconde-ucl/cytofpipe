## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", message = FALSE, warning = FALSE, knitr.table.format = "markdown")
assign("depthtrigger", 60, data.table:::.global)

## ----init, echo=F--------------------------------------------------------
library(openCyto)
library(data.table)

template_row <- data.table(alias = "TH", pop = "+-"
                          , parent = "cd3", dims = "cd4,cd8"
                          , gating_method = "mindensity", gating_args = ""
                          , collapseDataForGating = "TRUE", groupBy = "4"
                          , preprocessing_method = "", preprocessing_args = ""
                      )



## ----A+/-(before), echo = F, results = "asis"----------------------------
	this_row <- copy(template_row)	
	this_row[, pop := "+/-"]
	this_row[, alias := "*"]
	this_row[, dims := "cd4"]
	kable(this_row[, 1:6, with = F])

## ----A+/-(after), echo = F, results = "asis"-----------------------------
      kable(openCyto:::.preprocess_row(this_row)[, 1:6, with = F])
      

## ----A+B+(before), echo = F, results = "asis"----------------------------
	this_row <- copy(template_row)	
	this_row[, pop := "+-"]
	this_row[, alias := "T helper"]
	kable(this_row[, 1:6, with = F])

## ----A+B+(after), echo = F, results = "asis"-----------------------------
      kable(openCyto:::.preprocess_row(this_row)[, 1:6, with = F])

## ----A+/-B+/-(before), echo = F, results = "asis"------------------------
	this_row <- copy(template_row)	
	this_row[, pop := "+/-+/-"]
	this_row[, alias := "*"]
	kable(this_row[, 1:6, with = F])

## ----A+/-B+/-(after), echo = F, results = "asis"-------------------------
      kable(openCyto:::.preprocess_row(this_row)[, 1:6, with = F])

## ----multiPop(before), echo = F, results = "asis"------------------------
	this_row <- copy(template_row)	
	this_row[, pop := "*"]
	this_row[, alias := "CD4,CD8"]
	this_row[, gating_method := "curv2gate"]
	kable(this_row[, 1:6, with = F])

## ----multiPop(after), echo = F, results = "asis"-------------------------
      kable(openCyto:::.preprocess_row(this_row)[, 1:6, with = F])

## ----multiPop_noExpand(before), echo = F, results = "asis"---------------
	this_row[, alias := "*"]
	kable(this_row[, 1:6, with = F])

## ----finalizer, include=FALSE--------------------------------------------
assign("depthtrigger", 3, data.table:::.global)

