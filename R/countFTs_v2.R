#' @import methods
#' @import flowWorkspace
#' @import openCyto
#' @import ggcyto
#' @import parallel
#' @import RColorBrewer
#' @import gtools
#' @import gridExtra
#' @import grid
#' @import methods
#' @import flowCore
#' @import ncdfFlow
#' @import RcppArmadillo
#' @import BH
#' @import ggplot2
countFTs <- function(gs, FTindex, Tubeindex, well=FALSE) {

if (!well) {
	#extract tube number
	n <- strsplit(as.character(keyword(gs, "TUBE NAME")[,1]), "[^0-9]+")
	n <- (unlist(n))
	n <- n[even(1:length(n))]
	n <- as.numeric(n) #strips leading zeros if any

	#add tube number to pData to prepare for renaming
} else {
	n <- keyword(gs, "WELL ID")[,1]
}


pData(gs)$Tube <- as.character(n) 
pData(gs)$time <- fsApply(getData(gs), function(x){max(exprs(x)[,"Time"]) - min(exprs(x)[,"Time"])})


#extract FT pop counts
s <- getPopStats(gs, statistic="count", format="wide")
FTcounts <- s[rownames(s) %in% FTindex$Tag, , drop=F]
colnames(FTcounts) <- pData(gs)$Tube


#find matching names in FTindex and Tubeindex, relabel FTcounts and s if found

FTcountsRowNames <- sapply(rownames(FTcounts), function(X) {
					index <- which(FTindex[,"Tag"] == X)
					if (length(index) == 1) { FTindex[index, "Strain"]
					} else { X }
				})

FTcountsColNames <- sapply(colnames(FTcounts), function(X) {
					index <- which(Tubeindex[,"Tube"] == X)
					if (length(index) == 1) { Tubeindex[index, "Treatment"]
					} else { X }
				})


rownames(FTcounts) <- FTcountsRowNames
colnames(FTcounts) <- FTcountsColNames


#reorder alphabetically
FTcounts <- FTcounts[order(rownames(FTcounts)), order(colnames(FTcounts)), drop=F]
s <- s[, order(colnames(s)), drop=F]

parents <- s[rownames(s) == "nonMaxF", , drop=F]

FTpercent <- t(t(FTcounts)/colSums(FTcounts))


FTcounts <- rbind(FTcounts, parents)

res <- list("fullStats" = s, "counts" = FTcounts, "percents" = FTpercent, "pData" = pData(gs))

return(res)
}








