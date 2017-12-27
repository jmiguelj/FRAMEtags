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
countCTs <- function(gs, CTindex, Tubeindex, well=FALSE) {

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


#extract CT pop counts
s <- getPopStats(gs, statistic="count", format="wide")
CTcounts <- s[rownames(s) %in% CTindex$CellTag, , drop=F]
colnames(CTcounts) <- pData(gs)$Tube


#find matching names in CTindex and Tubeindex, relabel CTcounts and s if found

CTcountsRowNames <- sapply(rownames(CTcounts), function(X) {
					index <- which(CTindex[,"CellTag"] == X)
					if (length(index) == 1) { CTindex[index, "Strain"]
					} else { X }
				})

CTcountsColNames <- sapply(colnames(CTcounts), function(X) {
					index <- which(Tubeindex[,"Tube"] == X)
					if (length(index) == 1) { Tubeindex[index, "Treatment"]
					} else { X }
				})


rownames(CTcounts) <- CTcountsRowNames
colnames(CTcounts) <- CTcountsColNames


#reorder alphabetically
CTcounts <- CTcounts[order(rownames(CTcounts)), order(colnames(CTcounts)), drop=F]
s <- s[, order(colnames(s)), drop=F]

parents <- s[rownames(s) == "nonMaxF", , drop=F]

CTpercent <- t(t(CTcounts)/colSums(CTcounts))


CTcounts <- rbind(CTcounts, parents)

res <- list("fullStats" = s, "counts" = CTcounts, "percents" = CTpercent, "pData" = pData(gs))

return(res)
}








