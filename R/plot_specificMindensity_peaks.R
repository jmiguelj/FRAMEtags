### plot specificMindensity detailed peaks
#library("gridExtra")
#library("grid")
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

specificMindensityPlots <- function(gs, GTindex, whichFlowFrame=1, eventsToPlot=10000) {

#keep only rows that use "specificMindensity"
GTindex <- GTindex[GTindex[,"gating_method"]=="specificMindensity",]

GTline <- rownames(GTindex)
GTindex <- cbind(GTindex, GTline)

#get total gs ranges for overview plots
fr <- getData(gs, "nonMaxF")
fr <- as(fr,"flowFrame")
G.range <- range(exprs(fr)[,"G.n"], finite=T)
R.range <- range(exprs(fr)[,"R.n"], finite=T)
fR.range <- range(exprs(fr)[,"fRatio"], finite=T)
fT.range <- range(exprs(fr)[,"fTotal"], finite=T)

ranges <- list("G.n"=G.range, "R.n"=R.range, "fRatio"=fR.range, "fTotal"=fT.range)

plots <- apply(GTindex, 1, function(X) specificMindensityPlotRow(gs, whichFlowFrame, X["alias"], X["parent"], X["dims"], X["gating_args"], X["collapseDataForGating"], X["groupBy"], X["GTline"], eventsToPlot, ranges))

return(plots)

}

specificMindensityPlotRow <- function(gs, whichFlowFrame, alias, parent, dims, args, collapse, groupBy, GTline, eventsToPlot, ranges) {
	


#which row to analyze
#row <- 1

##manual parameters
#whichFlowFrame <- 1
#alias <- GTindex[row,"alias"]
#parent <- GTindex[row,"parent"]
#dims <- GTindex[row,"dims"]
#args <- GTindex[row,"gating_args"]
#collapse <- GTindex[row,"collapseDataForGating"]
#groupBy <- GTindex[row,"groupBy"]
#GTline <- GTindex[row,"GTline"]

args <- openCyto:::.argParser(args)

if (as.logical(collapse) && !is.na(groupBy)) {
	if (length(gs) < groupBy) groupBy <- length(gs)
		whichFlowFrame <- 1 : groupBy #GTindex[row,"groupBy"]
}


#extract subset flowset representing parent gate data, coerce as single flowframe
fr <- getData(gs[whichFlowFrame], parent)
fr <- as(fr,"flowFrame")

#combine parameters and turn ON detailed result
args$fr <- fr
args$channels <- dims
args$detailedRes <- TRUE


#execute mindensity with appropriate parameters
peakInfo <- do.call(.specificMindensity, args) #has "gate", "peaks" and "adjust"

#peaks is a data frame with at least, data.frame(x=NA, y=NA, valley1X=NA, valley2X=NA)
peaks <- peakInfo$peaks


#get additional data for plotting, only a sample to plot faster
if (eventsToPlot > nrow(fr)) eventsToPlot <- nrow(fr)
fr <- fr[sample(nrow(fr), eventsToPlot)]

frBg <- getData(gs[whichFlowFrame], "nonMaxF")
frBg <- as(frBg,"flowFrame")
if (eventsToPlot > nrow(frBg)) eventsToPlot <- nrow(frBg)
frBg <- frBg[sample(nrow(frBg), eventsToPlot)]

frAlias <- getData(gs[whichFlowFrame], alias)
frAlias <- as(frAlias,"flowFrame")
if (eventsToPlot > nrow(frAlias)) eventsToPlot <- nrow(frAlias)
frAlias <- frAlias[sample(nrow(frAlias), eventsToPlot)]

p <- ggplot(fr, aes_string(x = dims)) 
p <- p + geom_density(fill="black", adjust=peakInfo$args$adjust)
p <- p + geom_gate(peakInfo$gate)


for (i in 1:nrow(peaks)) {
	p <- p + annotate("rect", xmin = peaks[i,"valley1X"], xmax = peaks[i,"valley2X"], ymin = 0, ymax = max(peaks[,"y"]), alpha=0.5, fill=brewer.pal(9,"Set1")[i+1])
	p <- p + annotate("text", x = peaks[i,"x"], y = 1.1*max(peaks[,"y"]), label = paste(i, "\n(", peaks[i,"events"],  ")", sep=""), vjust=1, size=3)
}

xmin <- min(c(peaks[,"valley1X"],peaks[,"valley2X"]))
xmax <- max(c(peaks[,"valley1X"],peaks[,"valley2X"]))
p <- p + annotate("text", x = xmin - (0.1*(xmax-xmin)) , y = 1.1*max(peaks[,"y"]), label = paste("Peak # ->", "(events) ->", sep="\n"), hjust=0, vjust=1, size=3)
p <- p + labs(title = "Peaks found and Gate")

#p


### big picture plots

frBg_xmin <- ranges$G.n[1]
frBg_xmax <- ranges$G.n[2]
frBg_ymin <- ranges$R.n[1]
frBg_ymax <- ranges$R.n[2]

p2 <- ggplot(frBg, aes(x = G.n, y=R.n))
p2 <- p2 + xlim(frBg_xmin,frBg_xmax) + ylim(frBg_ymin,frBg_ymax)
p2 <- p2 + geom_hex(bins=128,  alpha=0.5)
p2 <- p2 + scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), trans = "sqrt")
p2 <- p2 + geom_hex(data=fr, bins=128, col="black")
p2 <- p2 + geom_hex(data=frAlias, bins=128, col="red")
p2 <- p2 + labs(title = "Overview: GFP v. RFP")
p2 <- p2 + annotate("label", x = frBg_xmin, y = frBg_ymax, label = parent, vjust=1, hjust=0, fill="black", col="white", size=5)
p2 <- p2 + annotate("label", x = frBg_xmin, y = frBg_ymax - (0.1*(frBg_ymax - frBg_ymin)), label = alias, vjust=1, hjust=0, fill="red", col="white", size=5)
#p2

frBg_xmin <- ranges$fRatio[1]
frBg_xmax <- ranges$fRatio[2]
frBg_ymin <- ranges$fTotal[1]
frBg_ymax <- ranges$fTotal[2]

frBg_xmin <- frBg_xmin + 0.2*(ranges$fRatio[2] - ranges$fRatio[1])
frBg_ymin <- frBg_ymin + 0.2*(ranges$fTotal[2] - ranges$fTotal[1])


p3 <- ggplot(frBg, aes(x = fRatio, y=fTotal))
p3 <- p3 + xlim(frBg_xmin,frBg_xmax) + ylim(frBg_ymin,frBg_ymax)
p3 <- p3 + geom_hex(bins=128, alpha=0.5)
p3 <- p3 + scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), trans = "sqrt")
p3 <- p3 + geom_hex(data=fr, bins=128, col="black")
p3 <- p3 + geom_hex(data=frAlias, bins=128, col="red")
p3 <- p3 + labs(title = "Overview: RFP/GFP v. RFP+GFP")
p3 <- p3 + annotate("label", x = frBg_xmin, y = frBg_ymax, label = parent, vjust=1, hjust=0, fill="black", col="white", size=5)
p3 <- p3 + annotate("label", x = frBg_xmin, y = frBg_ymax - (0.1*(frBg_ymax - frBg_ymin)), label = alias, vjust=1, hjust=0, fill="red", col="white", size=5)
#p3

title_string <- paste("Plot for Gating Template line "
			, GTline
			, "\nUsed to make gate \'"
			, alias
			, "\' (red)\nfrom parent gate \'"
			, parent
			,"\' (black)"
			, sep="")


args_string <- paste("     Gating parameters:\n       totalPeaks="
			, peakInfo$args$totalPeaks
			, "\n       selectedPeak="
			, peakInfo$args$selectedPeak
			, "\n       eventThreshold="
			, peakInfo$args$eventThreshold
			, "\n       selectFromRight="
			, peakInfo$args$selectFromRight
			, sep="")
direction <- "Left"
if (peakInfo$args$selectFromRight) direction <- "Right"

args_explained <- paste("     Peaks along \'"
			, dims
			,"\' dimension merged into "
			, peakInfo$args$totalPeaks
			, " total peaks\n     Gate made after peak "
			, peakInfo$args$selectedPeak
			, " counting from "
			, direction
			, sep="")

directions_string <- paste("If something looks wrong change line "
			, GTline
			, " in gating_template.csv"
			, sep="")

final_string <- paste(title_string, args_explained, args_string, directions_string, sep="\n\n")
t <- textGrob(final_string, x = unit(0.05, "npc"), y = unit(0.95, "npc"), just=c("left", "top"))


#p <- p + annotate("text", x = xmin - (0.1*(xmax-xmin)) , y = 1.3*max(peaks[,"y"]), label = final_string, hjust=0, vjust=1)
#grid.arrange(t, p, p2, p3, nrow=2, ncol=2)


cat("...plot for gating template line ", GTline, "\n", sep="")
return(list("t"=t, "p"=p, "p2"=p2, "p3"=p3))


}
