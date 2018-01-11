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


makeGT <- function (TagNames) {
  #collapse fcs files for ROC analysis
  
  collapseBy <- readline(prompt="Number of files to collapse (blank = don't collapse): ")
  if (collapseBy == "") {
    collapseBool <- "FALSE"
    groupByNum <- NA
  } else {
    collapseBool <- "TRUE"
    groupByNum <- collapseBy
  }
  
  eventThreshold <- readline(prompt="Ignore peaks below N events (blank = 20): ")
  if (eventThreshold == "") eventThreshold <- "20"

  
### automatically make gating template based on FTs in index file
# function for adding rows to gating template data frame to reduce syntax
addGtRow <- function(df, alias=NA, pop=NA, parent=NA, dims=NA, gmethod=NA, gargs=NA, collapse=collapseBool, group=groupByNum, ppmethod=NA, ppargs=NA) {
	row <- list(alias=alias
			, pop=pop
			, parent=parent
			, dims=dims
			, gating_method=gmethod
			, gating_args=gargs
			, collapseDataForGating=collapse
			, groupBy=group
			, preprocessing_method=ppmethod
			, preprocessing_args=ppargs
			)
	df[nrow(df)+1, names(row)] <- row
	return(df)
}

#detect which CellTags are present in primary and secondary quadrants
Q1 <- c("FT1","FT2","FT3","FT4"
	,"FT6","FT7","FT8","FT9"
	,"FT11","FT12","FT13","FT16","FT17")
Q1 <- Q1 %in% TagNames

Q2 <- c("FT20","FT21", "FT22") %in% TagNames
Q3 <- c("FT19","FT0") %in% TagNames
Q4 <- c("FT5","FT10", "FT15") %in% TagNames

#check that there are celltags
if(sum(Q1,Q2,Q3,Q4)<1) stop("No recognized CellTags in index file!")

#check for FT0 v FT19 
if(sum(Q3)>1) stop("FT0 and FT19 should not be used at the same time")
Q3isFT19 <- "FT19" %in% TagNames

q1 <- c("FT1","FT2","FT3"
	,"FT6","FT7","FT8"
	,"FT11","FT12")
q1 <- q1 %in% TagNames

q2 <- c("FT16","FT17") %in% TagNames
q3 <- c("FT13") %in% TagNames
q4 <- c("FT4","FT9") %in% TagNames


#ask for scatter and fluorescent cutoffs, default by no answer
fscMIN <- readline(prompt="Enter min FSC.A cutoff (default=6e4): ")
fscMAX <- readline(prompt="Enter max FSC.A cutoff (default=2.6e5): ")
sscMIN <- readline(prompt="Enter min SSC.A cutoff (default=1.5e4): ")
sscMAX <- readline(prompt="Enter max SSC.A cutoff (default=2.6e5): ")

G.nMIN <- readline(prompt="Enter min GFP cutoff in logicle scale (default=0.5): ")
G.nMAX <- readline(prompt="Enter max GFP cutoff in logicle scale (default=4.5): ")
R.nMIN <- readline(prompt="Enter min RFP cutoff in logicle scale (default=0.4): ")
R.nMAX <- readline(prompt="Enter max RFP cutoff in logicle scale (default=4.5): ")

if(fscMIN == "") fscMIN <- "6e4"
if(fscMAX == "") fscMAX <- "2.6e5"
if(sscMIN == "") sscMIN <- "1.5e4"
if(sscMAX == "") sscMAX <- "2.6e5"

if(G.nMIN == "") G.nMIN <- "0.5"
if(G.nMAX == "") G.nMAX <- "4.5"
if(R.nMIN == "") R.nMIN <- "0.4"
if(R.nMAX == "") R.nMAX <- "4.5"


#initialize gating template data frame
df <- data.frame(alias=character(),
                 pop=character(),
                 parent=character(),
                 dims=character(),
                 gating_method=character(),
                 gating_args=character(),
                 collapseDataForGating=character(),
                 groupBy=character(),
                 preprocessing_method=character(),
                 preprocessing_args=character(),
                 stringsAsFactors=FALSE)



#add scatter and fluor min max gates
scGateArgs <- paste("min = c("
			,fscMIN
			,","
			,sscMIN
			,"), max=c("
			,fscMAX
			,","
			,sscMAX
			,")"
			, sep="")

flGateArgs <- paste("min = c("
			,G.nMIN
			,","
			,R.nMIN
			,"), max=c("
			,G.nMAX
			,","
			,R.nMAX
			,")"
			, sep="")



df <- addGtRow(df,"nonMaxS", "nonMaxS+", "root",  "F.A,S.A", "boundary", scGateArgs)
df <- addGtRow(df,"nonMaxF", "nonMaxF+", "nonMaxS",  "G.n,R.n", "boundary", flGateArgs)


peakN <- NULL


### primary quadrants (Q) ###
#G.n Q ref split selecting RIGHT
alias <- "Q1Q2"
pop <- "G.n+"
parent <- "nonMaxF"
dims <- "G.n"
if (sum(Q1,Q4)>0 & sum(Q2,Q3)>0) {
	#calculate peakN, only options are 2, 3 
	#hedge on good separation between [Q2 Q3] : [q2 q3 FT15]
	#that regardless of potential peak split within [q1 q4 FT5 FT10]
	#assume [Q2 Q3] : [q2 q3 FT15] split dominates over splits internal to [q1 q4 FT5 FT10]
	#exception is when FT5 FT6 FT7 FT8 FT9 are missing
	
	#there are three independent lumps
	peakN <- sum(TRUE  # represents Q2 Q3 which has been checked
			, sum(q2, q3, "FT15" %in% TagNames)>0
			, sum(c(q2, q3, "FT5", "FT10") %in% TagNames)>0
			)
	
	FT.5.6.7.8.9 <- isTRUE(sum(c("FT5", "FT6", "FT7", "FT8", "FT9") %in% TagNames)>0)
	FT.1.2.3.4 <- isTRUE(sum(c("FT1", "FT2", "FT3", "FT4") %in% TagNames)>0)
	FT.11.12.10 <- isTRUE(sum(c("FT11", "FT12", "FT10") %in% TagNames)>0)
	
	if (!FT.5.6.7.8.9 && (FT.1.2.3.4 && FT.11.12.10)) peakN <- peakN + 1
	
	df <- addGtRow(df, alias, pop, parent,  dims, "specificMindensity", paste(
			"selectedPeak=1, totalPeaks="
			,peakN
			,", eventThreshold="
			,eventThreshold
			, sep=""
			)
	)
} else if (sum(Q1,Q4)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate")
} else if (sum(Q2,Q3)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate", "minAtmin=F")
}

#G.n Q ref split selecting LEFT, just copy above, change key names
alias <- "Q1Q2neg"
pop <- "G.n-"

tempRow <- df[nrow(df),]
tempRow$alias <- alias
tempRow$pop <- pop
df[nrow(df)+1, names(tempRow)] <- tempRow


#R.n Q ref split RIGHT of Q1Q2 split
alias <- "Q1Q4r"
pop <- "R.n+"
parent <- "Q1Q2"
dims <- "R.n"
if (sum(Q1)>0 & sum(Q4)>0) {
	#calculate peakN, only options are 2, 3, or 4 

	#there are three independent lumps
	peakN <- sum(TRUE  # represents Q4 which has been checked
			, sum(q3, q4)>0
			, sum(c("FT2", "FT7", "FT12", "FT17", "FT3", "FT8") %in% TagNames)>0
			, sum(c("FT1", "FT6", "FT11", "FT16") %in% TagNames)>0
			)

	df <- addGtRow(df, alias, pop, parent,  dims, "specificMindensity", paste(
			"selectedPeak=1, totalPeaks="
			,peakN
			,", eventThreshold="
			,eventThreshold
			, sep=""
			)
	)
} else if (sum(Q1)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate")
} else if (sum(Q4)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate", "minAtmin=F")
}

#R.n Q ref split LEFT of Q1Q2 split
alias <- "Q1Q4l"
pop <- "R.n+"
parent <- "Q1Q2neg"
dims <- "R.n"
if (sum(Q2)>0 & sum(Q3)>0) {
	#calculate peakN, only options are 2, 3, or 4 for the exact number of FTs in Q2 Q3
	peakN <- sum(Q2,Q3)
	
	df <- addGtRow(df, alias, pop, parent,  dims, "specificMindensity", paste(
			"selectedPeak=1, totalPeaks="
			,peakN
			,", eventThreshold="
			,eventThreshold
			, sep=""
			)
	)
} else if (sum(Q2)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate")
} else if (sum(Q3)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate", "minAtmin=F")
}


#Q gates
parent <- "nonMaxF"
dims <- "G.n,R.n"
method <- "refGate"
#args <- "Q1Q4:Q1Q2"
#Q1 gate
if (sum(Q1)>0) df <- addGtRow(df, "Q1", "G.n+R.n+", parent, dims, method, "Q1Q4r:Q1Q2")
#Q2 gate
if (sum(Q2)>0) df <- addGtRow(df, "Q2", "G.n-R.n+", parent, dims, method, "Q1Q4l:Q1Q2")
#Q3 gate
if (sum(Q3)>0) df <- addGtRow(df, "Q3", "G.n-R.n-", parent, dims, method, "Q1Q4l:Q1Q2")
#Q4 gate
if (sum(Q4)>0) df <- addGtRow(df, "Q4", "G.n+R.n-", parent, dims, method, "Q1Q4r:Q1Q2")


### secondary quadrants in Q1 (q) ###
#gate extreme wings first FT4 FT16
FT.4 <- "FT4" %in% TagNames
FT.3.9 <- isTRUE(sum(c("FT3","FT9") %in% TagNames) > 0)
FT.2.8 <- isTRUE(sum(c("FT2","FT8") %in% TagNames) > 0)
FT.1.7.13 <- isTRUE(sum(c("FT1","FT7", "FT13") %in% TagNames) > 0)
FT.6.12 <- isTRUE(sum(c("FT6","FT12") %in% TagNames) > 0)
FT.11.17 <- isTRUE(sum(c("FT11","FT17") %in% TagNames) > 0)
FT.16 <- "FT16" %in% TagNames

#gate Q1-FT4
hi <- isTRUE(sum(FT.3.9, FT.2.8, FT.1.7.13, FT.6.12, FT.11.17, FT.16)>0)
peakN <- sum(FT.4, FT.3.9, FT.2.8, FT.1.7.13, FT.6.12, FT.11.17, FT.16)
if (FT.4 & hi) {
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "FT4p", "fRatio-", "Q1", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "Q1-FT4", "fRatio+", "Q1", "fRatio", "specificMindensity", args)
} else if (FT.4) { 
	df <- addGtRow(df, "FT4p", "fRatio+", "Q1", "fRatio", "dummySubgate") 
} else if (hi) {
	df <- addGtRow(df, "Q1-FT4", "fRatio+", "Q1", "fRatio", "dummySubgate")
}

#gate Q1-FT4-FT16
low <- isTRUE(sum(FT.3.9, FT.2.8, FT.1.7.13, FT.6.12, FT.11.17)>0)
peakN <- sum(FT.3.9, FT.2.8, FT.1.7.13, FT.6.12, FT.11.17, FT.16)
if (FT.16 & low) {
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, ", selectFromRight=TRUE", sep="")
	df <- addGtRow(df, "FT16p", "fRatio+", "Q1-FT4", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "Q1-FT4-FT16", "fRatio-", "Q1-FT4", "fRatio", "specificMindensity", args)
} else if (FT.16) { 
	df <- addGtRow(df, "FT16p", "fRatio+", "Q1-FT4", "fRatio", "dummySubgate") 
} else if (low) {
	df <- addGtRow(df, "Q1-FT4-FT16", "fRatio+", "Q1-FT4", "fRatio", "dummySubgate")
}





#G.n q ref split
alias <- "q1q2"
pop <- "G.n+"
parent <- "Q1-FT4-FT16"
dims <- "G.n"

FT.9 <- "FT9" %in% TagNames
FT.17 <- "FT17" %in% TagNames

if (sum(q1,FT.9)>0 & sum(FT.17,q3)>0) {
	#calculate peakN, only options are 2, 3 or 4

	#there are three independent lumps
	peakN <- sum(TRUE  # represents FT.17 q3 which has been checked
			, sum(c("FT6", "FT7", "FT8", "FT9") %in% TagNames)>0
			, sum(c("FT11", "FT12") %in% TagNames)>0
			, sum(c("FT1", "FT2", "FT3") %in% TagNames)>0
			)

	df <- addGtRow(df, alias, pop, parent,  dims, "specificMindensity", paste(
			"selectedPeak=1, totalPeaks="
			,peakN
			,", eventThreshold="
			,eventThreshold
			, sep=""
			)
	)
} else if (sum(q1,FT.9)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate")
} else if (sum(FT.17,q3)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate", "minAtmin=F")
}


#R.n q ref split
alias <- "q1q4"
pop <- "R.n+"
parent <- "Q1-FT4-FT16"
dims <- "R.n"

if (sum(q1,FT.17)>0 & sum(q3,FT.9)>0) {
	#calculate peakN, only options are 2, 3 or 4

	#there are three independent lumps
	peakN <- sum(TRUE  # represents q3 FT.9 which has been checked
			, sum(c("FT2", "FT7", "FT12", "FT17") %in% TagNames)>0
			, sum(c("FT3", "FT8") %in% TagNames)>0
			, sum(c("FT1", "FT6", "FT11") %in% TagNames)>0
			)

	df <- addGtRow(df, alias, pop, parent,  dims, "specificMindensity", paste(
			"selectedPeak=1, totalPeaks="
			,peakN
			,", eventThreshold="
			,eventThreshold
			, sep=""
			)
	)
} else if (sum(q1,FT.17)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate")
} else if (sum(q3,FT.9)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate", "minAtmin=F")
}

#q gates
parent <- "Q1-FT4-FT16"
dims <- "G.n,R.n"
method <- "refGate"
args <- "q1q4:q1q2"
#q1 gate
if (sum(q1)>0) df <- addGtRow(df, "q1", "G.n+R.n+", parent, dims, method, args)
#q2 gate
if (sum(q2)>0) df <- addGtRow(df, "q2", "G.n-R.n+", parent, dims, method, args)
#q3 gate
if (sum(q3)>0) df <- addGtRow(df, "q3", "G.n-R.n-", parent, dims, method, args)
#q4 gate
if (sum(q4)>0) df <- addGtRow(df, "q4", "G.n+R.n-", parent, dims, method, args)


### subgates in Q4, Q2 ###
## gates in Q4 ##
low <- "FT15" %in% TagNames
hi <- isTRUE(sum(c("FT5","FT10") %in% TagNames) > 0)
if (low & hi) {
	peakN <- sum(c("FT5","FT10", "FT15") %in% TagNames)
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "FT15p", "G.n-", "Q4", "G.n", "specificMindensity", args)
	df <- addGtRow(df, "FT.5.10", "G.n+", "Q4", "G.n", "specificMindensity", args)
} else if (low) { 
	df <- addGtRow(df, "FT15p", "G.n+", "Q4", "G.n", "dummySubgate") 
} else if (hi) { 
	df <- addGtRow(df, "FT.5.10", "G.n+", "Q4", "G.n", "dummySubgate")
}

low <- "FT10" %in% TagNames
hi <- "FT5" %in% TagNames
if (low & hi) {
	df <- addGtRow(df, "FT10p", "G.n-", "FT.5.10", "G.n", "mindensity")
	df <- addGtRow(df, "FT5p", "G.n+", "FT.5.10", "G.n", "mindensity")
} else if (low) {
	df <- addGtRow(df, "FT10p", "G.n+", "FT.5.10", "G.n", "dummySubgate")
} else if (hi) { 
	df <- addGtRow(df, "FT5p", "G.n+", "FT.5.10", "G.n", "dummySubgate")
}


## gates in Q2 ##
low <- "FT22" %in% TagNames
hi <- isTRUE(sum(c("FT20","FT21") %in% TagNames) > 0)
if (low & hi) {
	peakN <- sum(c("FT20","FT21", "FT22") %in% TagNames)
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "FT22p", "R.n-", "Q2", "R.n", "specificMindensity", args)
	df <- addGtRow(df, "FT.20.21", "R.n+", "Q2", "R.n", "specificMindensity", args)
} else if (low) {
	df <- addGtRow(df, "FT22p", "R.n+", "Q2", "R.n", "dummySubgate")
} else if (hi) {
	df <- addGtRow(df, "FT.20.21", "R.n+", "Q2", "R.n", "dummySubgate")
}

low <- "FT21" %in% TagNames
hi <- "FT20" %in% TagNames
if (low & hi) {
	df <- addGtRow(df, "FT21p", "R.n-", "FT.20.21", "R.n", "mindensity")
	df <- addGtRow(df, "FT20p", "R.n+", "FT.20.21", "R.n", "mindensity")
} else if (low) {
	df <- addGtRow(df, "FT21p", "R.n+", "FT.20.21", "R.n", "dummySubgate")
} else if (hi) {
	df <- addGtRow(df, "FT20p", "R.n+", "FT.20.21", "R.n", "dummySubgate")
}


### subgates in q4, q2 ###
# these are just q4 == FT9p, q2 == FT17p


### subgates in q1 ###
## gate q1 low wing ##
FT.3 <- "FT3" %in% TagNames
FT.2.8 <- isTRUE(sum(c("FT2","FT8") %in% TagNames) > 0)
FT.1.7 <- isTRUE(sum(c("FT1","FT7") %in% TagNames) > 0)
FT.6.12 <- isTRUE(sum(c("FT6","FT12") %in% TagNames) > 0)
FT.11 <- "FT11" %in% TagNames
hi <- isTRUE(sum(FT.1.7, FT.6.12, FT.11)>0)

peakN <- sum(FT.3, FT.2.8, FT.1.7, FT.6.12, FT.11)

if ((FT.3 | FT.2.8) & hi) {
	if (FT.3 & FT.2.8) { 
		args <- paste("selectedPeak=2, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	} else {
		args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	}
	df <- addGtRow(df, "FT.2.3.8", "fRatio-", "q1", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "FT.1.6.7.11.12", "fRatio+", "q1", "fRatio", "specificMindensity", args)
} else if (FT.3 | FT.2.8) { 
	df <- addGtRow(df, "FT.2.3.8", "fRatio+", "q1", "fRatio", "dummySubgate") 
} else if (hi) { 
	df <- addGtRow(df, "FT.1.6.7.11.12", "fRatio+", "q1", "fRatio", "dummySubgate")
}


# subgates in q1 low wing #
if (FT.3 & FT.2.8) {
	args <- paste("selectedPeak=1, totalPeaks=2", ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "FT3p", "fRatio-", "FT.2.3.8", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "FT.2.8", "fRatio+", "FT.2.3.8", "fRatio", "specificMindensity", args)
} else if (FT.3) { 
	df <- addGtRow(df, "FT3p", "fRatio+", "FT.2.3.8", "fRatio", "dummySubgate") 
} else if (FT.2.8) { 
	df <- addGtRow(df, "FT.2.8", "fRatio+", "FT.2.3.8", "fRatio", "dummySubgate")
}

low <- "FT8" %in% TagNames
hi <- "FT2" %in% TagNames
if (low & hi) {
	df <- addGtRow(df, "FT8p", "fTotal-", "FT.2.8", "fTotal", "mindensity")
	df <- addGtRow(df, "FT2p", "fTotal+", "FT.2.8", "fTotal", "mindensity")
} else if (low) { 
	df <- addGtRow(df, "FT8p", "fTotal+", "FT.2.8", "fTotal", "dummySubgate")
} else if (hi) { 
	df <- addGtRow(df, "FT2p", "fTotal+", "FT.2.8", "fTotal", "dummySubgate")
}

## gate q1 high wing ##
hi <- isTRUE(sum(FT.6.12, FT.11) > 0)

if (FT.1.7 & hi) {
	peakN <- sum(FT.1.7, FT.6.12, FT.11)
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "FT.1.7", "fRatio-", "FT.1.6.7.11.12", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "FT.6.11.12", "fRatio+", "FT.1.6.7.11.12", "fRatio", "specificMindensity", args)
} else if (FT.1.7) { 
	df <- addGtRow(df, "FT.1.7", "fRatio+", "FT.1.6.7.11.12", "fRatio", "dummySubgate") 
} else if (hi) { 
	df <- addGtRow(df, "FT.6.11.12", "fRatio+", "FT.1.6.7.11.12", "fRatio", "dummySubgate")
}

# subgates in q1 high wing #
if (FT.6.12 & FT.11) {
	args <- paste("selectedPeak=1, totalPeaks=2", ", eventThreshold=", eventThreshold, ", selectFromRight=TRUE", sep="")
	df <- addGtRow(df, "FT.6.12", "fRatio-", "FT.6.11.12", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "FT11p", "fRatio+", "FT.6.11.12", "fRatio", "specificMindensity", args)
} else if (FT.6.12) { 
	df <- addGtRow(df, "FT.6.12", "fRatio+", "FT.6.11.12", "fRatio", "dummySubgate")
} else if (FT.11) { 
	df <- addGtRow(df, "FT11p", "fRatio+", "FT.6.11.12", "fRatio", "dummySubgate")
}


low <- "FT12" %in% TagNames
hi <- "FT6" %in% TagNames
if (low & hi) {
	df <- addGtRow(df, "FT12p", "fTotal-", "FT.6.12", "fTotal", "mindensity")
	df <- addGtRow(df, "FT6p", "fTotal+", "FT.6.12", "fTotal", "mindensity")
} else if (low) { 
	df <- addGtRow(df, "FT12p", "fTotal+", "FT.6.12", "fTotal", "dummySubgate")
} else if (hi) { 
	df <- addGtRow(df, "FT6p", "fTotal+", "FT.6.12", "fTotal", "dummySubgate")
}

## subgates in q1 mid wing ##
low <- "FT7" %in% TagNames
hi <- "FT1" %in% TagNames
if (low & hi) {
	df <- addGtRow(df, "FT7p", "fTotal-", "FT.1.7", "fTotal", "mindensity")
	df <- addGtRow(df, "FT1p", "fTotal+", "FT.1.7", "fTotal", "mindensity")
} else if (low) { 
	df <- addGtRow(df, "FT7p", "fTotal+", "FT.1.7", "fTotal", "dummySubgate")
} else if (hi) { 
	df <- addGtRow(df, "FT1p", "fTotal+", "FT.1.7", "fTotal", "dummySubgate")
}



### all final population gates ###
dims <- "G.n,R.n"
method <- "flowClust"
args <- "K=1, quantile=0.8"

if ("FT1" %in% TagNames) df <- addGtRow(df, "FT1", "FT1+", "FT1p", dims, method, args)
if ("FT2" %in% TagNames) df <- addGtRow(df, "FT2", "FT2+", "FT2p", dims, method, args)
if ("FT3" %in% TagNames) df <- addGtRow(df, "FT3", "FT3+", "FT3p", dims, method, args)
if ("FT4" %in% TagNames) df <- addGtRow(df, "FT4", "FT4+", "FT4p", dims, method, args)
if ("FT5" %in% TagNames) df <- addGtRow(df, "FT5", "FT5+", "FT5p", dims, method, args)
if ("FT6" %in% TagNames) df <- addGtRow(df, "FT6", "FT6+", "FT6p", dims, method, args)
if ("FT7" %in% TagNames) df <- addGtRow(df, "FT7", "FT7+", "FT7p", dims, method, args)
if ("FT8" %in% TagNames) df <- addGtRow(df, "FT8", "FT8+", "FT8p", dims, method, args)
if ("FT9" %in% TagNames) df <- addGtRow(df, "FT9", "FT9+", "q4", dims, method, args)
if ("FT10" %in% TagNames) df <- addGtRow(df, "FT10", "FT10+", "FT10p", dims, method, args)
if ("FT11" %in% TagNames) df <- addGtRow(df, "FT11", "FT11+", "FT11p", dims, method, args)
if ("FT12" %in% TagNames) df <- addGtRow(df, "FT12", "FT12+", "FT12p", dims, method, args)
if ("FT13" %in% TagNames) df <- addGtRow(df, "FT13", "FT13+", "q3", dims, method, args)
if ("FT15" %in% TagNames) df <- addGtRow(df, "FT15", "FT15+", "FT15p", dims, method, args)
if ("FT16" %in% TagNames) df <- addGtRow(df, "FT16", "FT16+", "FT16p", dims, method, args)
if ("FT17" %in% TagNames) df <- addGtRow(df, "FT17", "FT17+", "q2", dims, method, args)
if ("FT19" %in% TagNames) df <- addGtRow(df, "FT19", "FT19+", "Q3", dims, method, args)
if ("FT0" %in% TagNames) df <- addGtRow(df, "FT0", "FT0+", "Q3", dims, method, args)
if ("FT20" %in% TagNames) df <- addGtRow(df, "FT20", "FT20+", "FT20p", dims, method, args)
if ("FT21" %in% TagNames) df <- addGtRow(df, "FT21", "FT21+", "FT21p", dims, method, args)
if ("FT22" %in% TagNames) df <- addGtRow(df, "FT22", "FT22+", "FT22p", dims, method, args)

res <- list("gt"=df
	, "G.nRange"=c(as.numeric(G.nMIN), as.numeric(G.nMAX))
	, "R.nRange"=c(as.numeric(R.nMIN), as.numeric(R.nMAX))
)
return(res)

}

