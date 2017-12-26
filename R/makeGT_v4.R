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


makeGT <- function (CellTagName) {

### automatically make gating template based on CTs in index file
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
Q1 <- c("CT1","CT2","CT3","CT4"
	,"CT6","CT7","CT8","CT9"
	,"CT11","CT12","CT13","CT16","CT17")
Q1 <- Q1 %in% CellTagName

Q2 <- c("CT20","CT21", "CT22") %in% CellTagName
Q3 <- c("CT19","CT0") %in% CellTagName
Q4 <- c("CT5","CT10", "CT15") %in% CellTagName

#check that there are celltags
if(sum(Q1,Q2,Q3,Q4)<1) stop("No recognized CellTags in index file!")

#check for CT0 v CT19 
if(sum(Q3)>1) stop("CT0 and CT19 should not be used at the same time")
Q3isCT19 <- "CT19" %in% CellTagName

q1 <- c("CT1","CT2","CT3"
	,"CT6","CT7","CT8"
	,"CT11","CT12")
q1 <- q1 %in% CellTagName

q2 <- c("CT16","CT17") %in% CellTagName
q3 <- c("CT13") %in% CellTagName
q4 <- c("CT4","CT9") %in% CellTagName


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
	#hedge on good separation between [Q2 Q3] : [q2 q3 CT15]
	#that regardless of potential peak split within [q1 q4 CT5 CT10]
	#assume [Q2 Q3] : [q2 q3 CT15] split dominates over splits internal to [q1 q4 CT5 CT10]
	#exception is when CT5 CT6 CT7 CT8 CT9 are missing
	
	#there are three independent lumps
	peakN <- sum(TRUE  # represents Q2 Q3 which has been checked
			, sum(q2, q3, "CT15" %in% CellTagName)>0
			, sum(c(q2, q3, "CT5", "CT10") %in% CellTagName)>0
			)
	
	CT.5.6.7.8.9 <- isTRUE(sum(c("CT5", "CT6", "CT7", "CT8", "CT9") %in% CellTagName)>0)
	CT.1.2.3.4 <- isTRUE(sum(c("CT1", "CT2", "CT3", "CT4") %in% CellTagName)>0)
	CT.11.12.10 <- isTRUE(sum(c("CT11", "CT12", "CT10") %in% CellTagName)>0)
	
	if (!CT.5.6.7.8.9 && (CT.1.2.3.4 && CT.11.12.10)) peakN <- peakN + 1
	
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
			, sum(c("CT2", "CT7", "CT12", "CT17", "CT3", "CT8") %in% CellTagName)>0
			, sum(c("CT1", "CT6", "CT11", "CT16") %in% CellTagName)>0
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
	#calculate peakN, only options are 2, 3, or 4 for the exact number of CTs in Q2 Q3
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
#gate extreme wings first CT4 CT16
CT.4 <- "CT4" %in% CellTagName
CT.3.9 <- isTRUE(sum(c("CT3","CT9") %in% CellTagName) > 0)
CT.2.8 <- isTRUE(sum(c("CT2","CT8") %in% CellTagName) > 0)
CT.1.7.13 <- isTRUE(sum(c("CT1","CT7", "CT13") %in% CellTagName) > 0)
CT.6.12 <- isTRUE(sum(c("CT6","CT12") %in% CellTagName) > 0)
CT.11.17 <- isTRUE(sum(c("CT11","CT17") %in% CellTagName) > 0)
CT.16 <- "CT16" %in% CellTagName

#gate Q1-CT4
hi <- isTRUE(sum(CT.3.9, CT.2.8, CT.1.7.13, CT.6.12, CT.11.17, CT.16)>0)
peakN <- sum(CT.4, CT.3.9, CT.2.8, CT.1.7.13, CT.6.12, CT.11.17, CT.16)
if (CT.4 & hi) {
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "CT4p", "fRatio-", "Q1", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "Q1-CT4", "fRatio+", "Q1", "fRatio", "specificMindensity", args)
} else if (CT.4) { 
	df <- addGtRow(df, "CT4p", "fRatio+", "Q1", "fRatio", "dummySubgate") 
} else if (hi) {
	df <- addGtRow(df, "Q1-CT4", "fRatio+", "Q1", "fRatio", "dummySubgate")
}

#gate Q1-CT4-CT16
low <- isTRUE(sum(CT.3.9, CT.2.8, CT.1.7.13, CT.6.12, CT.11.17)>0)
peakN <- sum(CT.3.9, CT.2.8, CT.1.7.13, CT.6.12, CT.11.17, CT.16)
if (CT.16 & low) {
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, ", selectFromRight=TRUE", sep="")
	df <- addGtRow(df, "CT16p", "fRatio+", "Q1-CT4", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "Q1-CT4-CT16", "fRatio-", "Q1-CT4", "fRatio", "specificMindensity", args)
} else if (CT.16) { 
	df <- addGtRow(df, "CT16p", "fRatio+", "Q1-CT4", "fRatio", "dummySubgate") 
} else if (low) {
	df <- addGtRow(df, "Q1-CT4-CT16", "fRatio+", "Q1-CT4", "fRatio", "dummySubgate")
}





#G.n q ref split
alias <- "q1q2"
pop <- "G.n+"
parent <- "Q1-CT4-CT16"
dims <- "G.n"

CT.9 <- "CT9" %in% CellTagName
CT.17 <- "CT17" %in% CellTagName

if (sum(q1,CT.9)>0 & sum(CT.17,q3)>0) {
	#calculate peakN, only options are 2, 3 or 4

	#there are three independent lumps
	peakN <- sum(TRUE  # represents CT.17 q3 which has been checked
			, sum(c("CT6", "CT7", "CT8", "CT9") %in% CellTagName)>0
			, sum(c("CT11", "CT12") %in% CellTagName)>0
			, sum(c("CT1", "CT2", "CT3") %in% CellTagName)>0
			)

	df <- addGtRow(df, alias, pop, parent,  dims, "specificMindensity", paste(
			"selectedPeak=1, totalPeaks="
			,peakN
			,", eventThreshold="
			,eventThreshold
			, sep=""
			)
	)
} else if (sum(q1,CT.9)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate")
} else if (sum(CT.17,q3)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate", "minAtmin=F")
}


#R.n q ref split
alias <- "q1q4"
pop <- "R.n+"
parent <- "Q1-CT4-CT16"
dims <- "R.n"

if (sum(q1,CT.17)>0 & sum(q3,CT.9)>0) {
	#calculate peakN, only options are 2, 3 or 4

	#there are three independent lumps
	peakN <- sum(TRUE  # represents q3 CT.9 which has been checked
			, sum(c("CT2", "CT7", "CT12", "CT17") %in% CellTagName)>0
			, sum(c("CT3", "CT8") %in% CellTagName)>0
			, sum(c("CT1", "CT6", "CT11") %in% CellTagName)>0
			)

	df <- addGtRow(df, alias, pop, parent,  dims, "specificMindensity", paste(
			"selectedPeak=1, totalPeaks="
			,peakN
			,", eventThreshold="
			,eventThreshold
			, sep=""
			)
	)
} else if (sum(q1,CT.17)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate")
} else if (sum(q3,CT.9)>0) {
	df <- addGtRow(df, alias, pop, parent,  dims, "dummySubgate", "minAtmin=F")
}

#q gates
parent <- "Q1-CT4-CT16"
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
low <- "CT15" %in% CellTagName
hi <- isTRUE(sum(c("CT5","CT10") %in% CellTagName) > 0)
if (low & hi) {
	peakN <- sum(c("CT5","CT10", "CT15") %in% CellTagName)
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "CT15p", "G.n-", "Q4", "G.n", "specificMindensity", args)
	df <- addGtRow(df, "CT.5.10", "G.n+", "Q4", "G.n", "specificMindensity", args)
} else if (low) { 
	df <- addGtRow(df, "CT15p", "G.n+", "Q4", "G.n", "dummySubgate") 
} else if (hi) { 
	df <- addGtRow(df, "CT.5.10", "G.n+", "Q4", "G.n", "dummySubgate")
}

low <- "CT10" %in% CellTagName
hi <- "CT5" %in% CellTagName
if (low & hi) {
	df <- addGtRow(df, "CT10p", "G.n-", "CT.5.10", "G.n", "mindensity")
	df <- addGtRow(df, "CT5p", "G.n+", "CT.5.10", "G.n", "mindensity")
} else if (low) {
	df <- addGtRow(df, "CT10p", "G.n+", "CT.5.10", "G.n", "dummySubgate")
} else if (hi) { 
	df <- addGtRow(df, "CT5p", "G.n+", "CT.5.10", "G.n", "dummySubgate")
}


## gates in Q2 ##
low <- "CT22" %in% CellTagName
hi <- isTRUE(sum(c("CT20","CT21") %in% CellTagName) > 0)
if (low & hi) {
	peakN <- sum(c("CT20","CT21", "CT22") %in% CellTagName)
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "CT22p", "R.n-", "Q2", "R.n", "specificMindensity", args)
	df <- addGtRow(df, "CT.20.21", "R.n+", "Q2", "R.n", "specificMindensity", args)
} else if (low) {
	df <- addGtRow(df, "CT22p", "R.n+", "Q2", "R.n", "dummySubgate")
} else if (hi) {
	df <- addGtRow(df, "CT.20.21", "R.n+", "Q2", "R.n", "dummySubgate")
}

low <- "CT21" %in% CellTagName
hi <- "CT20" %in% CellTagName
if (low & hi) {
	df <- addGtRow(df, "CT21p", "R.n-", "CT.20.21", "R.n", "mindensity")
	df <- addGtRow(df, "CT20p", "R.n+", "CT.20.21", "R.n", "mindensity")
} else if (low) {
	df <- addGtRow(df, "CT21p", "R.n+", "CT.20.21", "R.n", "dummySubgate")
} else if (hi) {
	df <- addGtRow(df, "CT20p", "R.n+", "CT.20.21", "R.n", "dummySubgate")
}


### subgates in q4, q2 ###
# these are just q4 == CT9p, q2 == CT17p


### subgates in q1 ###
## gate q1 low wing ##
CT.3 <- "CT3" %in% CellTagName
CT.2.8 <- isTRUE(sum(c("CT2","CT8") %in% CellTagName) > 0)
CT.1.7 <- isTRUE(sum(c("CT1","CT7") %in% CellTagName) > 0)
CT.6.12 <- isTRUE(sum(c("CT6","CT12") %in% CellTagName) > 0)
CT.11 <- "CT11" %in% CellTagName
hi <- isTRUE(sum(CT.1.7, CT.6.12, CT.11)>0)

peakN <- sum(CT.3, CT.2.8, CT.1.7, CT.6.12, CT.11)

if ((CT.3 | CT.2.8) & hi) {
	if (CT.3 & CT.2.8) { 
		args <- paste("selectedPeak=2, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	} else {
		args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	}
	df <- addGtRow(df, "CT.2.3.8", "fRatio-", "q1", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "CT.1.6.7.11.12", "fRatio+", "q1", "fRatio", "specificMindensity", args)
} else if (CT.3 | CT.2.8) { 
	df <- addGtRow(df, "CT.2.3.8", "fRatio+", "q1", "fRatio", "dummySubgate") 
} else if (hi) { 
	df <- addGtRow(df, "CT.1.6.7.11.12", "fRatio+", "q1", "fRatio", "dummySubgate")
}


# subgates in q1 low wing #
if (CT.3 & CT.2.8) {
	args <- paste("selectedPeak=1, totalPeaks=2", ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "CT3p", "fRatio-", "CT.2.3.8", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "CT.2.8", "fRatio+", "CT.2.3.8", "fRatio", "specificMindensity", args)
} else if (CT.3) { 
	df <- addGtRow(df, "CT3p", "fRatio+", "CT.2.3.8", "fRatio", "dummySubgate") 
} else if (CT.2.8) { 
	df <- addGtRow(df, "CT.2.8", "fRatio+", "CT.2.3.8", "fRatio", "dummySubgate")
}

low <- "CT8" %in% CellTagName
hi <- "CT2" %in% CellTagName
if (low & hi) {
	df <- addGtRow(df, "CT8p", "fTotal-", "CT.2.8", "fTotal", "mindensity")
	df <- addGtRow(df, "CT2p", "fTotal+", "CT.2.8", "fTotal", "mindensity")
} else if (low) { 
	df <- addGtRow(df, "CT8p", "fTotal+", "CT.2.8", "fTotal", "dummySubgate")
} else if (hi) { 
	df <- addGtRow(df, "CT2p", "fTotal+", "CT.2.8", "fTotal", "dummySubgate")
}

## gate q1 high wing ##
hi <- isTRUE(sum(CT.6.12, CT.11) > 0)

if (CT.1.7 & hi) {
	peakN <- sum(CT.1.7, CT.6.12, CT.11)
	args <- paste("selectedPeak=1, totalPeaks=", peakN, ", eventThreshold=", eventThreshold, sep="")
	df <- addGtRow(df, "CT.1.7", "fRatio-", "CT.1.6.7.11.12", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "CT.6.11.12", "fRatio+", "CT.1.6.7.11.12", "fRatio", "specificMindensity", args)
} else if (CT.1.7) { 
	df <- addGtRow(df, "CT.1.7", "fRatio+", "CT.1.6.7.11.12", "fRatio", "dummySubgate") 
} else if (hi) { 
	df <- addGtRow(df, "CT.6.11.12", "fRatio+", "CT.1.6.7.11.12", "fRatio", "dummySubgate")
}

# subgates in q1 high wing #
if (CT.6.12 & CT.11) {
	args <- paste("selectedPeak=1, totalPeaks=2", ", eventThreshold=", eventThreshold, ", selectFromRight=TRUE", sep="")
	df <- addGtRow(df, "CT.6.12", "fRatio-", "CT.6.11.12", "fRatio", "specificMindensity", args)
	df <- addGtRow(df, "CT11p", "fRatio+", "CT.6.11.12", "fRatio", "specificMindensity", args)
} else if (CT.6.12) { 
	df <- addGtRow(df, "CT.6.12", "fRatio+", "CT.6.11.12", "fRatio", "dummySubgate")
} else if (CT.11) { 
	df <- addGtRow(df, "CT11p", "fRatio+", "CT.6.11.12", "fRatio", "dummySubgate")
}


low <- "CT12" %in% CellTagName
hi <- "CT6" %in% CellTagName
if (low & hi) {
	df <- addGtRow(df, "CT12p", "fTotal-", "CT.6.12", "fTotal", "mindensity")
	df <- addGtRow(df, "CT6p", "fTotal+", "CT.6.12", "fTotal", "mindensity")
} else if (low) { 
	df <- addGtRow(df, "CT12p", "fTotal+", "CT.6.12", "fTotal", "dummySubgate")
} else if (hi) { 
	df <- addGtRow(df, "CT6p", "fTotal+", "CT.6.12", "fTotal", "dummySubgate")
}

## subgates in q1 mid wing ##
low <- "CT7" %in% CellTagName
hi <- "CT1" %in% CellTagName
if (low & hi) {
	df <- addGtRow(df, "CT7p", "fTotal-", "CT.1.7", "fTotal", "mindensity")
	df <- addGtRow(df, "CT1p", "fTotal+", "CT.1.7", "fTotal", "mindensity")
} else if (low) { 
	df <- addGtRow(df, "CT7p", "fTotal+", "CT.1.7", "fTotal", "dummySubgate")
} else if (hi) { 
	df <- addGtRow(df, "CT1p", "fTotal+", "CT.1.7", "fTotal", "dummySubgate")
}



### all final population gates ###
dims <- "G.n,R.n"
method <- "flowClust"
args <- "K=1, quantile=0.8"

if ("CT1" %in% CellTagName) df <- addGtRow(df, "CT1", "CT1", "CT1p", dims, method, args)
if ("CT2" %in% CellTagName) df <- addGtRow(df, "CT2", "CT2", "CT2p", dims, method, args)
if ("CT3" %in% CellTagName) df <- addGtRow(df, "CT3", "CT3", "CT3p", dims, method, args)
if ("CT4" %in% CellTagName) df <- addGtRow(df, "CT4", "CT4", "CT4p", dims, method, args)
if ("CT5" %in% CellTagName) df <- addGtRow(df, "CT5", "CT5", "CT5p", dims, method, args)
if ("CT6" %in% CellTagName) df <- addGtRow(df, "CT6", "CT6", "CT6p", dims, method, args)
if ("CT7" %in% CellTagName) df <- addGtRow(df, "CT7", "CT7", "CT7p", dims, method, args)
if ("CT8" %in% CellTagName) df <- addGtRow(df, "CT8", "CT8", "CT8p", dims, method, args)
if ("CT9" %in% CellTagName) df <- addGtRow(df, "CT9", "CT9", "q4", dims, method, args)
if ("CT10" %in% CellTagName) df <- addGtRow(df, "CT10", "CT10", "CT10p", dims, method, args)
if ("CT11" %in% CellTagName) df <- addGtRow(df, "CT11", "CT11", "CT11p", dims, method, args)
if ("CT12" %in% CellTagName) df <- addGtRow(df, "CT12", "CT12", "CT12p", dims, method, args)
if ("CT13" %in% CellTagName) df <- addGtRow(df, "CT13", "CT13", "q3", dims, method, args)
if ("CT15" %in% CellTagName) df <- addGtRow(df, "CT15", "CT15", "CT15p", dims, method, args)
if ("CT16" %in% CellTagName) df <- addGtRow(df, "CT16", "CT16", "CT16p", dims, method, args)
if ("CT17" %in% CellTagName) df <- addGtRow(df, "CT17", "CT17", "q2", dims, method, args)
if ("CT19" %in% CellTagName) df <- addGtRow(df, "CT19", "CT19", "Q3", dims, method, args)
if ("CT0" %in% CellTagName) df <- addGtRow(df, "CT0", "CT0", "Q3", dims, method, args)
if ("CT20" %in% CellTagName) df <- addGtRow(df, "CT20", "CT20", "CT20p", dims, method, args)
if ("CT21" %in% CellTagName) df <- addGtRow(df, "CT21", "CT21", "CT21p", dims, method, args)
if ("CT22" %in% CellTagName) df <- addGtRow(df, "CT22", "CT22", "CT22p", dims, method, args)

res <- list("gt"=df
	, "G.nRange"=c(as.numeric(G.nMIN), as.numeric(G.nMAX))
	, "R.nRange"=c(as.numeric(R.nMIN), as.numeric(R.nMAX))
)
return(res)

}

