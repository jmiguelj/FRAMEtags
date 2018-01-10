
#####initialize
#library(flowWorkspace)
#library(openCyto)
#library(ggcyto)
#library(parallel)
#library(RColorBrewer)
#library(gtools)
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
#' @import flowClust
#' @export


getFRAMEtagData <- function() { # added to make script into function

#detect system variables
numCores <- detectCores()-1
#scriptDir <- dirname(sys.frame(1)$ofile)

# adjust param for plotGate()
flowWorkspace.par.set("plotGate", list(default.y = "S.A")) # set default Y-axis from "SSC-A" to "S.A"

#adjust minimum events to avoid flowClust pipeline error on incorrect gating
options("openCyto" = list(gating = list(minEvents = 5)))

#####retrieve FCS files
# ask for input folder
inFolder <- readline(prompt="Enter path of folder FCS files: ")

# create folder to output files
dir.create(file.path(inFolder, "output"))

##### MAKE GATING TEMPLATE
#detect/ask for FRAME-tag index file
ftFile <- list.files(path = inFolder, pattern = "Tag-Strain", full = TRUE, ignore.case = TRUE)
if (length(ftFile)==0) ftFile <- readline(prompt="Enter path of Tag-Strain.csv index file: ")
FTindex = read.csv(ftFile, stringsAsFactors = FALSE)





GTfile <- list.files(path = inFolder, pattern = "gating_template", full = TRUE, ignore.case = TRUE) # need to extract fluor limits into template var for prettyPlotGate
if (length(GTfile)==0) {
	cat("\nCreating gating template...")
	#load gating template function, takes vector of FT names and returns gating template as a data frame
	#source(file.path(scriptDir,"makeGT_v4.R"), local=TRUE)
	
	#extract FT names and corresponding strain identities
	TagName <- FTindex[["Tag"]]

	template <- makeGT(TagName) #create gating template dataframe
	df <- template$gt

	# make cvs from df
	GTfile = file.path(inFolder, "gating_template.csv")
	write.csv(df,GTfile,row.names=FALSE)

	cat("Done\n\n")
}


#####detect/ask for tube identities file
tubeFile <- list.files(path = inFolder, pattern = "Tube-Treatment", full = TRUE, ignore.case = TRUE)
if (length(tubeFile)==0) tubeFile <- readline(prompt="Enter path of Tube-Treatment.csv index file: ")
Tubeindex = read.csv(tubeFile, stringsAsFactors = FALSE)
if (tolower(colnames(Tubeindex)[1]) == "well") { #tolower input to ignore case
	isWell <- TRUE
	colnames(Tubeindex)[1] <- "Tube" # Relabel this precise case-sensitive label used in countFTs
} else {
	isWell <- FALSE
	colnames(Tubeindex)[1] <- "Tube" # Relabel this precise case-sensitive label used in countFTs
}

##### Register Custom Gating Functions and transforms
#register custom functions: minmaxTransform, specificMindensity, dummySubgate
#source(file.path(scriptDir,"customFunctions_v4.R"), local=TRUE)


#### detect/create analyzed files log and retrieve only new fcs files


logName = "files-analyzed_log.csv"
analyzedFile <- list.files(path = inFolder, pattern = logName, full = TRUE, ignore.case = TRUE)
if (length(analyzedFile) != 0) {
	analyzedIndex = read.csv(analyzedFile, stringsAsFactors = FALSE)
} else {
	setSize <- readline(prompt="Number of files in processing batch? (default=12):")
	if (setSize == "") setSize <- 12
	setSize <- as.numeric(setSize)
	
	analyzedIndex <- data.frame(entry=numeric()
				, files_analyzed=character()
				, set_size=numeric()
				, batch=numeric()
				,stringsAsFactors = FALSE)
	row <- list(entry=0, files_analyzed="current progress->", set_size=setSize, batch=0)
	analyzedIndex[nrow(analyzedIndex)+1, names(row)] <- row
}

setSize <- analyzedIndex[analyzedIndex[,"entry"]==0,"set_size"]
prevBatch <- analyzedIndex[analyzedIndex[,"entry"]==0,"batch"]

# get names of csv files to analyze, path defaulting to current working directory
fcsFileNames <- list.files(path = inFolder, pattern = "*.fcs", full = FALSE)
fcsFilePaths <- list.files(path = inFolder, pattern = "*.fcs", full = TRUE)

# only keep paths for new FCS files not previously analyzed (in log)
fcsFilePaths <- fcsFilePaths[fcsFileNames %in% setdiff(fcsFileNames, analyzedIndex[,"files_analyzed"])]

if (length(fcsFilePaths)==0) stop("All files in folder already analyzed")

##### Detect channel names and ask which to use, 
# load FCS files
cat("\nLoading fcs files (may take a minute)......")
fs  <- read.flowSet(fcsFilePaths, alter.names=TRUE)
cat("Done\n\n")

#extract parameter names and display
channelNames <- colnames(fs) # s4 object slot
cat("The detected parameters are:", channelNames,sep="\n")

# ask for correct channels
gfpName <- readline(prompt="Enter name of GFP channel: ")
rfpName <- readline(prompt="Enter name of RFP channel: ")
#bfpName <- readline(prompt="Enter name of BFP channel: ")
sscName <- readline(prompt="Enter name of SSC Area channel: ")
fscName <- readline(prompt="Enter name of FSC Area channel: ")
#fschName <- readline(prompt="Enter name of FSC Height channel: ")


if(gfpName == "") gfpName <- "FITC.A"
if(rfpName == "") rfpName <- "PE.Texas.Red.A"
#if(bfpName == "") bfpName <- "Cascade.Blue.A"
if(sscName == "") sscName <- "SSC.A"
if(fscName == "") fscName <- "FSC.A"


#set parameter names to match internal variable names
cat("Updating channel names...")
gfpIndex <- which(colnames(fs)==gfpName)
rfpIndex <- which(colnames(fs)==rfpName)
#bfpIndex <- which(colnames(fs)==bfpName)
sscIndex <- which(colnames(fs)==sscName)
fscIndex <- which(colnames(fs)==fscName)
#fschIndex <- which(colnames(fs)==fschName)


# rename columns for custom transform
colnames(fs)[gfpIndex] <- "G"
colnames(fs)[rfpIndex] <- "R"
#colnames(fs)[bfpIndex] <- "B"
colnames(fs)[sscIndex] <- "S.A"
colnames(fs)[fscIndex] <- "F.A" # can't use just "F", reserved for false
#colnames(fs)[fschIndex] <- "F.H"
cat("Done\n")

#fcsFilePaths[((i*12)+1):((i*12)+12)]
fsALL <- fs

######LOOP in chunks of FCS files


for (i in 0:((length(fsALL)/setSize))) {

	start <- (i*setSize)+1
	end <- (i*setSize)+setSize
	if (start > length(fsALL))break
	if (end > length(fsALL)) end <- (i*setSize) + (length(fsALL)%%setSize)

	fs <- fsALL[start:end]
	#assign("fs", fsALL[start:end], envir = .GlobalEnv)
	
	##### TRANSFORM
	#transform data to add calculated parameters: G.n, R.n, B.n, fTotal and fRatio
	cat("\nStarting analysis...\n")

	# load transform function, requires flowset with colnames: G.n, R.n, S.A, F.A
	#source(file.path(scriptDir,"transformGR.R"), local=FALSE)
	#fs <- transform(fs, G.n = scaleBy*(G/S.A), R.n = scaleBy*(R/S.A))
	#nanTrunc <- nanTransform(transformationId="nanTrunc", a=-1000)

	
	fs <- transformGR(fs)


	##### GATE FT POPULATIONS
	gt <- gatingTemplate(GTfile, autostart = 1L)
	gs <- GatingSet(fs)
	#err <- try(gating(gt, gs, parallel_type="multicore",mc.cores=numCores), silent = T)
	err <- gating(gt, gs, parallel_type="multicore",mc.cores=numCores)
	
	#if (class(err)=="try-error") {
		#break
		#errFiles <- unlist(strsplit(err[1], "failed at "))
		#errFiles <- errFiles[even(1:length(errFiles))]
		#sapply(errFiles, strsplit, split="\n")
	#}

	#print gating result
	cat("\nSaving summary plots...\n")
	#source(file.path(scriptDir,"prettyPlotGate.R"), local=TRUE)

	p <- prettyPlotGate(gs, FTindex$Tag)

	outGates <- paste((prevBatch+i+1), "_FT-gates.png", sep="")

	png(filename=file.path(inFolder, "output", outGates)
		, width = 2048, height = 2048)
	plot(p)

	dev.off()

	##### OUTPUT DATA FILES
	cat("\nSaving output files...\n")
	#source(file.path(scriptDir,"countFTs_v2.R"), local=TRUE)

	res <- countFTs(gs, FTindex, Tubeindex, well=isWell)

	outStats <- paste((prevBatch+i+1), "_full-stats.csv", sep="")
	outCounts <- paste((prevBatch+i+1), "_counts.csv", sep="")
	outPercents <- paste((prevBatch+i+1), "_percent.csv", sep="")
	#outRate <- paste((prevBatch+i+1), "_flowrate.csv", sep="")
	outGT <- paste((prevBatch+i+1), "_gating_template.csv", sep="")

	write.csv(res$fullStats,file=file.path(inFolder,"output", outStats),row.names=TRUE)
	write.csv(res$counts,file=file.path(inFolder,"output", outCounts),row.names=TRUE)
	write.csv(res$percents,file=file.path(inFolder,"output", outPercents),row.names=TRUE)
	#write.csv(res$flowrate,file=file.path(inFolder,"output", outRate),row.names=FALSE)
	file.copy(GTfile, file.path(inFolder,"output", outGT), overwrite=T)
	
	for (j in 1:length(colnames(res$fullStats))) {
		row <- list(entry=length(analyzedIndex$entry)
				, files_analyzed=colnames(res$fullStats)[j]
				, set_size=length(fs)
				, batch=prevBatch+i+1)
		analyzedIndex[nrow(analyzedIndex)+1, names(row)] <- row
	}
	analyzedIndex[analyzedIndex[,"entry"]==0,"batch"] <- prevBatch+i+1
	write.csv(analyzedIndex,file=file.path(inFolder, logName),row.names=FALSE)
	
}

####detailed plotting of detected peaks for specificMindensity
#ask if needed
whichToPlot <- readline(prompt="Do you want detailed plots of peaks/gates? Which samples? (#,#,# no spaces; blank = no plotting): ")

#only plot if input given
if (whichToPlot != "") {
	cat("Plotting detailed peaks/gates (will take several minutes) ...\n")
	#load functions
	#source(file.path(scriptDir,"plot_specificMindensity_peaks.R"), local=TRUE)
	
	#load csv gating template as data frame
	GTindex = read.csv(GTfile, stringsAsFactors = FALSE)
	
	#convert input to index list
	whichToPlot <- as.numeric(unlist(strsplit(whichToPlot, ",")))
	
	for (k in whichToPlot) {
		plotSampleName <- Tubeindex[Tubeindex[,"Tube"]==res$pData$Tube[k],"Treatment"]
		
		outPlotBaseName <- paste(k, "_", plotSampleName, "_plot_GT-line-", sep="")

		cat("\nSample ", k, "_", plotSampleName, "\ncalculating plots for...\n", sep="")

		plots <- specificMindensityPlots(gs, GTindex, whichFlowFrame=k, eventsToPlot=10000)

		cat("Saving plots...\n")
		for (j in 1:length(plots)) {
			
			outPlotName <- paste(outPlotBaseName, names(plots)[j], ".png", sep="")
			png(filename=file.path(inFolder, "output", outPlotName)
				, width = 1000, height = 800)
			grid.arrange(plots[[j]]$t, plots[[j]]$p, plots[[j]]$p2, plots[[j]]$p3, nrow=2, ncol=2)

			dev.off()
			cat("...", outPlotName, "\n", sep="")
		}
	}
	cat("Done!\n")

}


cat("############################## Finished!\n")

}
