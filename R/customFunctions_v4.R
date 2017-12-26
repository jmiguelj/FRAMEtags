
#####CUSTOM TRANSFORMS
#truncate for extreme negs and inf produced by SSC normalization/scaling to make logicle transform possible
minmaxTransform <- function(transformationId="minmaxTransform",a=1, b=262143)
{
    t <- new("transform", .Data=function(x){
        x[x<=a] <- a
        x[x>=b] <- b
        x
    })
    t@transformationId <- transformationId
    t
}

#remove NaNs produced by division to make logicle transform possible
nanTransform <- function(transformationId="nanTransform",a=-1000)
{
    t <- new("transform", .Data=function(x){
        x[is.nan(x)] <- a
        x
    })
    t@transformationId <- transformationId
    t
}




#####CUSTOM GATING METHODS
### find_valleysFull
## expanded valley finding function to include density value information in result
##

find_valleysFull <- function(x, y = NULL, num_valleys = NULL, adjust = 2, ...) {

  x <- as.vector(x)

  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find valleys.")
    return(NA)
  }
  
  if (is.null(y)) {
    dens <- density(x, adjust = adjust, ...)
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  # Discrete analogue to a second derivative applied to the KDE. See details.
  second_deriv <- diff(sign(diff(dens$y)))
  which_minima <- which(second_deriv == 2) + 1

  # The 'density' function can consider observations outside the observed range.
  # In rare cases, this can actually yield valleys outside this range. We remove
  # any such valleys.
  which_minima <- which_minima[findInterval(dens$x[which_minima], range(x)) == 1]

  # Next, we sort the valleys in descending order based on the density heights.
  which_minima <- which_minima[order(dens$y[which_minima], decreasing = FALSE)]

  # Returns the local minima. If there are none, we return 'NA' instead.
  if (length(which_minima) > 0) {
    valleys <- dens$x[which_minima]
    if (is.null(num_valleys) || num_valleys > length(valleys)) {
      num_valleys <- length(valleys)
    }
    valleys <- valleys[seq_len(num_valleys)]
  } else {
    valleys <- NA
  }
  valleys <- data.frame(x = valleys, y = dens$y[which_minima][seq_len(num_valleys)])
  valleysFull <- list("valleys" = valleys, "dens" = dens)
  valleysFull
}

#mindensity method were peaks to use are selected by left to right rank order
.specificMindensity <- function(fr, pp_res, channels=NA, filterId="specificMindensity", selectedPeak=1, totalPeaks=NULL, selectFromRight=FALSE, findAdjust=0.5, minAdjust=findAdjust, eventThreshold=20, separationThreshold=0.3, detailedRes=FALSE, ...){
	
	
	x <- exprs(fr)[, channels]

	#calculate full peak and valley information
	valleysFull <- find_valleysFull(x, adjust=findAdjust, ...)
	peaks <- openCyto:::.find_peaks(x=valleysFull$dens$x, y=valleysFull$dens$y ,adjust=findAdjust, ...) 
	valleys <- valleysFull$valleys
	
	#calculate density information
	dens <- valleysFull$dens
	totalEvents <- length(x)
	
	#plot(peaks)
	#points(valleys, col="green")
	
	#order by peaks and valleys by x
	peaks <- peaks[order(peaks[,"x"]),]
	valleys <- valleys[order(valleys[,"x"]),]

	####assure all valleys are bracketed by peaks by labeling, merging, ordering, and using diff() to find offending pairs
	peaks$type <- 1
	valleys$type <- 0
		
	peaks <- rbind(peaks, valleys)
	peaks <- peaks[order(peaks[,"x"]),]

	while(length(which(diff(peaks$type)==0)) > 0) {
		index <- min(which(diff(peaks$type)==0))
		if (peaks$type[index] == 1) { #its a peak
			if (peaks$y[index] >  peaks$y[index+1]) index <- index+1 #make index point to lower peak in offending pair
		} else { #its a valley
			if (peaks$y[index] <  peaks$y[index+1]) index <- index+1 #make index point to higher valley in offending pair
		}
		peaks <- peaks[-c(index),] #remove offending maxima
	}
	
	valleys <- peaks[peaks$type == 0,]
	peaks <- peaks[peaks$type == 1,]
	
	#remove extra extreme valleys
	if (peaks[1,"x"] > valleys[1,"x"]) valleys <- valleys[2:nrow(valleys),]
	if (peaks[nrow(peaks),"x"] < valleys[nrow(valleys),"x"]) valleys <- valleys[1:(nrow(valleys) - 1),]
	#now starts and ends with peaks
		

	#include valley information in peaks  ### ASSUMES valleys are always bracketed by peaks, aka nrow(peaks) = nrow(valleys) + 1
	peaks["valley1X"] <- c(dens$x[1],valleys[,"x"])
	peaks["valley1Y"] <- c(dens$y[1],valleys[,"y"])
	peaks["valley2X"] <- c(valleys[,"x"],dens$x[length(dens$x)])
	peaks["valley2Y"] <- c(valleys[,"y"],dens$y[length(dens$y)])

	#calculate events per peak
	peaks["events"] <- mapply(function(V1,V2){length(x[x > V1 & x < V2])}, peaks[,"valley1X"],peaks[,"valley2X"])

	peaks["ratio1"] <- 1 - (peaks["valley1Y"] / peaks["y"])
	peaks["ratio2"] <- 1 - (peaks["valley2Y"] / peaks["y"])

	#keep copy of original preprocessed peak info before merging
	unmergedPeaks <- peaks[,]
		
	###merge peaks, only replace ratios, replace valley info (only params used for final filtering)
	
	
	repeat {
		#select peak to merge based on events then ratio
		minRatio <- min(c(peaks$ratio1, peaks$ratio2))
		minEvents <- min(peaks$events)
		if (minEvents <= eventThreshold) {
			peakToMerge <- peaks[peaks$events == minEvents,]
			
			# at low event numbers crashes are expected so chose based on worst ratio
			minRatio <- min(c(peakToMerge$ratio1, peakToMerge$ratio2))
			peakToMerge <- peakToMerge[peakToMerge$ratio1 == minRatio | peakToMerge$ratio2 == minRatio,]

		} else {
			peakToMerge <- peaks[peaks$ratio1 == minRatio | peaks$ratio2 == minRatio,]
		}

		#in rare event ratios are identical, choose first one
		if (nrow(peakToMerge) != 1) peakToMerge <- peakToMerge[1,]

		
		#repeat until above separation threshold then break on expected peak number 
		if (peakToMerge$ratio1 > separationThreshold && peakToMerge$ratio2 > separationThreshold) {
			if (!is.null(totalPeaks)) {	
				if (nrow(peaks) <= totalPeaks) break
			} else {
				break #break at separation threshold if no peak number given
			}
		}

		#merge valley information left or right depending on which ratio is smaller
		if (peakToMerge$ratio1 < peakToMerge$ratio2) { #merge left
			
			 
			parentPeak <- peaks[which(peaks$x == peakToMerge$x) - 1,]
			if (which(peaks$x == peakToMerge$x) == 1) parentPeak <- peakToMerge #delete if already far left, merged onto itself and then deleted as usual (by x)
			parentPeak$valley2X <- peakToMerge$valley2X
			parentPeak$valley2Y <- peakToMerge$valley2Y
			parentPeak$ratio2 <- 1 - (parentPeak$valley2Y / parentPeak$y)
	
		} else { #merge right
			parentPeak <- peaks[which(peaks$x == peakToMerge$x) + 1,]
			if (which(peaks$x == peakToMerge$x) == nrow(peaks)) parentPeak <- peakToMerge #delete if already far right
			parentPeak$valley1X <- peakToMerge$valley1X
			parentPeak$valley1Y <- peakToMerge$valley1Y
			parentPeak$ratio1 <- 1 - (parentPeak$valley1Y / parentPeak$y)

		}
		
		#add events and copy this new peak back to peaks
		parentPeak$events <- parentPeak$events + peakToMerge$events
		peaks[peaks$x == parentPeak$x,] <- parentPeak

		#remove peakToMerge from peaks
		peaks <- peaks[peaks$x != peakToMerge$x,]

		
	}
	

	
	#order depending on direction of selection and chose selected "peak" by selecting valley
	peaks <- peaks[order(peaks[,"x"], decreasing=selectFromRight),]
	

	#use peak to select valley as this info was carried through merging (this really is where the gat needs to be)
	if (selectFromRight) targetValley <- peaks[selectedPeak, "valley1X"]
	if (!selectFromRight) targetValley <- peaks[selectedPeak, "valley2X"]	
		
	#use valley to select original bracketing peaks to make final gate, 
	#include closest match since some intervening peaks may have been removed in prefiltering
	#warning suppressed as returning -Inf or Inf for max() and min() is desireable
	min <- suppressWarnings(max(unmergedPeaks[unmergedPeaks[,"valley2X"] <= targetValley, "x"])) 
	max <- suppressWarnings(min(unmergedPeaks[unmergedPeaks[,"valley1X"] >= targetValley, "x"]))
	
	if (is.na(min)) min <- -Inf #this happens when targetValley is NA from selected peak being beyond index
	if (is.na(max)) max <- Inf
	
    my_gate <- openCyto:::.mindensity(fr, channel=channels, filterId = filterId, gate_range=c(min,max), adjust=minAdjust, ...)

    if (detailedRes) {
	args <- list("selectedPeak"=selectedPeak, "totalPeaks"=totalPeaks, "eventThreshold"=eventThreshold, "selectFromRight"=selectFromRight, "adjust" = findAdjust)
	res <- list("gate" = my_gate, "peaks" = peaks, "args" = args)
	return(res)
    } else {
	return(my_gate)
    }
}
registerPlugins(fun=.specificMindensity,methodName='specificMindensity',dep=NA)


#register gating method for when a CellTag is missing and want to gate whole range
.dummySubgate <- function(fr, pp_res, channels=NA, filterId="dummySubgate",minAtmin=T, ...){
    
	x <- exprs(fr)[, channels]
	if (minAtmin) {
		min <- min(x, na.rm=T) # ignore NaN as they can be produced by normalization
		max <- NULL
	} else {
		min <- max(x, na.rm=T)
		max <- NULL
	}
	my_gate <- openCyto:::.boundary(fr, channels=channels, min=min, max=max,...)
    	return(my_gate)
}
registerPlugins(fun=.dummySubgate,methodName='dummySubgate',dep=NA)


