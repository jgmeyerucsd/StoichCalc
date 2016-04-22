#############################################
####### checks if both heavy and light precursors are in at least one window together, otherwise sets window index ==1 and therefore values output are NA
####################################

ms2xic=function(file="151204_0008_SStoich_lact_bRP8.mzXML",
	wind.file="windoe.txt",
	windowtable="windoe.txt",
	precMz=674.3494,
	precZ = 3,
	fragMz=945.4499, 
	ppm=40,
	rtrange=c(45,48),
	type="h",
	Kcount=2,
	moddelta=4.025107){
	
  ###	requires mzR package
	require(mzR)
  
	###	read in window definitions from txt file
	windows<-read.delim(file=wind.file,head=F,stringsAsFactors=F)
	nwindows<-nrow(windows)
	
	###	open MS file and get headers
	ms <- openMSfile(filename=file)
	hd<-header(ms)
	nscans<-length(hd[[1]])

	### alternately, add a element to hd containing the swath window index
	ncycles<-length(hd[[1]])/(nwindows+1)
	hd[[22]]<-rep(seq(from=0,to=64),times=ncycles)
	
	###	which window(s) is precMz in?
	precWinIndex<-which(apply(as.matrix(windows),1,findInterval,x=precMz)==1)

	## a set of spectra of interest: MS2 scans containing the precursor m/z between retention time
	#### part to check if those in single window are in the same window h/l
	if(length(precWinIndex)==1){
	  ms2 <- which(hd$V22 == precWinIndex)
	  if(type=="h"){
	    precMzl<-precMz-(moddelta*Kcount)/precZ
	    precMzh<- precMz
  	}
  	if(type=="l"){
  	  precMzl<-precMz
  	  precMzh<- precMz+(moddelta*Kcount)/precZ
  	}
	  precWinIndexl<-which(apply(as.matrix(windows),1,findInterval,x=precMzl)==1)
	  precWinIndexh<-which(apply(as.matrix(windows),1,findInterval,x=precMzh)==1)

  	bothin<-intersect(precWinIndexh,precWinIndexl)
  	if(length(bothin)==1){
  	  precWinIndex<-bothin
  	}
  	if(length(bothin)!=1){
  	  precWinIndex<-c(1)
  	  print("heavy and light precursor in different SWATH windows")
  	  print("failed on single window")
  	  warning("heavy and light precursor in different SWATH windows")
  	  print("this should print if warn")
  	}
	}
	
	#### if precursor is in two windows, decide which to use based on type (h/l) and modmass
	if(length(precWinIndex)==2){
	  tempwin<-list(windows[precWinIndex[1],],windows[precWinIndex[2],])
	  if(type=="h"){
	    precMzl<-precMz-moddelta/precZ
	    precMzh<- precMz
	  }
	  if(type=="l"){
	    precMzl<-precMz
	    precMzh<- precMz+moddelta/precZ
	  }
	  inboth<-c()
	  for(i in 1:2){
	    inboth<-c(inboth,findInterval(precMzl,tempwin[[i]])==1 & findInterval(precMzh,tempwin[[i]])==1)
	  }
	  
    if(length(which(inboth=="TRUE"))>0){
      precWinIndex<-precWinIndex[which(inboth=="TRUE")]
    }
	  if(length(which(inboth=="TRUE"))==0){
	    precWinIndex<-c(1)
	    print("heavy and light precursor in different SWATH windows")
	    print("failed on double window")
	    warning("heavy and light precursor in different SWATH windows")
	    print("this should print if warned")
	  }
	}
	#if(length(precWinIndex)==1){
	  ms2 <- which(hd$V22 == precWinIndex[1])
	  rtsel <- hd$retentionTime[ms2] / 60 > rtrange[1] & hd$retentionTime[ms2] / 60 < rtrange[2]
	  rtsel[which(rtsel==T)]
	  ### list of scans within range that are ms2 scans
	  scanpeaks<-peaks(ms,scan=ms2[rtsel])
	  scanhead<-hd[ms2[rtsel],]
	  #### time values in minutes
	  rtvalues<-header(ms,scan=ms2[rtsel])[,"retentionTime"]/60
	  ### empty list of int values and vector of intensities
	  tempsum<-list()
	  
	  intvec<-rep(0,times=length(scanpeaks))
	  #### m/z range to extract
	  fragRange<- c(fragMz-(fragMz*ppm*1e-6),fragMz+(fragMz*ppm*1e-6))
	  #### loop to extract mz range from rt range
	  for(x in 1:length(scanpeaks)){
	    tempscans<-scanpeaks[[x]][,1]>fragRange[1] & scanpeaks[[x]][,1]< fragRange[2]
	    tempsum[[x]]<-scanpeaks[[x]][tempscans]
	    if(length(tempsum[[x]])>2){
	      #sum(tempsum[[x]] [((length(tempsum[[x]])/2)+1):length(tempsum[[x]])])
	      intvec[x]<-sum(tempsum[[x]] [((length(tempsum[[x]])/2)+1):length(tempsum[[x]])])
	    }
	    if(length(tempsum[[x]])<3){
	      intvec[x]<-tempsum[[x]][2]
	    }
	  }
	  print(fragMz)
	  return(data.frame(rtvalues,intvec))
	
	if(length(precWinIndex)==1){
	  return(data.frame())
	}
}