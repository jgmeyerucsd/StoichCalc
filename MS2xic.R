
### read in txt file with windows, chose
#txtlist<-list.files(pattern="txt")

#windows<-read.delim(file="windoe.txt",head=F,row.names=F)
#y8<-ms2xic(fragMz=1002.4713,ppm=50)
#y7<-ms2xic(fragMz=945.4499,ppm=50)
#y6<-ms2xic(fragMz=782.3865,ppm=50)
#y2<-ms2xic(fragMz=295.12988,ppm=50)
#b2<-ms2xic(fragMz=346.2275,ppm=50)
#b9<-ms2xic(fragMz=600.8228,ppm=50)


#test1<-ms2xic(file=mzxml,rtrange=peakbound,precMz=643.3304,fragMz = 968.4948)

ms2xic=function(file="151204_0008_SStoich_lact_bRP8.mzXML",
	wind.file="windoe.txt",
	windowtable="windoe.txt",
	precMz=674.3494,
	precZ = 3,
	fragMz=945.4499, 
	ppm=40,
	rtrange=c(45,48),
	type="h",
	moddelta=4.025107){

	###	requires mzR package
	require(mzR)

	###	read in window definitions from txt file
	windows<-read.delim(file=wind.file,head=F,stringsAsFactors=F)
	nwindows<-nrow(windows)

	###	open MS file and get headers
	ms <- openMSfile(filename=file)
	hd<-header(ms)

	###	make index of which ms2 scans are which window
	nscans<-length(hd[[1]])
	windowindexlist<-list()
	for(i in 1:nwindows){
		#print(i)
		windowindexlist[[i]]<-seq(from=i+1,to=nscans,by=65)
		}
	### alternately, add a element to hd containing the swath window index
	ncycles<-length(hd[[1]])/(nwindows+1)
	hd[[22]]<-rep(seq(from=0,to=64),times=ncycles)
	
	###	which window(s) is precMz in?
	precWinIndex<-c()
	for(i in 1:nwindows){
		if(findInterval(precMz,windows[i,])==1){
		precWinIndex<-c(precWinIndex,i)
			}
		}
	## a set of spectra of interest: MS2 scans containing the precursor m/z
	## between retention time
	if(length(precWinIndex)==1){
	  ms2 <- which(hd$V22 == precWinIndex)
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
	    print("heavy and light precursor in different SWATH windows")
	    stop("heavy and light precursor in different SWATH windows")
	    print("this should not print")
	  }
	  

	    
	}
	  
	ms2 <- which(hd$V22 == precWinIndex[1] | hd$V22 == precWinIndex[2])
	rtsel <- hd$retentionTime[ms2] / 60 > rtrange[1] & hd$retentionTime[ms2] / 60 < rtrange[2]
	rtsel[which(rtsel==T)]

	
	### list of scans within range that are ms2 scans
	scanpeaks<-peaks(ms,scan=ms2[rtsel])
	scanhead<-hd[ms2[rtsel],]

	#### time values in minutes
	rtvalues<-header(ms,scan=ms2[rtsel])[,"retentionTime"]/60
	
	### empty list of int values and vector of intensities
	tempsum<-list()
	intvec<-c()
	#### m/z range to extract
	fragRange<- c(fragMz-(fragMz*ppm*1e-6),fragMz+(fragMz*ppm*1e-6))
	#### loop to extract mz range from rt range
	for(x in 1:length(scanpeaks)){
		tempscans<-scanpeaks[[x]][,1]>fragRange[1] & scanpeaks[[x]][,1]< fragRange[2]
		tempsum[[x]]<-scanpeaks[[x]][tempscans]
		if(length(tempsum[[x]])>2){
			sum(tempsum[[x]] [((length(tempsum[[x]])/2)+1):length(tempsum[[x]])])
			intvec<-c(intvec,sum(tempsum[[x]] [((length(tempsum[[x]])/2)+1):length(tempsum[[x]])]))
			}
		if(length(tempsum[[x]])<3){
			intvec<-c(intvec,tempsum[[x]][2])
			}
		#print(scanpeaks[[x]][tempscans])
		}
	print(fragMz)
	return(data.frame(rtvalues,intvec))
	
	}



