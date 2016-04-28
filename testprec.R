


testprec=function(
  windows=read.delim(file="windoe.txt",head=F,stringsAsFactors=F),
  precMz=1003.8,
  precZ = 3,
  type="l",
  moddelta=4.025107
){
  nwindows<-nrow(windows)
  
  ###	open MS file and get headers
  
  #hd<-header(ms)
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
  if(precWinIndex==1){
    return(c())
  }
  if(precWinIndex!=1){
    return(c(precZ))
  }
  
}

