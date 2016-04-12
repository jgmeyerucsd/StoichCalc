
ms1xic=function(
	file="151204_0008_SStoich_lact_bRP8.mzXML",
	precMz=674.3494,
	ppm=20,
	rtrange=c(45,48))
	{

	require(mzR)
	ms <- openMSfile(filename=file)
	hd<-header(ms)
	## a set of spectra of interest: MS1 spectra eluted
	## between 30 and 35 minutes retention time
	ms1 <- which(hd$msLevel == 1)
	rtsel <- hd$retentionTime[ms1] / 60 > rtrange[1] & hd$retentionTime[ms1] / 60 < rtrange[2]
	### list of scans within range that are ms1 scans
	scanpeaks<-peaks(ms,scan=ms1[rtsel])
	#### time values in minutes
	rtvalues<-header(ms,scan=ms1[rtsel])[,"retentionTime"]/60
	### empty list of int values
	tempsum<-list()
	#intvec<-rep(0, times=length(scanpeaks)
	intvec<-c()
	#### m/z range to extract
	precRange<- c(precMz-(precMz*ppm*1e-6),precMz+(precMz*ppm*1e-6))
	#### loop to extract mz range from rt range
	for(x in 1:length(scanpeaks)){
		tempscans<-scanpeaks[[x]][,1]>precRange[1] & scanpeaks[[x]][,1]< precRange[2]
		tempsum[[x]]<-scanpeaks[[x]][tempscans]
		intvec<-c(intvec,tempsum[[x]][2])
		#print(scanpeaks[[x]][tempscans])
		}
	return(data.frame(rtvalues,intvec))
	}



m1<-ms1xic()
m2<-ms1xic(precMz=674.8508)
m3<-ms1xic(precMz=675.3515)



plot(m1,ty="l",col="blue",lwd=2)
lines(m2,ty="l",col="purple",lwd=2)
lines(m3,ty="l",col="brown",lwd=2)
#lines(rtvalues,intvec,col="black",lwd=2)
lines(y8,col="orange",lwd=2)
lines(y7,col="grey",lwd=2)
lines(y6,col="pink",lwd=2)
lines(y2,col=60,lwd=2)
lines(b2,col=69,lwd=2)
lines(b9,col="black",lwd=2)

tempsum
intvec
plot(rtvalues,intvec,ty="l")
peakscan25<-scanpeaks[[25]][,1]>precRange[1] & scanpeaks[[25]][,1]< precRange[2]


scanpeaks[[25]][peakscan25]

length(scanpeaks)
typeof(scanpeaks[1])
length(scanpeaks[[2]])

scanpeaks[,1]
peaks(ms,scan=ms1[rtsel][1])[,1]>674 & peaks(ms,scan=ms1[rtsel][1])[,1]<675



plot(M, aspect = 1, allTicks = FALSE)
plot3D(M)
?MSmap


