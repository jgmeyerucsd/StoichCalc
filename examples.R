
stoich1<-stoichwrapper(ppm=40,mzrange=c(100,1500))
stoich3<-stoichwrapper(ppm=20,mzrange=c(100,1500))
stoich4<-stoichwrapper(ppm=20,mzrange=c(100,1500),threshold=100)

stoich0<-stoichwrapper(mzxml="G:/tmp/BSAsucc/151023_0002_BSA_0pct_light_sw1.mzXML",
                       sky.report=set0pct,
                       ppm=20,mzrange=c(100,1500),
                       threshold=30)
stoich1<-stoichwrapper(mzxml="G:/tmp/BSAsucc/151023_0002_BSA_1pct_light_sw1.mzXML",
                       sky.report=set1pct,
                       ppm=20,mzrange=c(100,1500),
                       threshold=30)
stoich10<-stoichwrapper(mzxml="G:/tmp/BSAsucc/151023_0002_BSA_10pct_light_sw1.mzXML",
                        sky.report=set10pct,
                        ppm=20,mzrange=c(100,1500),
                        threshold=30)
stoich50<-stoichwrapper(mzxml="G:/tmp/BSAsucc/151023_0002_BSA_50pct_light_sw1.mzXML",
                        sky.report=set50pct,
                        ppm=20,mzrange=c(100,1500),
                        threshold=30)
stoich100<-stoichwrapper(mzxml="G:/tmp/BSAsucc/151023_0002_BSA_100pct_light_sw1.mzXML",
                        sky.report=set100pct,
                        ppm=20,mzrange=c(100,1500),
                        threshold=30)


return.lmrat=function(object=temp.peplist[[1]]){
  object@lm.ratio
}

return.median.ratio=function(object=temp.peplist[[1]]){
  object@median.ratio
}

return.rank1.ratio=function(object=temp.peplist[[1]]){
  object@rank1.ratio
}

return.top3.ratio=function(object=temp.peplist[[1]]){
  object@top3.ratio
}


return.both=function(object=temp.peplist[[1]]){
  return(c(object@median.ratio,object@lm.ratio,as.numeric(object@median.ratio)-as.numeric(object@lm.ratio)))
}

lapply(FUN=return.median.ratio,stoich50)->medians
lapply(FUN=return.lmrat,stoich50)->lmratios
lapply(FUN=return.rank1.ratio,stoich50)->rank1
lapply(FUN=return.top3.ratio,stoich50)->top3


ave(na.omit(as.numeric(unlist(medians))))
sd(na.omit(as.numeric(unlist(medians))))

ave(na.omit(as.numeric(unlist(lmratios))))
sd(na.omit(as.numeric(unlist(lmratios))))

ave(na.omit(as.numeric(unlist(rank1))))
sd(na.omit(as.numeric(unlist(rank1))))

ave(na.omit(as.numeric(unlist(top3))))
sd(na.omit(as.numeric(unlist(top3))))

dev.off()
par()
boxplot(na.omit(as.numeric(unlist(medians))),
        na.omit(as.numeric(unlist(lmratios))),
        na.omit(as.numeric(unlist(rank1))),
        na.omit(as.numeric(unlist(top3))),
        ylim=c(0,0.75),lwd=2
)
abline(h=0.5,lty=2)

temp.peplist<-stoich10
allratios<-list()
for(j in 1:length(temp.peplist)){
  allratios[[j]]<-temp.peplist[[j]]@median.ratio
}







#temp.peplist1<-temp.peplist()
y=rank(temp.peplist[[41]]@areas$position1$light)
x=temp.peplist[[41]]@areas$position1$heavy
#plot(x,y)
model<-lm(y~x)
#abline(model)
coef(model)[2]
data.frame(x,y+0)
model<-lm(x~y+0)
model[3]
plot(coef(model))
plot(model)
medians<-na.omit(as.numeric(unlist(allratios)))
ave(na.omit(as.numeric(unlist(allratios))))
plot(medians)
sd(na.omit(as.numeric(unlist(allratios))))

#plotmax<-max(as.numeric(na.omit(unlist(xicmat[,2:ncol(xicmat)]))))+max(as.numeric(na.omit(unlist(xicmat[,2:ncol(xicmat)]))))*0.05

counter=1
plot(xicmat[,1],xicmat[,2],ylim=c(0,plotmax),t="l",lwd=2,col=counter,ylab="intensity",xlab="time")

for(x in seq(from=3,to=ncol(xicmat))){
  counter=counter+1
  lines(xicmat[,1],xicmat[,x],lwd=2,col=counter) 
}
### part to plot those traces and save a file
counter=1
plot(xics[[1]][[1]],ylim=c(0,max(na.omit(unlist(xics)))+0.05*max(na.omit(unlist(xics)))),type="l",col="black",lwd=2)
lapply(FUN=lines,xics[[2]],type="l",col=counter)
}

}




p1$prec$LL[set10pct[testline,"Precursor.Charge"]]
p1$prec$LH[set10pct[testline,"Precursor.Charge"]]
p1$prec$HH[set10pct[testline,"Precursor.Charge"]]
