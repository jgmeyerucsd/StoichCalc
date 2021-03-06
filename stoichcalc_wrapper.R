s<-read.delim(file="stoichout1test1.tsv", head=T,stringsAsFactors=F)
s2<-read.delim(file="stoichout1test2.tsv", head=T,stringsAsFactors=F)

setwd("~/stoichiometry")
### read in skyline report
s<-read.delim(file="stoichout1test5.tsv", head=T,stringsAsFactors=F)


s<-unique(s)
peptides<-unique(s[,1])
head(s)

mzxmlfile="G:/tmp/BSAsucc/151023_0002_BSA_10pct_light_sw1.mzXML"

set10pct<-s[which(s[,"File.Name"]=="151023_0002_BSA_10pct_light_sw1.wiff"),]
testline<-which(set10pct[,1]==sequence)


x<-peptides[1]
p1<-getIons(sequence=peptides[9],modmass=100.016044,moddelta=4.025107)

setClass("Peptide", 
         slots=list(Kcount="integer",
                    sequence="character",
                    modpos="integer",
                    ionlist="list",
                    areas="ANY",
                    areas.filtered="ANY",
                    median.ratio="ANY",
                    peakbounds="vector",
                    prec.z="vector"
                    ))



stoichwrapper=function(mzxml="G:/tmp/BSAsucc/151023_0002_BSA_10pct_light_sw1.mzXML",
  sky.report=set10pct,
  wind.file="windoe.txt",
  windowtable="windoe.txt",
  ppm=20,
  threshold=10,
  mzmax=1500)
  {

  ###### get peptide info and make an object of peptide class with info
  peptides<-unique(sky.report[,"Peptide"])
  peplist<-list()
  for(k in peptides){
    #x<-peptides[15]
    print(k)
    temppep<-new(Class="Peptide")
    temppep@ionlist<-getIons(sequence=k)
    temppep@sequence<-unlist(strsplit(k,split=""))
    temppep@modpos<-which(temppep@sequence=="K")
    temppep@Kcount<-length(temppep@modpos)
    templines<-which(sky.report[,1]==k)
    temppep@prec.z<-sky.report[templines,"Precursor.Charge"]
    temppep@peakbounds<-c(sky.report[templines[1],"Start.Time"],sky.report[templines[1],"End.Time"])
    
    ####################################################################################################################
    ### build a list of listed XICs 
    ##### for case with single lys and double lys separtely
    ####################################################################################################################
    xics<-list()
    xicmat<-list()
    if(temppep@Kcount==1){
      ### light XICs
      counter=1
      for(y in temppep@prec.z){
        fragments<-na.omit(unlist(temppep@ionlist$L))
        fragments<-fragments[fragments<=mzmax]
        nfrag<-length(fragments)
        ### reset counter to 1
        counter=1
       # xics[[y]]<-y
        for(z in fragments){
          xics$light[[paste(y)]][[paste(z)]]<-ms2xic(file=mzxml,precMz=temppep@ionlist$prec$L[y], fragMz=z, type="l",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm=ppm)
          if(counter==1){
            ### part to plot while extracting only works if max height is known, fix later
            xicmat[[paste(y)]]<-data.frame(xics$light[[paste(y)]][[paste(z)]])
            colnames(xicmat[[paste(y)]])[2]<-paste("prec=",round(temppep@ionlist$prec$L[y],digits=1),", ","frag=",z,collapse="_",sep="")
            plot(xics$light[[paste(y)]][[paste(z)]][,1],xics$light[[paste(y)]][[paste(z)]][,2],ylim=c(0,40),type="l",col=counter,lwd=2)
            
          }
          if(counter!=1){
            xicmat[[paste(y)]]<-cbind(xicmat[[paste(y)]],xics$light[[paste(y)]][[paste(z)]][,2])
            colnames(xicmat[[paste(y)]])[counter+1]<-paste("prec=",round(temppep@ionlist$prec$L[y],digits=1),", ","frag=",z,collapse="_",sep="")
            lines(xics$light[[paste(y)]][[paste(z)]][,1],xics$light[[paste(y)]][[paste(z)]][,2],ylim=c(0,40),col=counter,lwd=2)
          }
          counter=counter+1
        }
      }
      tempcounter<-counter
      #######################################################################
      #### heavy XICs, single lysine
      ########################################################################
      for(y in temppep@prec.z){
        fragments<-na.omit(unlist(temppep@ionlist$H))
        fragments<-fragments[fragments<=1500]
        counter<-tempcounter
        for(z in fragments){
          xics$heavy[[paste(y)]][[paste(z)]]<-ms2xic(file=mzxml,precMz=temppep@ionlist$prec$H[y], fragMz=z, type="h",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm=ppm)
          xicmat[[paste(y)]]<-cbind(xicmat[[paste(y)]],xics$heavy[[paste(y)]][[paste(z)]][,2])
          colnames(xicmat[[paste(y)]])[counter+1]<-paste("prec=",round(temppep@ionlist$prec$H[y],digits=1),", ","frag=",z,collapse="_",sep="")
          counter=counter+1
        }
      }
    ### TODO: copy the part above to get a sense of the average noise per scan from each transition
    #################################################################################################################
    ### compute all the fragment to fragment ratios
    ###################################################################################################################
      
      sumlight<-c()
      sumheavy<-c()
      ratio<-c()
      for(i in temppep@prec.z){
        for(x in 1:nfrag){
          print(x)
          sumlight<-c(sumlight,sum(na.omit(xicmat[[paste(i)]][,x+1])))
          sumheavy<-c(sumheavy,sum(na.omit(xicmat[[paste(i)]][,x+nfrag+1])))
          ratio<-c(ratio,sum(na.omit(xicmat[[paste(i)]][,x+1]))/(sum(na.omit(xicmat[[paste(i)]][,x+1]))+sum(na.omit(xicmat[[paste(i)]][,x+nfrag+1]))))
        }
      }
      median(na.omit(ratio))
      #################################################
      ### filter based on which values are at least 10 
      sumheavy2<-sumheavy[which(sumlight>=threshold)]
      sumlight2<-sumlight[which(sumlight>=threshold)]
      ratios<-sumlight2/(sumlight2+sumheavy2)
      
      #### return the media value of all ratios    
      temppep@areas[["position1"]]<-list(light=sumlight,heavy=sumheavy,ratio=ratio)
      temppep@median.ratio[["position1"]]<-median(na.omit(ratios))
    
    
    print("singleK peptide finished")
    }
    
    
    ################################################
    ##### part for the peptides with 2 lysine
    ######################################################
    #####################################################################################################################################
    #### need to fix:  now there is a problem that sometimes the different precursor will have one more or less scan/RT, decide how to fix
    ###############################################################
    ###   either trim the newer one, or add zeros to new one, or revert to list form
    
    xics<-list()
    b_diff.xicmat<-list()
    if(temppep@Kcount==2){
      
      #### split into b_diff and y_diff
      #### first extract b_diff
      for(y in temppep@prec.z){
        ###   get fragments ions for b_diff to extract
        fragments<-na.omit(unlist(temppep@ionlist$LL$b_diff))
        fragments<-fragments[fragments<=mzmax]
        nfrag<-length(fragments)
        ### reset counter to 1
        counter=1
        #b_diff.xicmat[[paste(y)]]<-data.frame()
        for(z in fragments){
          #xicdataframe<-cbind(xicdataframe,ms2xic(file=mzxml,precMz=temppep@ionlist$prec$L[y], fragMz=z, rtrange=temppep@peakbounds))
          xics$light[[paste(y)]][[paste(z)]]<-ms2xic(file=mzxml,precMz=temppep@ionlist$prec$LL[y], fragMz=z, type="l",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount,ppm= ppm)
          if(counter==1){
            ### part to plot while extracting, only works if max height is known
            b_diff.xicmat[[paste(y)]]<-data.frame(xics$light[[paste(y)]][[paste(z)]])
            colnames(b_diff.xicmat[[paste(y)]])[2]<-paste("prec=",round(temppep@ionlist$prec$LL[y],digits=1),", ","frag=",z,collapse="_",sep="")
            plot(xics$light[[paste(y)]][[paste(z)]][,1],xics$light[[paste(y)]][[paste(z)]][,2],ylim=c(0,40),type="l",col=counter,lwd=2)
          }
          if(counter!=1){
            b_diff.xicmat[[paste(y)]]<-cbind(b_diff.xicmat[[paste(y)]],test=xics$light[[paste(y)]][[paste(z)]][,2],deparse.level = 2)
            colnames(b_diff.xicmat[[paste(y)]])[counter+1]<-paste("prec=",round(temppep@ionlist$prec$LL[y],digits=1),", ","frag=",z,collapse="_",sep="")
            lines(xics$light[[paste(y)]][[paste(z)]][,1],xics$light[[paste(y)]][[paste(z)]][,2],ylim=c(0,40),col=counter,lwd=2)
          }
          counter=counter+1
        }
      }
      b_diff.nlight<-(ncol(b_diff.xicmat[[1]])-1)*length(temppep@prec.z)
      tempcount<-counter
      for(y in temppep@prec.z){
        fragments<-na.omit(unlist(temppep@ionlist$HH$b_diff))
        # xics[[y]]<-y
        fragments<-fragments[fragments<=mzmax]
        counter=tempcount
        for(z in fragments){
          
          xics$heavy[[paste(y)]][[paste(z)]]<-ms2xic(file=mzxml,precMz=temppep@ionlist$prec$HH[y], fragMz=z, type="h",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount,ppm= ppm)
          b_diff.xicmat[[paste(y)]]<-cbind(b_diff.xicmat[[paste(y)]],xics$heavy[[paste(y)]][[paste(z)]][,2])
          colnames(b_diff.xicmat[[paste(y)]])[counter+1]<-paste("prec=",round(temppep@ionlist$prec$HH[y],digits=1),", ","frag=",z,collapse="_",sep="")
          counter=counter+1
        }
        
      }
      
      ##################################################################################################################
      #### calculate the occupancy at first site
      ### compute all the fragment to fragment ratios
      sumlight<-c()
      sumheavy<-c()
      ratio<-c()
      for(i in temppep@prec.z){
        for(x in 1:nfrag){
          print(x)
          sumlight<-c(sumlight,sum(na.omit(b_diff.xicmat[[paste(i)]][,x+1])))
          sumheavy<-c(sumheavy,sum(na.omit(b_diff.xicmat[[paste(i)]][,x+nfrag+1])))
          ratio<-c(ratio,sum(na.omit(b_diff.xicmat[[paste(i)]][,x+1]))/(sum(na.omit(b_diff.xicmat[[paste(i)]][,x+1]))+sum(na.omit(b_diff.xicmat[[paste(i)]][,x+nfrag+1]))))
          }
        }
      
      #################################################
      ### filter based on which values are at least 10 
      sumheavy2<-sumheavy[which(sumlight>=threshold)]
      sumlight2<-sumlight[which(sumlight>=threshold)]
      ratios<-sumlight2/(sumlight2+sumheavy2)
      #### return the media value of all ratios    
      temppep@areas[["position1"]]<-list(light=sumlight,heavy=sumheavy,ratio=ratio)
      temppep@median.ratio[["position1"]]<-median(na.omit(ratios))
      
      
      
      #### next extract y_diff
      
      xics<-list()
      y_diff.xicmat<-list()
      #### split into b_diff and y_diff
      #### first extract b_diff
      counter=1
      for(y in temppep@prec.z){
        ###   get fragments ions for b_diff to extract
        fragments<-na.omit(unlist(temppep@ionlist$LL$y_diff))
        fragments<-fragments[fragments<=mzmax]
        nfrag<-length(fragments)
        counter=1
        for(z in fragments){
          #xicdataframe<-cbind(xicdataframe,ms2xic(file=mzxml,precMz=temppep@ionlist$prec$L[y], fragMz=z, rtrange=temppep@peakbounds))
          xics$light[[paste(y)]][[paste(z)]]<-ms2xic(file=mzxml,precMz=temppep@ionlist$prec$LL[y], fragMz=z, type="l",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm= ppm)
          if(counter==1){
            ### part to plot while extracting, only works if max height is known
            y_diff.xicmat[[paste(y)]]<-data.frame(xics$light[[paste(y)]][[paste(z)]])
            colnames(y_diff.xicmat[[paste(y)]])[2]<-paste("prec=",round(temppep@ionlist$prec$LL[y],digits=1),", ","frag=",z,collapse="_",sep="")
            plot(xics$light[[paste(y)]][[paste(z)]][,1],xics$light[[paste(y)]][[paste(z)]][,2],ylim=c(0,40),type="l",col=counter,lwd=2)
            
          }
          if(counter!=1){
            y_diff.xicmat[[paste(y)]]<-cbind(y_diff.xicmat[[paste(y)]],xics$light[[paste(y)]][[paste(z)]][,2])
            colnames(y_diff.xicmat[[paste(y)]])[counter+1]<-paste("prec=",round(temppep@ionlist$prec$LL[y],digits=1),", ","frag=",z,collapse="_",sep="")
            lines(xics$light[[paste(y)]][[paste(z)]][,1],xics$light[[paste(y)]][[paste(z)]][,2],ylim=c(0,40),col=counter,lwd=2)
          }
          counter=counter+1
        }
      }
      #y_diff.nlight<-ncol(y_diff.xicmat)-1
      tempcount<-counter
      for(y in temppep@prec.z){
        fragments<-na.omit(unlist(temppep@ionlist$HH$y_diff))
        # xics[[y]]<-y
        fragments<-fragments[fragments<=1500]
        counter<-tempcount
        for(z in fragments){
          
          xics$heavy[[paste(y)]][[paste(z)]]<-ms2xic(file=mzxml,precMz=temppep@ionlist$prec$HH[y], fragMz=z, type="h",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm= ppm)
          y_diff.xicmat[[paste(y)]]<-cbind(y_diff.xicmat[[paste(y)]],xics$heavy[[paste(y)]][[paste(z)]][,2])
          colnames(y_diff.xicmat[[paste(y)]])[counter+1]<-paste("prec=",round(temppep@ionlist$prec$HH[y],digits=1),", ","frag=",z,collapse="_",sep="")
          counter=counter+1
        }
        
      }
      
      ##################################################################################################################
      #### calculate the occupancy at first site
      ### compute all the fragment to fragment ratios
      sumlight<-c()
      sumheavy<-c()
      ratio<-c()
      for(i in temppep@prec.z){
        for(x in 1:nfrag){
          print(x)
          sumlight<-c(sumlight,sum(na.omit(y_diff.xicmat[[paste(i)]][,x+1])))
          sumheavy<-c(sumheavy,sum(na.omit(y_diff.xicmat[[paste(i)]][,x+nfrag+1])))
          ratio<-c(ratio,sum(na.omit(y_diff.xicmat[[paste(i)]][,x+1]))/(sum(na.omit(y_diff.xicmat[[paste(i)]][,x+1]))+sum(na.omit(y_diff.xicmat[[paste(i)]][,x+nfrag+1]))))
        }
      }
      
      #################################################
      ### filter based on which values are at least 10 
      sumheavy2<-sumheavy[which(sumlight>=threshold)]
      sumlight2<-sumlight[which(sumlight>=threshold)]
      ratios<-sumlight2/(sumlight2+sumheavy2)
      #### return the media value of all ratios    
      temppep@areas[["position2"]]<-list(light=sumlight,heavy=sumheavy,ratio=ratio)
      temppep@median.ratio[["position2"]]<-median(na.omit(ratios))
      print(k)
      print("double K peptide finished")
      
      
    }
    if(temppep@Kcount==3){
      print(k)
      print("3+ lys not yet supported")
      temppep@median.ratio<-"three lys not yet supported"
      }
    ### store all the peptides in a list and return that list
    print(k)
    peplist[[paste(k)]]<-temppep
  }
    temp.peplist<-peplist
    unlist(temp.peplist)@median.ratio
    length(temp.peplist)
    allratios<-list()
    for(j in 1:length(temp.peplist)){
      allratios[[j]]<-temp.peplist[[j]]@median.ratio
    }
    temp.peplist[[7]]@sequence
    
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
