setwd("~/stoichiometry")
### read in skyline report
s<-read.delim(file="stoichout1test7.tsv", head=T,stringsAsFactors=F)

head(s)


#### function to fill up basic character of peptide objects



skyline2peptide=function(filename="stoichout1test8.tsv"){
  s<-read.delim(file=filename, head=T,stringsAsFactors=F)
  sreport<-new(Class="skylineReport")
  sreport@peptides<-list()
  pepvec<-unique(s[,"Peptide"])
  for(k in pepvec){
    #k<-pepvec[1]
    print(k)
    pepcount=pepcount+1
    print(pepcount)
    temppep<-new(Class="Peptide")
    temppep@ionlist<-getIons(sequence=k)
    temppep@sequence<-unlist(strsplit(k,split=""))
    temppep@modpos<-which(temppep@sequence=="K")
    temppep@Kcount<-length(temppep@modpos)

    temppep@prec.z<-unique(s[templines,"Precursor.Charge"])
    for(x in temppep@prec.z){
      temppep@reportlines[[x]]<-which(s[,"Peptide"]==k & s[,"Precursor.Charge"]==x)
    }
    
    

    temppep@peakbounds<-c(sky.report[templines[1],"Start.Time"],sky.report[templines[1],"End.Time"])
    temppep@protein.name<-sky.report[templines[1],"Protein"]
    temppep@protein.position<-sky.report[templines[1],"Begin.Pos"]+temppep@modpos-1
    
    ### determine which ion is rank 1 for each precursor charge
    #### currently just takes the first line from the skyline report that has rank 1 for that precursor charge
    
    for(x in temppep@prec.z){
      temppep@reportlines[[x]][which(as.numeric(s[temppep@reportlines[[x]],"Library.Rank"])==1)]
      temppep@rank1ion[[x]]<-s[temppep@reportlines[[x]][which(as.numeric(s[temppep@reportlines[[x]],"Library.Rank"])==1)][1],]
    }
    sreport@peptides[[paste(k)]]<-temppep
  }
  
  
  
}
  
s<-unique(s)
peptides<-unique(s[,1])
head(s)

mzxmlfile="G:/tmp/BSAsucc/151023_0002_BSA_10pct_light_sw1.mzXML"
set100pct<-s[which(s[,"File.Name"]=="151023_0002_BSA_100pct_light_sw1.mzXML"),]
set50pct<-s[which(s[,"File.Name"]=="151023_0002_BSA_50pct_light_sw1.mzXML"),]
set10pct<-s[which(s[,"File.Name"]=="151023_0002_BSA_10pct_light_sw1.mzXML"),]
set1pct<-s[which(s[,"File.Name"]=="151023_0002_BSA_1pct_light_sw1.mzXML"),]
set0pct<-s[which(s[,"File.Name"]=="151023_0002_BSA_0pct_light_sw1.mzXML"),]
testline<-which(set10pct[,1]==sequence)

##### almost works for any peptide, finish function after testing
stoich1<-stoichwrapper(mzxml="G:/tmp/BSAsucc/151023_0002_BSA_1pct_light_sw1.mzXML",
                       sky.report=set1pct,
                       ppm=20,mzrange=c(100,1500),
                       threshold=30)

stoichwrapper=function(mzxml="G:/tmp/BSAsucc/151023_0002_BSA_10pct_light_sw1.mzXML",
                       sky.report=set10pct,
                       windowtab=read.delim(file="windoe.txt",head=F,stringsAsFactors=F),
                       ppm=30,
                       threshold=30,
                       mzrange=c(500,1500))
{
  #system.time(ms2xic(file=mzxml,precMz=temppep@ionlist$prec$LL[y], fragMz=600, type="l",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm=10))
  ###### get peptide info and make an object of peptide class with info
  peptides<-unique(sky.report[,"Peptide"])
  peplist<-list()
  #windowtab<-read.delim(file=windowtable,head=F,stringsAsFactors=F)
  require(mzR)
  mslink <- openMSfile(filename=mzxml)
  hdlink <-	header(mslink)
  pepcount=0
  for(k in peptides){
    #k<-peptides[6]
    print(k)
    pepcount=pepcount+1
    print(pepcount)
    temppep<-new(Class="Peptide")
    temppep@ionlist<-getIons(sequence=k)
    temppep@sequence<-unlist(strsplit(k,split=""))
    temppep@modpos<-which(temppep@sequence=="K")
    temppep@Kcount<-length(temppep@modpos)
    templines<-which(sky.report[,1]==k)
    temppep@prec.z<-sky.report[templines,"Precursor.Charge"]
    temppep@peakbounds<-c(sky.report[templines[1],"Start.Time"],sky.report[templines[1],"End.Time"])
    temppep@protein.name<-sky.report[templines[1],"Protein"]
    temppep@protein.position<-sky.report[templines[1],"Begin.Pos"]+temppep@modpos-1
    
    ##### check if at least one precursor fall is in the same window
    ##### replace temppep@prec.z with those that work
    
    use<-c()
    if(temppep@Kcount<3){
      for(x in temppep@prec.z){
        print(temppep@ionlist$prec[[1]][x])
        print(x)
        use<-c(use,testprec(precMz=temppep@ionlist$prec[[1]][x],precZ = x))
      }
    }
    temppep@prec.z<-use
    
    ####################################################################################################################
    ### build a list of listed XICs 
    ##### for case with single lys and double lys separtely
    ####################################################################################################################
    xics<-list()
    xicmat<-list()
    nfrag<-c()
    if(temppep@Kcount==1){
      ### light XICs
      counter=1
      for(y in temppep@prec.z){
        fragments<-na.omit(unlist(temppep@ionlist$L))
        fragments<-fragments[fragments<=mzrange[2]]
        nfrag[y]<-length(fragments)
        bfrag<-unlist(temppep@ionlist$L$b_diff)
        bfrag.ord<-which(bfrag!="NA")[1:(length(which(bfrag!="NA"))/2)]
        yfrag<-unlist(temppep@ionlist$L$y_diff)
        yfrag.ord<-which(yfrag!="NA")[1:(length(which(yfrag!="NA"))/2)]
        temppep@iontypes<-list(bions=bfrag.ord,yions=yfrag.ord)
        #nbfrag<-length(bfrag)
        ### reset counter to 1
        counter=1
        # xics[[y]]<-y
        for(z in fragments){
          xics$light[[paste(y)]][[paste(z)]]<-ms2xic(ms=mslink,hd=hdlink,precMz=temppep@ionlist$prec$L[y], fragMz=z, type="l",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm=ppm)
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
      
      #######################################################################
      #### heavy XICs, single lysine
      ########################################################################
      for(y in temppep@prec.z){
        fragments<-na.omit(unlist(temppep@ionlist$H))
        fragments<-fragments[ fragments<=mzrange[2]]
        counter<-length(xicmat[[paste(y)]])
        for(z in fragments){
          xics$heavy[[paste(y)]][[paste(z)]]<-ms2xic(ms=mslink,hd=hdlink,file=mzxml,precMz=temppep@ionlist$prec$H[y], fragMz=z, type="h",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm=ppm)
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
      sumlightnames<-c()
      sumheavynames<-c()
      ratio<-c()
      for(i in temppep@prec.z){
        tempnames<-names(xicmat[[paste(i)]])
        for(x in 1:nfrag[i]){
          print(x)
          sumlight<-c(sumlight,sum(na.omit(xicmat[[paste(i)]][,x+1])))
          sumheavy<-c(sumheavy,sum(na.omit(xicmat[[paste(i)]][,x+nfrag[i]+1])))
          sumlightnames<-c(sumlightnames,tempnames[x+1])
          sumheavynames<-c(sumheavynames,tempnames[x+nfrag[i]+1])
          ratio<-c(ratio,sum(na.omit(xicmat[[paste(i)]][,x+1]))/(sum(na.omit(xicmat[[paste(i)]][,x+1]))+sum(na.omit(xicmat[[paste(i)]][,x+nfrag[i]+1]))))
        }
      }
      #################################################
      ### filter based on which values are at least threshold
      sumheavy2<-sumheavy[which(sumheavy>=threshold)]
      sumlight2<-sumlight[which(sumheavy>=threshold)]
      
      ### vector of ion names for reporting
      nprec<-length(temppep@prec.z)
      ionnames<-rep(c(paste("b",temppep@iontypes$bions,sep=""),paste("b",temppep@iontypes$bions,sep="","++"),paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,"++",sep="")),times=nprec)
      names(sumheavy)<-sumheavynames
      names(sumlight)<-sumlightnames
      
      #### get and store ion names and masses
      names(sumlight2)<-sumlightnames[which(sumheavy>=threshold)]
      names(sumheavy2)<-sumheavynames[which(sumheavy>=threshold)]
      ordinals<-ionnames[which(sumheavy>=threshold)]
      temppep@iontypes[["position1"]]<-list(ordinals=ordinals,massL=names(sumlight2),massH=names(sumheavy2))
      ratios<-sumlight2/(sumlight2+sumheavy2)
      #### return the media value of all ratios
      temppep@areas[["position1"]]<-list(light=sumlight,heavy=sumheavy,ratio=ratio)
      temppep@areas.filtered[["position1"]]<-list(light=sumlight2,heavy=sumheavy2,ratio=ratios)
      
      temppep@median.ratio[["position1"]]<-median(round(na.omit(ratios),digits=4))
      temppep@median.ratio.sd<-sd(na.omit(ratios))
      if(length(ratio)>=1){
        temppep@median.ratio.unfilt[["position1"]]<-median(round(na.omit(ratio),digits=4))
      }
      ### unused part for top3 and linear model stoich calc
      if(length(sumheavy)>2){
        y=sumlight
        x=sumheavy
        model<-lm(y~x)
        temppep@lm.ratio[["position1"]]<-coef(model)[2]
        temppep@top3.ratio[["position1"]]<-ave(sumlight[order(sumheavy,decreasing=T)][1:3]/(sumheavy[order(sumheavy,decreasing=T)][1:3]+sumlight[order(sumheavy,decreasing=T)][1:3]))[1]
      }
      #### heavy rank1 ratios
      temppep@heavy.rank1.ratio[["position1"]]<-list(filtered=sumlight2[sumheavy2==max(sumheavy2)]/(sumheavy2[sumheavy2==max(sumheavy2)]+sumlight2[sumheavy2==max(sumheavy2)]))
      temppep@heavy.rank1.ratio[["position1"]][["unfiltered"]]<-sumlight[sumheavy==max(sumheavy)]/(sumheavy[sumheavy==max(sumheavy)]+sumlight[sumheavy==max(sumheavy)])
      ##### light rank1 ratios
      temppep@light.rank1.ratio[["position1"]]<-list(filtered=sumlight2[sumlight2==max(sumlight2)]/(sumheavy2[sumlight2==max(sumlight2)]+sumlight2[sumlight2==max(sumlight2)]))
      temppep@light.rank1.ratio[["position1"]][["unfiltered"]]<-sumlight[sumlight==max(sumlight)]/(sumheavy[sumlight==max(sumlight)]+sumlight[sumlight==max(sumlight)])
      
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
        ### test if precursor is in two windows, otherwise skip
        ###   get fragments ions for b_diff to extract
        fragments<-na.omit(unlist(temppep@ionlist$LL$b_diff))
        fragments<-fragments[fragments<=mzrange[2]]
        nfrag[y]<-length(fragments)
        nfrag[y]<-length(fragments)
        bfrag<-unlist(temppep@ionlist$LL$b_diff)
        bfrag.ord<-which(bfrag!="NA")[1:(length(which(bfrag!="NA"))/2)]
        yfrag<-unlist(temppep@ionlist$LL$y_diff)
        yfrag.ord<-which(yfrag!="NA")[1:(length(which(yfrag!="NA"))/2)]
        temppep@iontypes<-list(bions=bfrag.ord,yions=yfrag.ord)
        ### set counter to 1
        counter=1
        #b_diff.xicmat[[paste(y)]]<-data.frame()
        for(z in fragments){
          #xicdataframe<-cbind(xicdataframe,ms2xic(file=mzxml,precMz=temppep@ionlist$prec$L[y], fragMz=z, rtrange=temppep@peakbounds))
          xics$light[[paste(y)]][[paste(z)]]<-ms2xic(ms=mslink,hd=hdlink,file=mzxml,precMz=temppep@ionlist$prec$LL[y], fragMz=z, type="l",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount,ppm= ppm)
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
      
      for(y in temppep@prec.z){
        fragments<-na.omit(unlist(temppep@ionlist$HH$b_diff))
        # xics[[y]]<-y
        fragments<-fragments[fragments<=mzrange[2]]
        counter<-length(b_diff.xicmat[[paste(y)]])
        for(z in fragments){
          
          xics$heavy[[paste(y)]][[paste(z)]]<-ms2xic(ms=mslink,hd=hdlink,file=mzxml,precMz=temppep@ionlist$prec$HH[y], fragMz=z, type="h",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount,ppm= ppm)
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
      sumlightnames<-c()
      sumheavynames<-c()
      ratio<-c()
      for(i in temppep@prec.z){
        tempnames<-names(b_diff.xicmat[[paste(i)]])
        if(length(b_diff.xicmat[[paste(i)]])>0){
          for(x in 1:nfrag[i]){
            print(x)
            sumlight<-c(sumlight,sum(na.omit(b_diff.xicmat[[paste(i)]][,x+1])))
            sumheavy<-c(sumheavy,sum(na.omit(b_diff.xicmat[[paste(i)]][,x+nfrag[i]+1])))
            sumlightnames<-c(sumlightnames,tempnames[x+1])
            sumheavynames<-c(sumheavynames,tempnames[x+nfrag[i]+1])
            ratio<-c(ratio,sum(na.omit(b_diff.xicmat[[paste(i)]][,x+1]))/(sum(na.omit(b_diff.xicmat[[paste(i)]][,x+1]))+sum(na.omit(b_diff.xicmat[[paste(i)]][,x+nfrag[i]+1]))))
          }
        }
      }
      
      #################################################
      ### filter based on which values are at least 10 
      sumheavy2<-sumheavy[which(sumheavy>=threshold)]
      sumlight2<-sumlight[which(sumheavy>=threshold)]
      
      ### vector of ion names for reporting
      nprec<-length(temppep@prec.z)
      ionnames<-rep(c(paste("b",temppep@iontypes$bions,sep=""),paste("b",temppep@iontypes$bions,sep="","++"),paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,"++",sep="")),times=nprec)
      names(sumheavy)<-sumheavynames
      names(sumlight)<-sumlightnames
      
      #### get and store ion names and masses
      names(sumlight2)<-sumlightnames[which(sumheavy>=threshold)]
      names(sumheavy2)<-sumheavynames[which(sumheavy>=threshold)]
      ordinals<-ionnames[which(sumheavy>=threshold)]
      temppep@iontypes[["position1"]]<-list(ordinals=ordinals,massL=names(sumlight2),massH=names(sumheavy2))
      
      ratios<-sumlight2/(sumlight2+sumheavy2)
      #### return the median value of all ratios    
      temppep@areas[["position1"]]<-list(light=sumlight,heavy=sumheavy,ratio=ratio)
      temppep@areas.filtered[["position1"]]<-list(light=sumlight2,heavy=sumheavy2,ratio=ratios)
      temppep@median.ratio[["position1"]]<-round(median(na.omit(ratios)),digits = 4)
      length(ratio)
      if(length(ratio)>=1){
        temppep@median.ratio.unfilt[["position1"]]<-median(round(na.omit(ratio),digits=4))
      }
      #temppep@median.ratio.unfilt[["position1"]]<-median(round(na.omit(ratio),digits=4))
      #if(length(sumheavy)>2){
      #  
      #  y=sumlight
      #  x=sumheavy
      #  model<-lm(y~x)
      #  temppep@lm.ratio[["position1"]]<-coef(model)[2]
      #  temppep@top3.ratio[["position1"]]<-ave(sumlight[order(sumheavy,decreasing=T)][1:3]/(sumheavy[order(sumheavy,decreasing=T)][1:3]+sumlight[order(sumheavy,decreasing=T)][1:3]))[1]
      #  
      #}
      #### heavy rank 1 ratios
      temppep@heavy.rank1.ratio[["position1"]]<-list(filtered=sumlight2[sumheavy2==max(sumheavy2)]/(sumheavy2[sumheavy2==max(sumheavy2)]+sumlight2[sumheavy2==max(sumheavy2)]))
      temppep@heavy.rank1.ratio[["position1"]][["unfiltered"]]<-sumlight[sumheavy==max(sumheavy)]/(sumheavy[sumheavy==max(sumheavy)]+sumlight[sumheavy==max(sumheavy)])
      ##### light rank1 ratios
      temppep@light.rank1.ratio[["position1"]]<-list(filtered=sumlight2[sumlight2==max(sumlight2)]/(sumheavy2[sumlight2==max(sumlight2)]+sumlight2[sumlight2==max(sumlight2)]))
      temppep@light.rank1.ratio[["position1"]][["unfiltered"]]<-sumlight[sumlight==max(sumlight)]/(sumheavy[sumlight==max(sumlight)]+sumlight[sumlight==max(sumlight)])
      
      #### next extract y_diff
      
      xics<-list()
      y_diff.xicmat<-list()
      #### split into b_diff and y_diff
      #### first extract b_diff
      counter=1
      for(y in temppep@prec.z){
        ###   get fragments ions for b_diff to extract
        fragments<-na.omit(unlist(temppep@ionlist$LL$y_diff))
        fragments<-fragments[ fragments<=mzrange[2]]
        nfrag[y]<-length(fragments)
        counter=1
        for(z in fragments){
          #xicdataframe<-cbind(xicdataframe,ms2xic(file=mzxml,precMz=temppep@ionlist$prec$L[y], fragMz=z, rtrange=temppep@peakbounds))
          xics$light[[paste(y)]][[paste(z)]]<-ms2xic(ms=mslink,hd=hdlink,file=mzxml,precMz=temppep@ionlist$prec$LL[y], fragMz=z, type="l",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm= ppm)
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
        fragments<-fragments[ fragments<=mzrange[2]]
        counter<-length(y_diff.xicmat[[paste(y)]])
        for(z in fragments){
          
          xics$heavy[[paste(y)]][[paste(z)]]<-ms2xic(ms=mslink,hd=hdlink,file=mzxml,precMz=temppep@ionlist$prec$HH[y], fragMz=z, type="h",precZ=y, rtrange=temppep@peakbounds, Kcount=temppep@Kcount, ppm= ppm)
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
      sumlightnames<-c()
      sumheavynames<-c()
      ratio<-c()
      for(i in temppep@prec.z){
        tempnames<-names(y_diff.xicmat[[paste(i)]])
        if(length(y_diff.xicmat[[paste(i)]])>0){
          for(x in 1:nfrag[i]){
            print(x)
            sumlight<-c(sumlight,sum(na.omit(y_diff.xicmat[[paste(i)]][,x+1])))
            sumheavy<-c(sumheavy,sum(na.omit(y_diff.xicmat[[paste(i)]][,x+nfrag[i]+1])))
            sumlightnames<-c(sumlightnames,tempnames[x+1])
            sumheavynames<-c(sumheavynames,tempnames[x+nfrag[i]+1])
            ratio<-c(ratio,sum(na.omit(y_diff.xicmat[[paste(i)]][,x+1]))/(sum(na.omit(y_diff.xicmat[[paste(i)]][,x+1]))+sum(na.omit(y_diff.xicmat[[paste(i)]][,x+nfrag[i]+1]))))
          }
        }
      }
      
      #################################################
      ### filter based on which values are at least 10 
      sumheavy2<-sumheavy[which(sumheavy>=threshold)]
      sumlight2<-sumlight[which(sumheavy>=threshold)]
      
      ### vector of ion names for reporting
      nprec<-length(temppep@prec.z)
      ionnames<-rep(c(paste("b",temppep@iontypes$bions,sep=""),paste("b",temppep@iontypes$bions,sep="","++"),paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,"++",sep="")),times=nprec)
      names(sumheavy)<-sumheavynames
      names(sumlight)<-sumlightnames
      
      #### get and store ion names and masses
      names(sumlight2)<-sumlightnames[which(sumheavy>=threshold)]
      names(sumheavy2)<-sumheavynames[which(sumheavy>=threshold)]
      ordinals<-ionnames[which(sumheavy>=threshold)]
      temppep@iontypes[["position2"]]<-list(ordinals=ordinals,massL=names(sumlight2),massH=names(sumheavy2))
      
      ratios<-sumlight2/(sumlight2+sumheavy2)
      #### return the media value of all ratios    
      temppep@areas[["position2"]]<-list(light=sumlight,heavy=sumheavy,ratio=ratio)
      temppep@areas.filtered[["position2"]]<-list(light=sumlight2,heavy=sumheavy2,ratio=ratios)
      temppep@median.ratio[["position2"]]<-round(median(na.omit(ratios)), digits=4)
      if(length(ratio)>=1){
        temppep@median.ratio.unfilt[["position2"]]<-median(round(na.omit(ratio),digits=4))
      }
      print(k)
      print("double K peptide finished")
      ### unused part to calc stoich from top3 or from linear model
      #if(length(sumheavy)>2){
      #  
      #  y=sumlight
      #  x=sumheavy
      # model<-lm(y~x)
      #  temppep@lm.ratio[["position2"]]<-coef(model)[2]
      #  temppep@top3.ratio[["position2"]]<-ave(sumlight[order(sumheavy,decreasing=T)][1:3]/(sumheavy[order(sumheavy,decreasing=T)][1:3]+sumlight[order(sumheavy,decreasing=T)][1:3]))[1]
      #  
      #}
      temppep@heavy.rank1.ratio[["position2"]]<-list(filtered=sumlight2[sumheavy2==max(sumheavy2)]/(sumheavy2[sumheavy2==max(sumheavy2)]+sumlight2[sumheavy2==max(sumheavy2)]))
      temppep@heavy.rank1.ratio[["position2"]][["unfiltered"]]<-sumlight[sumheavy==max(sumheavy)]/(sumheavy[sumheavy==max(sumheavy)]+sumlight[sumheavy==max(sumheavy)])
      ##### light rank1 ratios
      temppep@light.rank1.ratio[["position2"]]<-list(filtered=sumlight2[sumlight2==max(sumlight2)]/(sumheavy2[sumlight2==max(sumlight2)]+sumlight2[sumlight2==max(sumlight2)]))
      temppep@light.rank1.ratio[["position2"]][["unfiltered"]]<-sumlight[sumlight==max(sumlight)]/(sumheavy[sumlight==max(sumlight)]+sumlight[sumlight==max(sumlight)])
      
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
  peplist  
}
