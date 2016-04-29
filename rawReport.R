

k<-paste(object[[6]]@sequence,collapse="")

temppep<-object[[9]]
i=9

testreport<-fragmentRawReport(object=stoich0,output="testout2.csv")

fragmentRawReport=function(object=stoich0,
                           skyline.report=set0pct,
                           output="testout1.csv"){
  peptides<-skyline.report[,1]
  npep<-length(object)
  unshared.win.count=0
  fragreport<-data.frame()
  for(i in 1:npep){
    temppep<-object[[i]]
    #skyreportlines<-which(peptides==paste(temppep@sequence,collapse=""))
    #skyline.report[]
    nprec<-length(temppep@prec.z)
    nmod<-length(temppep@modpos)
    print(i)
    ##########################################################################################
    ##########      single K peptides
    ##########################################################################################

    if(length(temppep@prec.z)<1){
      print("unshared")
      unshared.win.count=unshared.win.count+1
    }
    if(nmod==1 &  length(temppep@prec.z)>=1){
      
      ### LIGHT precursor and frag vectors
      #light.prec.vec<-substr(names(temppep@areas$position1$light),"\\=","\\,")
      #light.prec.vec<-sub(x=names(temppep@areas$position1$light),"^*\\=","")

      light.prec.vec<-sub(x=names(temppep@areas$position1$light),".*prec=(.*),.*","\\1")
      light.frag.vec<-sub(x=names(temppep@areas$position1$light),".*frag=(.*).*$","\\1")
      ### LIGHT charge vec
      light.prec.z.vec<-c()
      for(x in temppep@prec.z){
        light.prec.z.vec<-c(light.prec.z.vec,rep(x,times=length(light.prec.vec)/nprec))
      }
      
      ### HEAVY prec and frag vectors
      heavy.prec.vec<-substr(names(temppep@areas$position1$heavy),start=6,stop=10)
      heavy.frag.vec<-substr(names(temppep@areas$position1$heavy),start=18,stop=25)
      
      ### HEAVY charge vec
      heavy.prec.z.vec<-c()
      for(x in temppep@prec.z){
        heavy.prec.z.vec<-c(heavy.prec.z.vec,rep(x,times=length(heavy.prec.vec)/nprec))
      }
      
      ###  area vectors
      light.area.vec<-temppep@areas$position1$light
      heavy.area.vec<-temppep@areas$position1$heavy
      
      ### part to determine ion rank or just order them by rank
      heavy.ranks<-as.numeric(rank(-heavy.area.vec,ties.method = "first"))
      
      ### combined heavy/light vectors
      precursor.mz<-c(light.prec.vec,heavy.prec.vec)
      precursor.z<-c(light.prec.z.vec,heavy.prec.z.vec)
      fragment.mz<-c(light.frag.vec,heavy.frag.vec)
      area<-c(light.area.vec,heavy.area.vec)
      label.type<-c(rep("L",times=length(light.area.vec)),rep("H",times=length(heavy.prec.vec)))
      observed.rank<-rep(heavy.ranks,times=2)

      
      temp.reportline<-data.frame(Protein=temppep@protein.name,Peptide=paste(temppep@sequence,collapse=""),peptide.mod.pos=temppep@modpos,protein.mod.pos=temppep@protein.position)
      addlines<-(length(light.area.vec)*2)
      temp.reportfront<-data.frame()
      for(j in 1:addlines){
        temp.reportfront<-rbind(temp.reportfront,temp.reportline)
      }
      
      #### add part as below for nonsense with fragments above 1500
      yfrag<-temppep@ionlist$LL$y_diff[na.omit(temppep@ionlist$L$b_diff<1500)]
      iontype<-rep(rep(c(c(paste("b",temppep@iontypes$bions,sep=""),paste("b",temppep@iontypes$bions,sep="","++"))[na.omit(temppep@ionlist$L$b_diff<1500)],c(paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,"++",sep=""))[na.omit(temppep@ionlist$L$y_diff<1500)]),times=nprec),times=2)

      #### build report
      temp.fragreport<-cbind(temp.reportfront,precursor.mz,precursor.z,fragment.mz,iontype,area,observed.rank,label.type)
    }
    
    ##########################################################################################
    ##########      single K peptides
    ##########################################################################################
    
    if(nmod==2 &  length(temppep@prec.z)>=1){
      #pos1.light.prec.vec<-as.numeric(substr(names(temppep@areas$position1$light),start=6,stop=10))
      #pos1.light.frag.vec<-substr(names(temppep@areas$position1$light),start=18,stop=25)
      #pos2.light.prec.vec<-as.numeric(substr(names(temppep@areas$position2$light),start=6,stop=10))
      #pos2.light.frag.vec<-substr(names(temppep@areas$position2$light),start=18,stop=25)
      
      pos1.light.prec.vec<-sub(x=names(temppep@areas$position1$light),".*prec=(.*),.*","\\1")
      pos1.light.frag.vec<-sub(x=names(temppep@areas$position1$light),".*frag=(.*).*$","\\1")
      pos2.light.prec.vec<-sub(x=names(temppep@areas$position2$light),".*prec=(.*),.*","\\1")
      pos2.light.frag.vec<-sub(x=names(temppep@areas$position2$light),".*frag=(.*).*$","\\1")
      
      ### LIGHT charge vec
      
      pos1.light.prec.z.vec<-c()
      pos2.light.prec.z.vec<-c()
      for(x in temppep@prec.z){
        pos1.light.prec.z.vec<-c(pos1.light.prec.z.vec,rep(x,times=length(pos1.light.prec.vec)/nprec))
        pos2.light.prec.z.vec<-c(pos2.light.prec.z.vec,rep(x,times=length(pos2.light.prec.vec)/nprec))
      }
      
      #pos1.heavy.prec.vec<-substr(names(temppep@areas$position1$heavy),start=6,stop=10)
      #pos1.heavy.frag.vec<-substr(names(temppep@areas$position1$heavy),start=18,stop=25)
      #pos2.heavy.prec.vec<-substr(names(temppep@areas$position2$heavy),start=6,stop=10)
      #pos2.heavy.frag.vec<-substr(names(temppep@areas$position2$heavy),start=18,stop=25)
      
      pos1.heavy.prec.vec<-sub(x=names(temppep@areas$position1$heavy),".*prec=(.*),.*","\\1")
      pos1.heavy.frag.vec<-sub(x=names(temppep@areas$position1$heavy),".*frag=(.*).*$","\\1")
      pos2.heavy.prec.vec<-sub(x=names(temppep@areas$position2$heavy),".*prec=(.*),.*","\\1")
      pos2.heavy.frag.vec<-sub(x=names(temppep@areas$position2$heavy),".*frag=(.*).*$","\\1")
      
      ### HEAVY prec charge vectors
      pos1.heavy.prec.z.vec<-c()
      pos2.heavy.prec.z.vec<-c()
      for(x in temppep@prec.z){
        pos1.heavy.prec.z.vec<-c(pos1.heavy.prec.z.vec,rep(x,times=length(pos1.heavy.prec.vec)/nprec))
        pos2.heavy.prec.z.vec<-c(pos2.heavy.prec.z.vec,rep(x,times=length(pos2.heavy.prec.vec)/nprec))
      }
      
      pos1.light.area.vec<-temppep@areas$position1$light
      pos1.heavy.area.vec<-temppep@areas$position1$heavy
      pos2.light.area.vec<-temppep@areas$position2$light
      pos2.heavy.area.vec<-temppep@areas$position2$heavy
      
      
      ### part to determine ion rank or just order them by rank
      pos1.heavy.ranks<-as.numeric(rank(-pos1.heavy.area.vec,ties.method = "first"))
      pos2.heavy.ranks<-as.numeric(rank(-pos2.heavy.area.vec,ties.method = "first"))
      
      ### combined heavy/light vecotrs
      pos1.precursor.mz<-c(pos1.light.prec.vec,pos1.heavy.prec.vec)
      pos1.fragment.mz<-c(pos1.light.frag.vec,pos1.heavy.frag.vec)
      pos1.area<-c(pos1.light.area.vec,pos1.heavy.area.vec)
      pos1.label.type<-c(rep("L",times=length(pos1.light.area.vec)),rep("H",times=length(pos1.heavy.prec.vec)))
      pos1.precursor.z<-c(pos1.light.prec.z.vec,pos1.heavy.prec.z.vec)
      pos1.observed.rank<-rep(pos1.heavy.ranks,times=2)
      
      pos2.precursor.mz<-c(pos2.light.prec.vec,pos2.heavy.prec.vec)
      pos2.fragment.mz<-c(pos2.light.frag.vec,pos2.heavy.frag.vec)
      pos2.area<-c(pos2.light.area.vec,pos2.heavy.area.vec)
      pos2.label.type<-c(rep("L",times=length(pos2.light.area.vec)),rep("H",times=length(pos2.heavy.prec.vec)))
      pos2.precursor.z<-c(pos2.light.prec.z.vec,pos2.heavy.prec.z.vec)
      pos2.observed.rank<-rep(pos2.heavy.ranks,times=2)
      
      
      pos1.temp.reportline<-data.frame(Protein=temppep@protein.name,Peptide=paste(temppep@sequence,collapse=""),peptide.mod.pos=temppep@modpos[1],protein.mod.pos=temppep@protein.position[1])
      pos2.temp.reportline<-data.frame(Protein=temppep@protein.name,Peptide=paste(temppep@sequence,collapse=""),peptide.mod.pos=temppep@modpos[2],protein.mod.pos=temppep@protein.position[2])
      
      
      pos1.addlines<-(length(pos1.light.area.vec)*2)
      pos2.addlines<-(length(pos2.light.area.vec)*2)
      
      pos1.temp.reportfront<-data.frame()
      pos2.temp.reportfront<-data.frame()
      
      for(j in 1:pos1.addlines){
        pos1.temp.reportfront<-rbind(pos1.temp.reportfront,pos1.temp.reportline)
      }
      
      for(j in 1:pos2.addlines){
        pos2.temp.reportfront<-rbind(pos2.temp.reportfront,pos2.temp.reportline)
      }
      
      bfrag<-temppep@ionlist$LL$b_diff[temppep@ionlist$LL$b_diff<1500]
      yfrag<-temppep@ionlist$LL$y_diff[temppep@ionlist$LL$y_diff<1500]
      
      
      pos1.iontype<-rep(rep(c(paste("b",temppep@iontypes$bions,sep=""),paste("b",temppep@iontypes$bions,sep="","++"))[na.omit(temppep@ionlist$LL$b_diff<1500)],times=nprec),times=2)
      pos2.iontype<-rep(rep(c(paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,sep="","++"))[na.omit(temppep@ionlist$LL$y_diff<1500)],times=nprec),times=2)
      
      
      #pos2.iontype<-rep(rep(c(paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,"++",sep="")),times=nprec),times=2)
      
      #### build report
      #cbind(temp.reportfront,precursor.mz,fragment.mz,iontype,area,label.type)
      pos1.temp.fragreportlines<-cbind(pos1.temp.reportfront,precursor.mz=pos1.precursor.mz,precursor.z=pos1.precursor.z,fragment.mz=pos1.fragment.mz,iontype=pos1.iontype,area=pos1.area,observed.rank=pos1.observed.rank,label.type=pos1.label.type)
      pos2.temp.fragreportlines<-cbind(pos2.temp.reportfront,precursor.mz=pos2.precursor.mz,precursor.z=pos2.precursor.z,fragment.mz=pos2.fragment.mz,iontype=pos2.iontype,area=pos2.area,observed.rank=pos2.observed.rank,label.type=pos2.label.type)
      temp.fragreport<-rbind(pos1.temp.fragreportlines,pos2.temp.fragreportlines)
    }
    
  fragreport<-rbind(fragreport,temp.fragreport)
  
  }
  print("proportion with no shared prec window")
  print(unshared.win.count/npep)
  write.csv(file=output,fragreport,row.names = F)
  fragreport
}

