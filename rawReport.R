

k<-paste(object[[6]]@sequence,collapse="")

temppep<-object[[5]]
i=5

testreport<-fragmentRawReport()

fragmentRawReport=function(object=stoich0,skyline.report=set0pct){
  peptides<-skyline.report[,1]
  unique(skyline.report[,1])
  npep<-length(object)
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
    
    if(nmod==1){
      
      light.prec.vec<-as.numeric(substr(names(temppep@areas$position1$light),start=6,stop=10))
      light.frag.vec<-substr(names(temppep@areas$position1$light),start=18,stop=25)
      
      heavy.prec.vec<-substr(names(temppep@areas$position1$heavy),start=6,stop=10)
      heavy.frag.vec<-substr(names(temppep@areas$position1$heavy),start=18,stop=25)
      
      light.area.vec<-temppep@areas$position1$light
      heavy.area.vec<-temppep@areas$position1$heavy
      
      precursor.mz<-c(light.prec.vec,heavy.prec.vec)
      fragment.mz<-c(light.frag.vec,heavy.frag.vec)
      area<-c(light.area.vec,heavy.area.vec)
      label.type<-c(rep("L",times=length(light.area.vec)),rep("H",times=length(heavy.prec.vec)))
      
      
      
      temp.reportline<-data.frame(Protein=temppep@protein.name,Peptide=paste(temppep@sequence,collapse=""),peptide.mod.pos=temppep@modpos,protein.mod.pos=temppep@protein.position)
      addlines<-(length(light.area.vec)*2)
      temp.reportfront<-data.frame()
      for(j in 1:addlines){
        temp.reportfront<-rbind(temp.reportfront,temp.reportline)
      }
      
      #### add part as below for nonsense with fragments above 1500
      
      iontype<-rep(rep(c(paste("b",temppep@iontypes$bions,sep=""),paste("b",temppep@iontypes$bions,sep="","++"),paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,"++",sep="")),times=nprec),times=2)

      #### build report
      temp.fragreportlines<-cbind(temp.reportfront,precursor.mz,fragment.mz,iontype,area,label.type)
    }
    
    ##########################################################################################
    ##########      single K peptides
    ##########################################################################################
    
    if(nmod==2){
      pos1.light.prec.vec<-as.numeric(substr(names(temppep@areas$position1$light),start=6,stop=10))
      pos1.light.frag.vec<-substr(names(temppep@areas$position1$light),start=18,stop=25)
      pos2.light.prec.vec<-as.numeric(substr(names(temppep@areas$position2$light),start=6,stop=10))
      pos2.light.frag.vec<-substr(names(temppep@areas$position2$light),start=18,stop=25)
      
      pos1.heavy.prec.vec<-substr(names(temppep@areas$position1$heavy),start=6,stop=10)
      pos1.heavy.frag.vec<-substr(names(temppep@areas$position1$heavy),start=18,stop=25)
      pos2.heavy.prec.vec<-substr(names(temppep@areas$position2$heavy),start=6,stop=10)
      pos2.heavy.frag.vec<-substr(names(temppep@areas$position2$heavy),start=18,stop=25)
      
      pos1.light.area.vec<-temppep@areas$position1$light
      pos1.heavy.area.vec<-temppep@areas$position1$heavy
      pos2.light.area.vec<-temppep@areas$position2$light
      pos2.heavy.area.vec<-temppep@areas$position2$heavy
      
      pos1.precursor.mz<-c(pos1.light.prec.vec,pos1.heavy.prec.vec)
      pos1.fragment.mz<-c(pos1.light.frag.vec,pos1.heavy.frag.vec)
      pos1.area<-c(pos1.light.area.vec,pos1.heavy.area.vec)
      pos1.label.type<-c(rep("L",times=length(pos1.light.area.vec)),rep("H",times=length(pos1.heavy.prec.vec)))
      
      pos2.precursor.mz<-c(pos2.light.prec.vec,pos2.heavy.prec.vec)
      pos2.fragment.mz<-c(pos2.light.frag.vec,pos2.heavy.frag.vec)
      pos2.area<-c(pos2.light.area.vec,pos2.heavy.area.vec)
      pos2.label.type<-c(rep("L",times=length(pos2.light.area.vec)),rep("H",times=length(pos2.heavy.prec.vec)))
      
      
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
      
      length(na.omit(temppep@ionlist$LL$b_diff<1500))
      
      pos1.iontype<-rep(rep(c(paste("b",temppep@iontypes$bions,sep=""),paste("b",temppep@iontypes$bions,sep="","++"))[na.omit(temppep@ionlist$LL$b_diff<1500)],times=nprec),times=2)
      pos2.iontype<-rep(rep(c(paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,sep="","++"))[na.omit(temppep@ionlist$LL$y_diff<1500)],times=nprec),times=2)
      
      
      pos2.iontype<-rep(rep(c(paste("y",temppep@iontypes$yions,sep=""),paste("y",temppep@iontypes$yions,"++",sep="")),times=nprec),times=2)
      
      #### build report
      cbind(temp.reportfront,precursor.mz,fragment.mz,iontype,area,label.type)
      pos1.temp.fragreportlines<-cbind(pos1.temp.reportfront,precursor.mz=pos1.precursor.mz,fragment.mz=pos1.fragment.mz,iontype=pos1.iontype,area=pos1.area,label.type=pos1.label.type)
      
      pos2.temp.fragreportlines<-cbind(pos2.temp.reportfront,precursor.mz=pos2.precursor.mz,fragment.mz=pos2.fragment.mz,iontype=pos2.iontype,area=pos2.area,label.type=pos2.label.type)
      temp.fragreport<-rbind(pos1.temp.fragreportlines,pos2.temp.fragreportlines)
      
    }
  fragreport<-rbind(fragreport,temp.fragreport)
    

  }
  fragreport
  
}

