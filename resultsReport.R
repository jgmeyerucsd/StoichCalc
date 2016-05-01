

k<-paste(object[[6]]@sequence,collapse="")

temppep<-object[[9]]
i=9

results0<-fragmentResultsReport(object=stoich0,output="0pctv2.csv")
results1<-fragmentResultsReport(object=stoich1,output="1pctv2.csv")
results10<-fragmentResultsReport(object=stoich10,output="10pctv2.csv")
results50<-fragmentResultsReport(object=stoich50,output="50pctv2.csv")
results100<-fragmentResultsReport(object=stoich100,output="100pctv2.csv")



###   results report should have:
###   protein
### peptide
### site#
### ion e.g. y7
### L/L+H value
### diff vs. common ion (always diff first)
### observed rank (NA in median calc)
### fraction # eventually


fragmentResultsReport=function(object=stoich0,
                           skyline.report=set0pct,
                           output="testout1.csv"){
  
  peptides<-skyline.report[,1]
  npep<-length(object)
  unshared.win.count=0
  report<-data.frame()
  for(i in 1:npep){
    temppep<-object[[i]]

    nmod<-length(temppep@modpos)
    print(i)
    ##########################################################################################
    ##########      single K peptides
    ##########################################################################################
    
    if(length(temppep@prec.z)<1){
      print("unshared")
      unshared.win.count=unshared.win.count+1
    }

    if(length(temppep@prec.z)>=1){
      
      templine<-data.frame()
      for(x in 1:nmod){
        print(x)
        protein<-temppep@protein.name
        peptide<-paste(temppep@sequence,sep="",collapse="")
        protein.site<-temppep@protein.position[x]
        heavy.rank1.stoich<-as.numeric(temppep@heavy.rank1.ratio[[x]]$filtered)
        heavy.rank1.stoich.unfiltered<-as.numeric(temppep@heavy.rank1.ratio[[x]]$unfiltered)
        light.rank1.stoich<-as.numeric(temppep@light.rank1.ratio[[x]]$filtered)
        light.rank1.stoich.unfiltered<-as.numeric(temppep@light.rank1.ratio[[x]]$unfiltered)
        if(length(heavy.rank1.stoich)==0){
          heavy.rank1.stoich<-NA
        }
        if(length(light.rank1.stoich)==0){
          light.rank1.stoich<-NA
        }
        median.stoich<-temppep@median.ratio[[x]]
        median.stoich.unfiltered<-temppep@median.ratio.unfilt[[x]]
        if(length( median.stoich.unfiltered)==0){
          median.stoich.unfiltered<-NA
        }
        if(x==1){
          templine<-data.frame(protein,peptide,protein.site,median.stoich,median.stoich.unfiltered,heavy.rank1.stoich,heavy.rank1.stoich.unfiltered,light.rank1.stoich,light.rank1.stoich.unfiltered)
        }
        if(x==2){
          templine<-rbind(templine,data.frame(protein,peptide,protein.site,median.stoich,median.stoich.unfiltered,heavy.rank1.stoich,heavy.rank1.stoich.unfiltered,light.rank1.stoich,light.rank1.stoich.unfiltered))
        }
      }
      report<-rbind(report,templine)
    }
  }
  
  names(report)
  
  write.csv(file=paste("summaryReport",output,sep="_",collapse=""),report,row.names = F)
  report
}

boxplot(
  na.omit(results0[,"median.ratio"]),
  na.omit(results1[,"median.ratio"]),
  na.omit(results10[,"median.ratio"]),
  na.omit(results50[,"median.ratio"]),
  na.omit(results100[,"median.ratio"])
  
)
boxplot(
  na.omit(results0[,"median.ratio.unfiltered"]),
  na.omit(results1[,"median.ratio.unfiltered"]),
  na.omit(results10[,"median.ratio.unfiltered"]),
  na.omit(results50[,"median.ratio.unfiltered"]),
  na.omit(results100[,"median.ratio.unfiltered"])
  
)
boxplot(
  na.omit(results0[,"light.rank1.stoich.unfiltered"]),
  na.omit(results1[,"light.rank1.stoich.unfiltered"]),
  na.omit(results10[,"light.rank1.stoich.unfiltered"]),
  na.omit(results50[,"light.rank1.stoich.unfiltered"]),
  na.omit(results100[,"light.rank1.stoich.unfiltered"])
)
### unfiltered heavy rank1
boxplot(
  na.omit(results0[,"heavy.rank1.stoich.unfiltered"]),
  na.omit(results1[,"heavy.rank1.stoich.unfiltered"]),
  na.omit(results10[,"heavy.rank1.stoich.unfiltered"]),
  na.omit(results50[,"heavy.rank1.stoich.unfiltered"]),
  na.omit(results100[,"heavy.rank1.stoich.unfiltered"])
)

boxplot(
  na.omit(results0[,"rank1.ratio.unfiltered"]),
  na.omit(results1[,"rank1.ratio.unfiltered"]),
  na.omit(results10[,"rank1.ratio.unfiltered"]),
  na.omit(results50[,"rank1.ratio.unfiltered"]),
  na.omit(results100[,"rank1.ratio.unfiltered"])
)

boxplot(
  na.omit(results0[,"median.ratio"]),
  na.omit(results0[,"rank1.ratio"]),
  na.omit(results1[,"median.ratio"]),
  na.omit(results1[,"rank1.ratio"]),
  na.omit(results10[,"median.ratio"]),
  na.omit(results10[,"rank1.ratio"]),
  na.omit(results50[,"median.ratio"]),
  na.omit(results50[,"rank1.ratio"])
)

boxplot(
  na.omit(results0[,"median.ratio"]),
  na.omit(results0[,"rank1.ratio"]),
  na.omit(results0[,"rank1.ratio.unfiltered"]),
  na.omit(results1[,"median.ratio"]),
  na.omit(results1[,"rank1.ratio"]),
  na.omit(results1[,"rank1.ratio.unfiltered"]),
  na.omit(results10[,"median.ratio"]),
  na.omit(results10[,"rank1.ratio"]),
  na.omit(results10[,"rank1.ratio.unfiltered"]),
  na.omit(results50[,"median.ratio"]),
  na.omit(results50[,"rank1.ratio"]),
  na.omit(results50[,"rank1.ratio.unfiltered"])
)
abline(h=0.5)
abline(h=0.1)
abline(h=0.01)
abline(h=0.0)
