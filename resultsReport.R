

k<-paste(object[[6]]@sequence,collapse="")

temppep<-object[[9]]
i=9

results0<-fragmentResultsReport(object=stoich0,output="0pct.csv")
results1<-fragmentResultsReport(object=stoich1,output="1pct.csv")
results10<-fragmentResultsReport(object=stoich10,output="10pct.csv")
results50<-fragmentResultsReport(object=stoich50,output="50pct.csv")
results100<-fragmentResultsReport(object=stoich100,output="100pct.csv")



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
        rank1.ratio<-as.numeric(temppep@rank1.ratio[[x]]$filtered)
        rank1.ratio.unfiltered<-as.numeric(temppep@rank1.ratio[[x]]$unfiltered)
        if(length(rank1.ratio)==0){
          rank1.ratio<-NA
        }
        median.ratio<-temppep@median.ratio[[x]]
        if(x==1){
          templine<-data.frame(protein,peptide,protein.site,median.ratio,rank1.ratio,rank1.ratio.unfiltered)
        }
        if(x==2){
          templine<-rbind(templine,data.frame(protein,peptide,protein.site,median.ratio,rank1.ratio,rank1.ratio.unfiltered))
        }
      }
      report<-rbind(report,templine)
    }
  }
  
  write.csv(file=paste("summaryReport",output,sep="_",collapse=""),report,row.names = F)
  report
}

boxplot(
  na.omit(results0[,"median.ratio"]),
  na.omit(results1[,"median.ratio"]),
  na.omit(results10[,"median.ratio"]),
  na.omit(results50[,"median.ratio"])
)

boxplot(
  na.omit(results0[,"rank1.ratio"]),
  na.omit(results1[,"rank1.ratio"]),
  na.omit(results10[,"rank1.ratio"]),
  na.omit(results50[,"rank1.ratio"])
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
