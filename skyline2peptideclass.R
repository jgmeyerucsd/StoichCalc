#setwd("~/stoichiometry")
### read in skyline report
#s<-read.delim(file="stoichout1test7.tsv", head=T,stringsAsFactors=F)


skyline2peptide=function(filename="stoichout1test8.tsv"){
  s<-read.delim(file=filename, head=T,stringsAsFactors=F)
  sreport<-new(Class="skylineReport")
  sreport@peptides<-list()
  pepvec<-unique(s[,"Peptide"])
  for(k in pepvec){
    #k<-pepvec[2]
    print(k)
    pepcount=pepcount+1
    print(pepcount)
    temppep<-new(Class="Peptide")
    temppep@ionlist<-getIons(sequence=k)
    temppep@sequence<-unlist(strsplit(k,split=""))
    temppep@modpos<-which(temppep@sequence=="K")
    temppep@Kcount<-length(temppep@modpos)
    templines<-which(s[,"Peptide"]==k)
    temppep@prec.z<-unique(s[templines,"Precursor.Charge"])
    for(x in temppep@prec.z){
      temppep@reportlines[[x]]<-templines[which(s[templines,"Precursor.Charge"]==x)]
    }
    temppep@peakbounds<-c(sky.report[templines[1],"Start.Time"],sky.report[templines[1],"End.Time"])
    temppep@protein.name<-sky.report[templines[1],"Protein"]
    temppep@protein.position<-sky.report[templines[1],"Begin.Pos"]+temppep@modpos-1
    
    ### determine which ion is rank 1 for each precursor charge
    #### currently just takes the first line from the skyline report that has rank 1 for that precursor charge
    
    for(x in temppep@prec.z){
      temppep@reportlines[[x]][which(as.numeric(s[temppep@reportlines[[x]],"Library.Rank"])==1)]
      bestrank<-min(na.omit(as.numeric(s[temppep@reportlines[[x]],"Library.Rank"])))
      temppep@rank1ion[[x]]<-s[temppep@reportlines[[x]][which(as.numeric(s[temppep@reportlines[[x]],"Library.Rank"])==bestrank)][1],c("Precursor.Charge","Product.Mz","Fragment.Ion","Product.Charge")]
    }
    sreport@peptides[[paste(k)]]<-temppep
  }
  
}

