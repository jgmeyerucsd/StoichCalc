


getwd()
setwd("C:/Users/jmeyer/Documents/Stoichiometry/")
files<-list.files()
files


####  read in acK custom skyline fragment-level ratio to global standards report
raw<-read.delim(file="TransitionLvl4Stoich_spectralMS2only_succfedfracs.tsv.txt", header=T,sep="\t",stringsAsFactors=F)
raw[1,]
ncol(raw)
raw[1,15]
#### subtract background values from area values

#### 18 samples, subtract 1 from 2, 3 from 4, etc
backgroundvec<-seq(from=15, to=50, by=2)
areavec<-seq(from=16, to=50, by=2)
i=1
bgsub<-list()
substr(colnames(raw),start=22,stop=30)->names
for(i in 1:length(backgroundvec)){
	bgsub[[i]]<-as.numeric(raw[,areavec[i]])-as.numeric(raw[,backgroundvec[i]])
	names(bgsub)[i]<-names[areavec[i]]
	}

raw[1,]
#### which lines are precursors for MS1 only calculation
precursor.lines<-which(raw[,"Fragment.Ion.Type"]=="precursor")
precursor.lines<-which(raw[,"Fragment.Ion.Type"]=="precursor")
fragment.lines<-which(raw[,"Fragment.Ion.Type"]!="precursor")
light.lines<-which(raw[,"Isotope.Label.Type"]=="light")
label.lines<-which(raw[,"Isotope.Label.Type"]!="light")

all.lines<-seq(from=1,to=nrow(raw))

#### which lines belong to which peptide sequence?
peptides<-unique(raw[,"Peptide.Modified.Sequence"])


#### loop through each peptide and make a named list where the peptide name slot contains which lines correspond to it
### empty lists
pep.all.lines<-list()
pep.all.light.lines<-list()
pep.all.heavy.lines<-list()
pep.prec.light.lines<-list()
pep.prec.heavy.lines<-list()
pep.frag.light.lines<-list()
pep.frag.heavy.lines<-list()

for(x in peptides){
	print(x)
	#### fill up all lines matching that peptide
	pep.all.lines[[x]]<-which(raw[,"Peptide.Modified.Sequence"]==x)
	### fill up all light and heavy lines matching that peptide
	pep.all.light.lines[[x]]<-pep.all.lines[[x]][which(raw[pep.all.lines[[x]],"Isotope.Label.Type"]=="light")]
	pep.all.heavy.lines[[x]]<-pep.all.lines[[x]][which(raw[pep.all.lines[[x]],"Isotope.Label.Type"]!="light")]
	#### fill up all light precursor lines
	pep.prec.light.lines[[x]]<-pep.all.light.lines[[x]][which(raw[pep.all.light.lines[[x]],"Fragment.Ion.Type"]=="precursor")]
	pep.prec.heavy.lines[[x]]<-pep.all.heavy.lines[[x]][which(raw[pep.all.heavy.lines[[x]],"Fragment.Ion.Type"]=="precursor")]
	### fill up all light and heavy fragment lines
	pep.frag.light.lines[[x]]<-pep.all.light.lines[[x]][which(raw[pep.all.light.lines[[x]],"Fragment.Ion.Type"]!="precursor")]
	pep.frag.heavy.lines[[x]]<-pep.all.heavy.lines[[x]][which(raw[pep.all.heavy.lines[[x]],"Fragment.Ion.Type"]!="precursor")]
	}


#### now compute MS1-level stoichiometries

stoich.value.list<-list()

for(x in peptides){
	### first compute prec. stoich.
	templight<-bgsub[[1]][pep.prec.light.lines[[x]]]
	templight[templight<0]<-0
	templight.sum<-sum(templight)
	tempheavy<-bgsub[[1]][pep.prec.heavy.lines[[x]]]
	tempheavy.sum<-sum(tempheavy)
	temp.stoich<-sum(templight)/(sum(templight)+sum(tempheavy))
	stoich.value.list[[x]][["prec. stoich"]]<-temp.stoich
	#### now for fragment-level
	templight<-bgsub[[1]][pep.frag.light.lines[[x]]]
	templight[templight<0]<-0
	templight.sum<-sum(templight)
	tempheavy<-bgsub[[1]][pep.frag.heavy.lines[[x]]]
	tempheavy.sum<-sum(tempheavy)
	temp.stoich<-sum(templight)/(sum(templight)+sum(tempheavy))
	stoich.value.list[[x]][["frag. stoich"]]<-temp.stoich
	}

##### have computed above using the unfractionated
raw[,"fraction."]
### make vector of which fraction the peptide is in
fraction.num<-list()

for(x in peptides){
	min(pep.all.lines[[x]])
	fraction.num[[x]]<-raw[min(pep.all.lines[[x]]),"fraction."]
	}

#### do the same as above but put in the fraction of that peptide for the bgsub
names(bgsub)
pos.vec<-c(3,4,5,6,7,8,9,10)

for(x in peptides){
	list.pos<-pos.vec[fraction.num[[x]]]
	### first compute prec. stoich.
	templight<-bgsub[[list.pos]][pep.prec.light.lines[[x]]]
	templight[templight<0]<-0
	templight.sum<-sum(templight)
	tempheavy<-bgsub[[list.pos]][pep.prec.heavy.lines[[x]]]
	tempheavy.sum<-sum(tempheavy)
	temp.stoich<-sum(templight)/(sum(templight)+sum(tempheavy))
	stoich.value.list[[x]][["prec. stoich bRP"]]<-temp.stoich
	#### now for fragment-level
	templight<-bgsub[[list.pos]][pep.frag.light.lines[[x]]]
	templight[templight<0]<-0
	templight.sum<-sum(templight)
	tempheavy<-bgsub[[list.pos]][pep.frag.heavy.lines[[x]]]
	tempheavy.sum<-sum(tempheavy)
	temp.stoich<-sum(templight)/(sum(templight)+sum(tempheavy))
	stoich.value.list[[x]][["frag. stoich bRP"]]<-temp.stoich
	}
	

##### adjust the above code to look for outliers in the data


write.table(stoich.value.list,file="testwrite.txt",sep="\t")






















proteins.kac<-substr(as.character(Kac[,"Protein.Name"]),start=4,stop=9)
unique(proteins.kac)->unique.prot.kac
unique.prot.kac<-unique.prot.kac[nchar(unique.prot.kac)>=2]


####filter to remove iRT peptides
kac.filt<-Kac[which(nchar(proteins.kac)>=2),]
proteins.kac<-substr(as.character(kac.filt[,"Protein.Name"]),start=4,stop=9)
uniq.prot.kac<-unique(proteins.kac)



#### 	fix acK report to have acK site appended to protein name
kac.pos<-mod.pos.protein(table=kac.filt)
#sapply(FUN=paste,X=kac.pos[1:10,2],kac.pos[1:10,1],sep="_")
protname<-paste(kac.pos[,2],kac.pos[,1],sep="_")
kac.table<-cbind(protname,kac.filt[,2],kac.filt[,5],kac.filt[,7:66])
colnames(kac.table)[2]<-"Peptide"
colnames(kac.table)[1]<-"Protein"
colnames(kac.table)[3]<-"FragMz"

### output test table for mapDIA
write.table(kac.table,file="ratios.table.test1.txt",row.names=F,quote=F,sep="\t")

#### map file names to replicate names
filemapping.Kac<-read.delim(file="Kac replicate and file name.csv", header=T,sep=",")
colnames(kac.table)[4:63]<-substr(as.character(filemapping.Kac[,"File.Name"]),start=13,stop=17)

#######


#######  addd a small number to each ratio
kac.table<-kac.table[,order(colnames(kac.table))]
addsmallnum=function(x){
	x=as.numeric(x)+0.0000001
	}
kac<-sapply(FUN=addsmallnum, kac.table[,1:60])
kac.ratio.table<-cbind(kac,kac.table[,61:63])

#### write table of Kac ratios for QC later
write.table(kac.ratio.table,file="20151106.ratios.table.txt",row.names=F,quote=F,sep="\t")

##############################################################################################
##############################################################################################
###		now Kac measurement table is ready with ratios to global standards
##############################################################################################
###		prepare prot level ratios and adjust the site-level ratios
##############################################################################################


files
### read in prot-level quant values from spectronaut
protlvl<-read.delim(file="20150209_ProteinQuant_filtered.txt", header=T,sep="\t",colClasses=c("character",rep("numeric",times=60)))
protlvl[1,]


### fix column names protlvl
colnamestemp<-substr(colnames(protlvl)[2:61],start=24,stop=30)
write.table(colnamestemp,file="protlvl_colnamestemp2.tsv",sep="\t",quote=F)
colnames(protlvl)[2:61]<-as.character(read.table(file="protlvl_colnamestemp.tsv",sep="\t")[,1])

### 	get list of proteins and unique proteins
protlvl.protlist<-substr(protlvl[,1],start=4,stop=9)
protlvl.uniq.prot<-unique(protlvl.protlist)
protlvl.uniq.prot<-protlvl.uniq.prot[nchar(protlvl.uniq.prot)>=2]

#### make prot.table object
prot.table<-matrix(ncol=60,nrow=0)
i=1
protlvl.uniq.prot[1]->x
for(x in protlvl.uniq.prot){
	protlvl[which(protlvl.protlist==x),]
	colSums(x=protlvl[which(protlvl.protlist==x),2:61])
	prot.table<-rbind(prot.table,colSums(x=protlvl[which(protlvl.protlist==x),2:61]))
	rownames(prot.table)[i]<-as.character(x)
	i=i+1
	}

#### 	compute ratios to sets 2w chow water control animal
prot.table<-prot.table[,order(colnames(prot.table))]
prot.table[1,]
protlvl.protlist[1]
nrow(prot.table)
protratio.table<-matrix(ncol=0,nrow=nrow(prot.table))
for(i in 1:5){
	setrowstemp<-seq(from=i+5,to=60,by=5)
	tempratios<-prot.table[,setrowstemp]/prot.table[,i]
	protratio.table<-cbind(tempratios,protratio.table)
	}
protratio.table[1,]
ncol(protratio.table)
protratio.table<-protratio.table[,order(colnames(protratio.table))]



#### how many proteins with acK sites are quantified?
na.omit(match(rownames(prot.table),unique.prot.kac))

protratio.table1[1,]


########################################################################
###	write a table with protein-level ratios for any followup needed
write.table(protratio.table,file="20151106.protlvlonly.ratios.txt",row.names=T,quote=F,sep="\t")

######################################################################
####	protein ratios ready
######################################################################
#### look for acetyl proteins in protein ratio table
#### divide the fragment ratios by protein ratios
#### then make a new matrix for mapDIA
######################################################################

kac[1,]
kac.ratio.table[1,]
kac.ratio.table<-kac
kac.ratio.table<-kac.ratio.table[,order(colnames(kac.ratio.table))]
kac.ratio.table[1,]
kac.prot<-substr(kac.ratio.table[,"Protein"],start=4,stop=9)
unique.prot.kac<-unique(substr(kac.ratio.table[,"Protein"],start=4,stop=9))

protratio.table[1,]
unique.prot.kac<-unique(proteins.kac.ratio.table)

#### add five columns to protlvl ratio table with ratio ==1 for controls
protratio.table[1,]
nrow(protratio.table)
for(x in 1:5){
	protratio.table<-cbind(rep(1,times=777),protratio.table)
	}
colnames(protratio.table)[1:5]<-c("02w01","02w02","02w03","02w04", "02w05")

####  TAKE RATIOS
adjratios.table<-matrix(ncol=63,nrow=0)

x<-unique.prot.kac[1]

for(x in unique.prot.kac){
	protlvl.temprow<-which(rownames(protratio.table)==x)
	print(length(protlvl.temprow))
	if(length(protlvl.temprow)>0){
		kac.temprows<-which(proteins.kac==x)
		#protratio.table[protlvl.temprow,]

		adjratios.table<-rbind(adjratios.table,cbind(t(t(kac.ratio.table[kac.temprows,1:60])/protratio.table[protlvl.temprow,]),kac.ratio.table[kac.temprows,61:63]))
		}
	}

	mat<-matrix(1,ncol=2,nrow=2,TRUE)
	dev<-c(5,10)
	t(t(mat) / dev)
	kac.ratio.table[1:2,]
	t(kac.ratio.table[kac.temprows,1:60])
	protratio.table[protlvl.temprow,]
	t(t(kac.ratio.table[kac.temprows,1:60])/protratio.table[protlvl.temprow,])
	kac.ratio.table[1:2,1:60]/protratio.table[protlvl.temprow,]
adjratios.table[1,]
kac.ratio.table[1,]
protratio.table[protlvl.temprow,]
ncol(kac.ratio.table[kac.temprows,1:60])
ncol(protratio.table[protlvl.temprow,])

####
kac[1,1:5]
nrow(kac)
nrow(adjratios.table)
nrow(kac.ratio.table)
which(

kac.ratio.table[1,]
adjratios.table[1,]

adjratios.table[100,]
min(na.omit(adjratios.table[,1:60]))*3000000
#####  filter the mapDIA output table to only include those proteins adjusted for quantity

final.table<-adjratios.table

#unique(substr(as.character(adjratios.table[,"Protein"]),start=4,stop=9))

columns2w<-matrix(ncol=0,nrow=nrow(adjratios.table))
for(i in 1:5){
	columns2w<-cbind(columns2w,rep(1,times=nrow(adjratios.table)))
	}




final.table<-cbind(final.table[,63],final.table[,62],final.table[,61],final.table[,1:60]*3000000)
colnames(final.table)

colnames(final.table)<-c("protein","Peptide","FragMz",as.character(filemapping.Kac[,"Replicate"]))
final.table[1,]
typeof(as.character(filemapping.Kac[,"Replicate"]))






write.table(final.table,file="20151106.Kac.protlvladj.60files.txt",row.names=F,quote=F,sep="\t")
?apply


pctmissing<-apply(final.table[,4:63], 1, function(x) sum(is.na(x))) / ncol(final.table[,4:63]) * 100
numNAs <- apply(final.table[,4:63], 1, function(x) sum(is.na(x)))

which(pctmissing==100)

final.table

##### remove rows with NA
##### did this in excel temporarily





countNA<-
na.omit(


apply(FUN=na.omit,1,X=final.table)

min(adjratios.table[,1:55])*10000000





#####################   read in mapDIA output and make volcano plots and heatmaps
mapdia<-read.delim(file="C:/mapDIA_v2.0.3/windows64binary/final_acK/peptide_level.txt", header=T,sep="\t",stringsAsFactors=F)

###	average the individuals
averaged.table<-cbind(mapdia[,"Protein"],rowSums(mapdia[,3:7])/5,rowSums(mapdia[,8:12])/5,
	rowSums(mapdia[,13:17])/5,rowSums(mapdia[,18:22])/5,
	rowSums(mapdia[,23:27])/5,rowSums(mapdia[,28:32])/5,
	rowSums(mapdia[,33:37])/5,rowSums(mapdia[,38:42])/5,
	rowSums(mapdia[,43:47])/5,rowSums(mapdia[,48:52])/5,
	rowSums(mapdia[,53:57])/5,rowSums(mapdia[,58:62])/5)


typeof(averaged.table[1:2,])
as.numeric(averaged.table[1:2,2:13][1,])+as.numeric(averaged.table[1:2,2:13][2,])
	


as.data.frame(averaged.table,stringsAsFactors=F,colClasses=c("character",rep("numeric",times=12)))->mapdia.grouped
typeof(mapdia.grouped[which(sites==x),2:13])
as.numeric(mapdia.grouped[which(sites==x),2:13])

matrix(as.numeric(unlist(mapdia.grouped[,2:13])),nr=nrow(mapdia.grouped[,2:13]),dimnames=list(mapdia.grouped[,1]))->value.matrix


mapdia[1,]
uniq.sites<-unique(mapdia[,"Protein"])
sites<-mapdia[,"Protein"]
newmat<-matrix(nc=12,nr=0)
#### sum of all site signals 
uniq.sites[1]->x
for(x in uniq.sites){
	#cbind(mapdia[min(which(sites==x)),1],mapdia.grouped[which(sites==x),2:13]
	#mapdia.grouped[which(sites==x),2:13]
	if(length(which(sites==x))>1){
		newmat<-rbind(newmat,colSums(value.matrix[which(sites==x),1:12]))
		}
	if(length(which(sites==x))==1){
		newmat<-rbind(newmat,value.matrix[which(sites==x),1:12])
		}
	
	}

dimnames(newmat)<-list(uniq.sites)
heatmap(newmat,Colv=NA,labRow=NA)


################## volcano plots
mapdia<-read.delim(file="C:/mapDIA_v2.0.3/windows64binary/final_acK/peptide_level.txt", header=T,sep="\t",stringsAsFactors=F)
output<-read.delim(file="C:/mapDIA_v2.0.3/windows64binary/final_acK/analysis_output.txt", header=T,sep="\t",stringsAsFactors=F)


output$comparisons<-factor(output[,"Label2"])

output[1,]
?sample
sample(x=c("TRUE","FALSE"),10,replace=TRUE)

#### filter rows to contain only those with FDR<0.01
out.filt1<-output[which(output[,"FDR"]<=0.01),]
nrow(output)
nrow(out.filt1)
out.filt2<-output[which(output[,"log2FC"]<=-1 | output[,"log2FC"]>=1),]
out.filt2[,"log2FC"]

####	plot all significant changes
plot(x=-output[,"log2FC"],y=output[,"log_oddsDE"],pch=20,main="all significant changes",xlab="log2(Fold Change)", ylab="odds differential expression")
points(x=-out.filt2[,"log2FC"],y=out.filt2[,"log_oddsDE"],pch=20,col="red")


####	generate volcano plots of all individual comparisons


volcano.plots=function(output=output){
	output$comparisons<-factor(output[,"Label2"])
	#### subsets filtered by FDR and log fold change
	out.filt<-output[which(output[,"FDR"]<=0.01),]
	out.filt<-output[which(out.filt[,"log2FC"]<=-1 | out.filt[,"log2FC"]>=1),]
	
	num.compare<-length(levels(output$comparisons))
	

	##### loop through the comparison sets and plot

	for(x in levels(output$comparisons)){


		tiff(filename = "grouped.adjusted.tif",
			width = 2000, height = 2000, units = "px", pointsize = 12,
			compression = c("none"),bg = "white", res = 300, restoreConsole = TRUE,
			type = c("windows"))

		plot(x=-output[,"log2FC"],y=output[,"log_oddsDE"],pch=20,main="all significant changes",xlab="log2(Fold Change)", ylab="odds differential expression")
		points(x=-out.filt2[,"log2FC"],y=out.filt2[,"log_oddsDE"],pch=20,col="red")
		
		
		}
	
	}

levels(output$comparisons)


output[which(output$comparisons==levels(output$comparisons)[4]),]






tiff(filename = "grouped.adjusted.tif",
     width = 2000, height = 2000, units = "px", pointsize = 12,
     compression = c("none"),bg = "white", res = 300, restoreConsole = TRUE,
     type = c("windows"))
heatmap(x=as.matrix(averaged.table[,2:13]),Colv=NA,labRow=NA)
dev.off()

na.omit(match(unique.prot.kac,allprot.list))
mapDIA<-read.delim(file=files[120], header=T,sep="\t")
mapDIA[1,]
