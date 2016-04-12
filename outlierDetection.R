raw<-read.delim(file="Transitionlvl_BSA_testoutliers.tsv", header=T,sep="\t",stringsAsFactors=F)
raw[1,]
ncol(raw)

### how many samples and datapoints per sample and where does each start
nsamples<-c(5)
npoints<-c(7)
names(raw)
startcol<-15
sstarts<-seq(from=startcol,to=npoints*nsamples+startcol,by=npoints)

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


x=1

onepepmat<-raw[pep.all.lines[[1]],sstarts[1]:(sstarts[1]+npoints-1)]


pep.all.lines[[1]]



#### plot the lines of interest
x<-
y<-
z<-
plot3D(x, y, z)
colors<-rep(1, times=25)
colors[1:2]<-2
colors[11:25]<-3
colors
plot3d(onepepmat[,4:6],col=colors)



