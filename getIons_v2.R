#### function to return named list of fragment masses given mod mass and heavy/light delta mass
#### returns +1 and +2 charge fragments

getIons=function(sequence=peptides[16],modmass=100.016044,moddelta=4.025107){
	options(digits=10)
	### number and position of lysine residues
	seqsplit<-unlist(strsplit(sequence,split=""))
	seqrev<-rev(seqsplit)
	peplen<-length(seqsplit)
	Kpos<-which(seqsplit=="K")
	Kcount<-length(Kpos)	
	proton=1.00727647
	water=18.010565
	###
	### handle peptides with diff K numbers differently
	if(Kcount==1){
	  Bdiffions<-seq(from=Kpos[1],to=peplen-1)
	  Ydiffions<-seq(from=peplen-Kpos[1]+1,to=peplen-1)
    ionlist_l<-list()
    ionnamelist<-list()
    #for()
    
    for(x in Bdiffions){
      if(x>=2 & x!=peplen){
        ionlist_l[["b_diff"]][[paste("b",x,sep="")]]<-round(sum(unlist(AA.masses[seqsplit[1:x]]))+modmass+proton,digits=4)
        ionnamelist[["b_diff"]]
      }
    }
    for(x in Ydiffions){
      if(x>=2 & x!=peplen){
        ionlist_l[["y_diff"]][[paste("y",x,sep="")]]<-round(sum(unlist(AA.masses[seqrev[1:x]]))+modmass+proton+water,digits=4)
      }
    }
    ### determine ion masses for heavy/heavy
    ionlist_h<-lapply(FUN=function(x) x+moddelta,ionlist_l[1:2])
    
    ### add doubly charged
    for(x in names(ionlist_l)){
      ionlist_l[[x]]<-c(ionlist_l[[x]],(ionlist_l[[x]]+proton)/2)
      names(ionlist_l[[x]])[(length(ionlist_l[[x]])/2+1):length(ionlist_l[[x]])]<-paste("2",names(ionlist_l[[x]])[(length(ionlist_l[[x]])/2+1):length(ionlist_l[[x]])],sep="")
      names(ionlist_l[[x]])[1:length(ionlist_l[[x]])/2]<-paste("1",names(ionlist_l[[x]])[1:length(ionlist_l[[x]])/2],sep="")
      ionlist_h[[x]]<-c(ionlist_h[[x]],(ionlist_h[[x]]+proton)/2)
      names(ionlist_h[[x]])[(length(ionlist_h[[x]])/2+1):length(ionlist_h[[x]])]<-paste("2",names(ionlist_h[[x]])[(length(ionlist_h[[x]])/2+1):length(ionlist_h[[x]])],sep="")
      names(ionlist_h[[x]])[1:length(ionlist_h[[x]])/2]<-paste("1",names(ionlist_h[[x]])[1:length(ionlist_h[[x]])/2],sep="")
    }
    
    ### determine precursor for light
    prec<-list()
     for(i in 1:5){
      prec$L[[i]]<-round(sum(unlist(AA.masses[seqsplit]))+modmass*Kcount+proton*i+water,digits=4)/i
    }
    ### determine heavy precursor masses
    for(i in 1:5){
      prec$H[[i]]<-round(sum(unlist(AA.masses[seqsplit]))+moddelta+modmass*Kcount+proton*i+water,digits=4)/i
    }
    
    ionlist<-list(prec=prec,L=ionlist_l,H=ionlist_h)
    
	}
	if(Kcount==2){
		#### determine ordinals for diff ions and those containing both
		Bdiffions<-seq(from=Kpos[1],to=Kpos[2]-1)
		Ydiffions<-seq(from=peplen-Kpos[2]+1,to=peplen-Kpos[1])
		Bbothions<-seq(from=Kpos[2],to=peplen)
		Ybothions<-seq(from=peplen-Kpos[1]+1,to=peplen)

		#### determine ion masses for light/light
		#### ions cannot be ==peplen or 1
		ionlist_ll<-list()
		iontype_ll<-list()
		for(x in Bdiffions){
			if(x>=2 & x!=peplen){
				ionlist_ll[["b_diff"]][[x]]<-round(sum(unlist(AA.masses[seqsplit[1:x]]))+modmass+proton,digits=4)
				
				}
			}
		for(x in Ydiffions){
			if(x>=2 & x!=peplen){
				ionlist_ll[["y_diff"]][[x]]<-round(sum(unlist(AA.masses[seqrev[1:x]]))+modmass+proton+water,digits=4)
				}
			}
		for(x in Bbothions){
			if(x>=2 & x!=peplen){
				ionlist_ll[["b_both"]][[x]]<-round(sum(unlist(AA.masses[seqsplit[Kpos[1]:x]]))+modmass+proton,digits=4)
				}
			}
		for(x in Ybothions){
			if(x>=2 & x!=peplen){
				ionlist_ll[["y_both"]][[x]]<-round(sum(unlist(AA.masses[seqrev[1:x]]))+modmass+proton+water,digits=4)
				}
			}
		
		### determine ion masses for heavy/heavy
		ionlist_hh<-c(lapply(FUN=function(x) x+moddelta,ionlist_ll[1:2]),lapply(FUN=function(x) x+moddelta*2,ionlist_ll[3:4]))
				
		### determine ion masses for light/heavy
		ionlist_lh<-c(ionlist_ll[1],lapply(FUN=function(x) x+moddelta,ionlist_ll[2]),lapply(FUN=function(x) x+moddelta,ionlist_ll[3:4]))
		
		### determine ion masses for light/heavy
		ionlist_hl<-c(lapply(FUN=function(x) x+moddelta,ionlist_ll[1]),ionlist_ll[2],lapply(FUN=function(x) x+moddelta,ionlist_ll[3:4]))
		
		### add +2 charge fragments to each slot
		#### and fix names to be [precZ:b/y:ordinal]
		for(x in names(ionlist_ll)){
		  print(ionlist_ll[[x]])
		  ionlist_ll[[x]]<-c(ionlist_ll[[x]],(ionlist_ll[[x]]+proton)/2)
		  ionlist_hh[[x]]<-c(ionlist_hh[[x]],(ionlist_hh[[x]]+proton)/2)
		}
		
		### finally, determine precursor masses TODO: clean this up into a loop
		### light/light
		p1ll<-round(sum(unlist(AA.masses[seqsplit]))+modmass*Kcount+proton*1+water,digits=4)/1
		p2ll<-round(sum(unlist(AA.masses[seqsplit]))+modmass*Kcount+proton*2+water,digits=4)/2
		p3ll<-round(sum(unlist(AA.masses[seqsplit]))+modmass*Kcount+proton*3+water,digits=4)/3
		p4ll<-round(sum(unlist(AA.masses[seqsplit]))+modmass*Kcount+proton*4+water,digits=4)/4
		p5ll<-round(sum(unlist(AA.masses[seqsplit]))+modmass*Kcount+proton*5+water,digits=4)/5
		### light/heavy
		p1lh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta+modmass*Kcount+proton*1+water,digits=4)/1
		p2lh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta+modmass*Kcount+proton*2+water,digits=4)/2
		p3lh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta+modmass*Kcount+proton*3+water,digits=4)/3
		p4lh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta+modmass*Kcount+proton*4+water,digits=4)/4
		p5lh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta+modmass*Kcount+proton*5+water,digits=4)/5
		### heavy/heavy
		p1hh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta*2+modmass*Kcount+proton*1+water,digits=4)/1
		p2hh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta*2+modmass*Kcount+proton*2+water,digits=4)/2
		p3hh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta*2+modmass*Kcount+proton*3+water,digits=4)/3
		p4hh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta*2+modmass*Kcount+proton*4+water,digits=4)/4
		p5hh<-round(sum(unlist(AA.masses[seqsplit]))+moddelta*2+modmass*Kcount+proton*5+water,digits=4)/5
		#### list of all precursors and fragments
		ionlist<-list(prec=list(LL=c(p1ll,p2ll,p3ll,p4ll,p5ll),LH=c(p1lh,p2lh,p3lh,p4lh,p5lh),HH=c(p1hh,p2hh,p3hh,p4hh,p5hh)),LL=ionlist_ll,LH=ionlist_lh,HL=ionlist_hl,HH=ionlist_hh)
  	}
	if(Kcount>=3){
	  ionlist<-list() 
	}
	#### have lists to extract
	return(ionlist)
	}
