


#### make function to return named list of fragment masses


getIons=function(sequence=peptides[21],modmass=100.016044,moddelta=4.025107){
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
	if(Kcount==2){
		#### determine ordinals for diff ions and those containing both
		Bdiffions<-seq(from=Kpos[1],to=Kpos[2]-1)
		Ydiffions<-seq(from=peplen-Kpos[2]+1,to=peplen-Kpos[1])
		Bbothions<-seq(from=Kpos[2],to=peplen)
		Ybothions<-seq(from=peplen-Kpos[1]+1,to=peplen)

		#### determine ion masses for light/light
		#### ions cannot be ==peplen or 1
		ionlist_ll<-list()
		for(x in Bdiffions){
			if(x>=2 & x!=peplen){
				ionlist_ll[["b_diff"]][[x]]<-round(sum(unlist(AA.masses[seqsplit[Kpos[1]:x]]))+modmass+proton,digits=4)
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

		#ionlist_ll
		#ionlist_hh
		#ionlist_lh
		#ionlist_hl

		}

	#### have lists to extract
	ionlists<-list(LL=ionlist_ll,LH=ionlist_lh,HL=ionlist_hl,HH=ionlist_hh)
	return(ionlist)
	}



getPepPrecMz=function(pepmass=M,chargemax=5){
		listcharge<-c()
	for(i in 1:chargemax){
		listcharge[i] = (pepmass+i*1.00727)/i
		}
	return(listcharge)
	}