
MS1specplot<-function(dbkey,base.dir,finlist,SubDB,TMP=NULL,lc.err=NULL,
				mz.err=NULL,single=NULL,prcx=NULL){
	if (is.null(lc.err)) lc.err=0.02
	if (is.null(TMP)) TMP<-SelectSubDB(base.dir,finlist,SubDB)
	if (is.null(mz.err)){
	 if (TMP[[2]]=="FT"){
	  mz.err=0.001
	 } else if (TMP[[2]]=="QTOF"){
	  mz.err=0.02
	 } else {
	  print("Unknown MS Analyzer, set SubDB equal to FTneg,
					 FTpos, QTOFneg or QTOFpos")
	 }
	}
	if(is.null(single)){
	 oldpar<-par(no.readonly=T)
	 par(mfrow=c(1,1))
	}
	if(is.null(prcx))prcx=0.7
	subDB<-TMP[[1]]
	dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	exp<-subDB[[1]][dbkey_row,4]
	procfc<-LoadProc(exp,TMP)
	proc<-procfc[[1]]
	if(class(proc)=="character"){
	 print("No nodes.txt file present")
	 return()
	}
	fc<-procfc[[2]]
	dMS1<-MS1plot(dbkey,exp,subDB,proc,fc,finlist,lc.err,mz.err,prcx)
	return(dMS1)
}

# MS1specplot(dbkey,base.dir,finlist,SubDB="FTneg")