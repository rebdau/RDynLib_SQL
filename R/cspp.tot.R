# No rows representing the different COMPIDs of the experiment have to be
# added to compound_add.txt, this is done while running the function. In case
# no conversions are found for a particular COMPID, a row of NAs is added
# preceeded by the name of the COMPID, the function jumps then to the next
# COMPID. However, at the end of the COMPID list belonging to the experiment,
# some final expected rows might be missing whenever no conversions were
# found for the few last COMPIDs. This is not a problem as these last
# COMPIDs followed by NAs will be added upon running the next experiment.
# Note: Conversions for a particular COMPID are added to the row with row
# number equal to the COMPID.

cspp.tot<-function(base.dir,finlist,SubDB,Prod.exp,inp.x=NULL,peakwidth=NULL,
                   mzerr=NULL,comp_add=NULL){
  if (is.null(peakwidth)) peakwidth<-0.2
  TMP<-SelectSubDB(base.dir,finlist,SubDB)
  subDB<-TMP[[1]]
  AnalMS<-TMP[[2]]
  if (is.null(mzerr)) {
    if (AnalMS=="FT") {
      mzerr<-0.01
    } else if (AnalMS=="QTOF"){	
      mzerr<-0.015
    }
  }
  if (is.null(inp.x)) {
    setwd(TMP[[3]])
    inp.x<-read.table("compound.csv",header=T,sep="\t",stringsAsFactor=F)
    # has problems when primes are used, replace them by pr
  }
  if (is.null(comp_add)) {
    setwd(TMP[[4]])
    comp_add<-read.table("compound_add.txt",header=T,sep="\t",
                         stringsAsFactors=FALSE)
    # generating a data frame with conversion type, mass difference
    # (mzdiff), elution order (direc) and column to write the results
    # to (conv.col).
  }
  conv.table<-read.table("cspp.txt",header=T,sep="\t",stringsAsFactors=FALSE)
  conv.type<-as.character(conv.table[,1])
  conv.mz<-as.numeric(conv.table[,3])
  conv.direc<-as.integer(conv.table[,5])
  conv.colmn<-as.integer(conv.table[,6])
  conver<-data.frame(conv.type,conv.mz,conv.direc,conv.colmn)
  k=1
  repeat{
    mzdiff=conver[k,2]
    direc=conver[k,3]
    cspp.df<-conv.CSPP(inp.x,mzdiff,direc,peakwidth,mzerr,Prod.exp,
                       subDB,AnalMS)
    conv.col=conver[k,4]
    comp_add<-rank.cspp(cspp.df,conv.col,comp_add)
    if(k==dim(conver)[1])break
    k=k+1
  }
  comp_add[,1]<-inp.x[1:(dim(comp_add)[1]),1]
  return(comp_add)
}

# cspp_add<-cspp.tot(base.dir,finlist,SubDB="FTneg",Prod.exp=2)
# write.table(cspp_add,"compound_add.txt",sep="\t",row.names=F)