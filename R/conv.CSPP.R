# direc refers to elution order of substrate and product. It is specific for
# a certain conversion: 1, product elute earlier than substrate; 2, product
# elute later than substrate; 3, product might elute earlier or later, but not
# at the same time as the substrate.
# Prod.exp refers to the considered experiment: EXPID
# peakwidth is necessary to create a minimum retention time difference between
# substrate and product.
# min relates to the minimum intensity of the product ions.

conv.CSPP<-function(inp.x,mzdiff,direc,peakwidth,mzerr,Prod.exp,subDB,AnalMS){
  Prod.dat<-inp.x[inp.x[,9]==Prod.exp,]
  Prod.dat<-Prod.dat[order(Prod.dat[,5]),]
  Sub.dat<-Prod.dat
  cspp.df<-data.frame(COMPID.sub=integer(),MZ.sub=numeric(),
                      IONS.sub=integer(),COMPID.prod=integer(),
                      MZ.prod=numeric(),IONS.prod=integer(),
                      COMMON_IONS=integer(),DOT_IONS=numeric(),
                      COMMON_LOSS=integer(),DOT_LOSS=numeric(),
                      FORW_IONS=numeric(),REV_IONS=numeric(),
                      FORW_LOSS=numeric(),REV_LOSS=numeric(),
                      stringsAsFactors=FALSE)
  tR.df<-data.frame(tR.sub=numeric(),tR.prod=numeric(),stringsAsFactors=FALSE)
  i=1
  repeat{
    mz.sub<-Sub.dat[i,5]
    mz.prod<-mz.sub+mzdiff
    prd.low<-mz.prod-mzerr
    prd.high<-mz.prod+mzerr
    j=1
    while(j<=dim(Prod.dat)[1]){
      if(Prod.dat[j,5]<prd.low){
        Prod.dat<-Prod.dat[-j,]
        next
      }
      if((Prod.dat[j,5]>=prd.low)&(Prod.dat[j,5]<=prd.high)){
        out<-targMS2comp(Sub.dat[i,1],Prod.dat[j,1],subDB,AnalMS)
        if((direc==1)&(Prod.dat[j,8]<Sub.dat[i,8]-peakwidth)){
          cspp.df<-rbind(cspp.df,out)
          tR<-c(Sub.dat[i,8],Prod.dat[j,8])
          tR.df<-rbind(tR.df,tR)
        }
        if((direc==2)&(Prod.dat[j,8]>Sub.dat[i,8]+peakwidth)){
          cspp.df<-rbind(cspp.df,out)
          tR<-c(Sub.dat[i,8],Prod.dat[j,8])
          tR.df<-rbind(tR.df,tR)
        }
        if((direc==3)&((Prod.dat[j,8]<Sub.dat[i,8]-peakwidth)|
                       (Prod.dat[j,8]>Sub.dat[i,8]+peakwidth))){
          cspp.df<-rbind(cspp.df,out)
          tR<-c(Sub.dat[i,8],Prod.dat[j,8])
          tR.df<-rbind(tR.df,tR)
        }
        j=j+1
        next
      }
      if(Prod.dat[j,5]>prd.high){
        break
      }
    }
    if(dim(Prod.dat)[1]==0)break
    i=i+1
  }
  return(cspp.df)
}
