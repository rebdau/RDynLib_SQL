NeutLossMatch<-function(neutprod.num,neutloss,err){
	i=1
	repeat{
	 lowmass<-neutloss[i,1]-err
	 highmass<-neutloss[i,1]+err
	 j=1
	 repeat{
	  if((neutprod.num[j]>=lowmass)&(neutprod.num[j]<=highmass)){
	   infostr<-paste(neutprod.num[j],"Da --",neutloss[i,2],
				"--",neutloss[i,3],sep=" ")
	   print(infostr)
	  }
	  if(j==length(neutprod.num))break
	  j=j+1
	 }
	 if(i==dim(neutloss)[1])break
	 i=i+1
	}
	writeLines("\n")
}

