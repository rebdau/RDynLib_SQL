# change the name of someone's compound.csv file to compound2.csv in the
# windows explorer. The same for the compound_name.txt file.
# The same rows in the compound.csv, compound2.csv, compound_name.txt and
# compound_name2.txt files refer to the same COMPID; check this in advance

Egal<-function(base.dir,finlist,SubDB,file2=NULL,file4=NULL){
	if (is.null(file2)) file2="compound2.csv"
	if (is.null(file4)) file4="compound_name2.txt"
	options(stringsAsFactors=F)
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	setwd(TMP[[4]])
	file5<-read.table("compound_name.txt",header=T,sep="\t",stringsAsFactors=F)
	file5<-file5[order(file5[,1]),]
	file6<-read.table(file4,header=T,sep="\t",stringsAsFactors=F)
	file6<-file6[order(file6[,1]),]
	setwd(TMP[[3]])
	file1<-read.table("compound.csv",header=T,sep="\t",stringsAsFactors=F)
	file1<-file1[order(file1[,1]),]
	file3<-read.table(file2,header=T,sep="\t",stringsAsFactors=F)
	file3<-file3[order(file3[,1]),]
	sel1<-which(file1[,3]!=file3[,3])		# COMPNAME
	sel2<-which(file1[,4]!=file3[,4])		# FORMULA
	sel3<-which(file1[,7]!=file3[,7])		# PPM_DEVIATION
	sel4<-which(file1[,13]!=file3[,13])		# SMILES
	sel5<-which(file5[,2]!=file6[,2])		# CONSERVATIVE_NAME

	sl<-data.frame(file1[sel1,1],file1[sel1,5],file1[sel1,8],file1[sel1,3],file3[sel1,3])
	write.table(sl,"egal_COMPNAME.csv",sep="\t",row.names=F)
	print(sl)
	print("data frame has been stored as egal_COMPNAME.csv in the CSV directory.")
	print("use this .csv file to check all and afterwards remove it from the CSV directory.")
	print("Give entry numbers that should be removed.")
	print("Separate them by comma, e.g.: 1,3,7")
	print("In case no entries have to be removed: enter 0.")
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Entries:"))
	if(rm_sel!="0"){
	 rm_sel<-as.integer(unlist(strsplit(rm_sel,split=",")))
	 sel1<-sel1[-rm_sel]
	}
	sl<-data.frame(file1[sel1,1],file1[sel1,5],file1[sel1,8],file1[sel1,3],file3[sel1,3])
	print(sl)
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Satisfied? (y/n)"))
	if (rm_sel=="n") return()

	sl<-data.frame(file1[sel2,1],file1[sel2,5],file1[sel2,8],file1[sel2,3:4],file3[sel2,3:4])
	write.table(sl,"egal_FORMULA.csv",sep="\t",row.names=F)
	print(sl)
	print("data frame has been stored as egal_FORMULA.csv in the CSV directory.")
	print("use this .csv file to check all and afterwards remove it from the CSV directory.")
	print("Give entry numbers that should be removed.")
	print("Separate them by comma, e.g.: 1,3,7")
	print("In case no entries have to be removed: enter 0.")
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Entries:"))
	if(rm_sel!="0"){
	 rm_sel<-as.integer(unlist(strsplit(rm_sel,split=",")))
	 sel2<-sel2[-rm_sel]
	}
	sl<-data.frame(file1[sel2,1],file1[sel2,5],file1[sel2,8],file1[sel2,3:4],file3[sel2,3:4])
	print(sl)
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Satisfied? (y/n)"))
	if (rm_sel=="n") return()

	sl<-data.frame(file1[sel3,1],file1[sel3,5],file1[sel3,8],file1[sel3,7],file3[sel3,7])
	write.table(sl,"egal_PPM.csv",sep="\t",row.names=F)
	print(sl)
	print("data frame has been stored as egal_PPM.csv in the CSV directory.")
	print("use this .csv file to check all and afterwards remove it from the CSV directory.")
	print("Give entry numbers that should be removed.")
	print("Separate them by comma, e.g.: 1,3,7")
	print("In case no entries have to be removed: enter 0.")
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Entries:"))
	if(rm_sel!="0"){
	 rm_sel<-as.integer(unlist(strsplit(rm_sel,split=",")))
	 sel3<-sel3[-rm_sel]
	}
	sl<-data.frame(file1[sel3,1],file1[sel3,5],file1[sel3,8],file1[sel3,7],file3[sel3,7])
	print(sl)
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Satisfied? (y/n)"))
	if (rm_sel=="n") return()

	sl<-data.frame(file1[sel4,1],file1[sel4,5],file1[sel4,8],file1[sel4,13],file3[sel4,13])
	write.table(sl,"egal_SMILES.csv",sep="\t",row.names=F)
	print(sl)
	print("data frame has been stored as egal_SMILES.csv in the CSV directory.")
	print("use this .csv file to check all and afterwards remove it from the CSV directory.")
	print("Give entry numbers that should be removed.")
	print("Separate them by comma, e.g.: 1,3,7")
	print("In case no entries have to be removed: enter 0.")
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Entries:"))
	if(rm_sel!="0"){
	 rm_sel<-as.integer(unlist(strsplit(rm_sel,split=",")))
	 sel4<-sel4[-rm_sel]
	}
	sl<-data.frame(file1[sel4,1],file1[sel4,5],file1[sel4,8],file1[sel4,13],file3[sel4,13])
	print(sl)
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Satisfied? (y/n)"))
	if (rm_sel=="n") return()

	sl<-data.frame(file1[sel5,1],file1[sel5,5],file1[sel5,8],file5[sel5,1],file5[sel5,2],file6[sel5,2])
	write.table(sl,"egal_CONSERVAT.csv",sep="\t",row.names=F)
	print(sl)
	print("data frame has been stored as egal_CONSERVAT.csv in the CSV directory.")
	print("use this .csv file to check all and afterwards remove it from the CSV directory.")
	print("Give entry numbers that should be removed.")
	print("Separate them by comma, e.g.: 1,3,7")
	print("In case no entries have to be removed: enter 0.")
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Entries:"))
	if(rm_sel!="0"){
	 rm_sel<-as.integer(unlist(strsplit(rm_sel,split=",")))
	 sel5<-sel5[-rm_sel]
	}
	sl<-data.frame(file1[sel5,1],file1[sel5,5],file1[sel5,8],file1[sel4,13],file3[sel4,13])
	print(sl)
	rm_sel=""
	while (rm_sel=="") rm_sel<-c(readline("Satisfied? (y/n)"))
	if (rm_sel=="n") return()

	if(length(sel1)!=0){
	 for (i in 1:length(sel1)) {
	  file1[sel1[i],3]<-file3[sel1[i],3]
	 }
	}
	if(length(sel2)!=0){
	 for (i in 1:length(sel2)) {
	  file1[sel2[i],4]<-file3[sel2[i],4]
	 }
	}
	if(length(sel3)!=0){
	 for (i in 1:length(sel3)) {
	  file1[sel3[i],7]<-file3[sel3[i],7]
	 }
	}
	if(length(sel4)!=0){
	 for (i in 1:length(sel4)) {
	  file1[sel4[i],13]<-file3[sel4[i],13]
	 }
	}
	if(length(sel5)!=0){
	 for (i in 1:length(sel5)) {
	  file5[sel5[i],2]<-file6[sel5[i],2]
	 }
	}
	write.table(file1,"compound.csv",sep="\t",row.names=F)
	setwd(TMP[[4]])
	write.table(file5,"compound_name.txt",sep="\t",row.names=F)	
	return(list(file1,file5))
}

# file<-Egal(base.dir,finlist,SubDB="QTOFneg")


