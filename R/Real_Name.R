Real_Name<-function(HasName,sel,Assoc_line,Names_lst,SBDB){
	Nr.Name<-length(which(HasName%in%sel))
	Mult.names<-as.character(Assoc_line[1])
	if (Nr.Name==1) {
	 sel.f<-which(HasName%in%sel)	# Assoc column for which the Assoc_line entry has also a name
	 sel.r<-sel[-sel.f]		# Assoc column(s) for which the aligned Assoc_line entries have NULL or ! as name
	 Real.Compid<-Names_lst[[HasName[sel.f]]][Assoc_line[HasName[sel.f]],1] 
	 	# Yields in the Assoc_line the COMPID of the sub-database for which a name exists
		# [1] 11026
	 Real.subDB<-HasName[sel.f]
	 	# sel.f is index of HasName, but HasName[sel.f] is necessary as index 
	 	# for subdb selection.
		# [1] 1
	 Real.hist<-paste(as.character(SBDB[Real.subDB]),as.character(Real.Compid),sep=" ")
		# [1] "FTneg 11026" 
		# this is later entered in the SUBSID columns of the other sub-databases to which the name was transferred.
	 Real.Name<-Names_lst[[HasName[sel.f]]][Assoc_line[HasName[sel.f]],2]
		# [1] "(epi)catechin-4_8pr-(epi)catechin-4pr_8prpr-(epi)catechin 3 AN KM"
	 Real.smiles<-Names_lst[[HasName[sel.f]]][Assoc_line[HasName[sel.f]],4]
		# [1] "C12=CC(=CC(=C1C(C(C(O2)C=3C=C(C(=CC3)O)O)O)C=4C5=C(C(=CC4O)O)C(C(C(O5)C=6C=C(C(=CC6)O)O)O)C=7C8=C(C(=CC7O)O)CC(C(O8)C=9C=CC(=C(C9)O)O)O)O)O"
	 for (k in sel.r) {
		# All data of the COMPID 11026 in subDB 1, i.e. FTneg, is entered into the other sub-databases for which a corresponding COMPID could be aligned.
	  Names_lst[[k]][Assoc_line[k],2]<-Real.Name
	  Names_lst[[k]][Assoc_line[k],3]<-Real.hist
	  Names_lst[[k]][Assoc_line[k],4]<-Real.smiles
	 }
	}else{
	 for (k in 1:3) {
	  Mult.names<-append(Mult.names,Names_lst[[k]][Assoc_line[k],2]) # names of the aligned COMPIDS in each of the sub-databases
		# Mult.names
	  	# [1] "11026"                                                            
		# [2] "(epi)catechin-4_8pr-(epi)catechin-4pr_8prpr-(epi)catechin 3 AN KM"
		# [3] NA                                                                 
		# [4] "(epi)catechin-4_8pr-(epi)catechin-4pr_8prpr-(epi)catechin 3 AN KM"
	 }
	}
#	print(Mult.names)
	return(list(Names_lst,Mult.names))
}
