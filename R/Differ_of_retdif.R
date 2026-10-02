#' @title Select matches based on retention-time difference consistency
#'
#' @description
#' Helper function that evaluates all possible combinations of candidate
#' FT-QTOF feature matches and selects the combination with the most consistent
#' retention-time differences. Consistency is assessed using the standard
#' deviation of the retention-time differences across the candidate matches.
#'
#' @param Lst A list of candidate FT-QTOF feature matches. Each list element
#'   corresponds to an FT feature and contains one or more candidate QTOF
#'   matches. The sixth column contains the retention-time difference used
#'   to evaluate each combination.
#' @param res.F An object to which the selected FT-QTOF matches are appended.
#'
#' @return The updated `res.F` object containing, for each element of `Lst`,
#'   the candidate belonging to the combination with the lowest standard
#'   deviation of retention-time differences.
#'
#' @noRd
Differ_of_retdif<-function(Lst,res.F){
	lv<-c()
	sd.vector<-c()
	index.matrix<-c()
		#lv contains the number of rows in each of the Lst entries.
	for (i in 1:length(Lst)){
	 a<-length(Lst[[i]][,6])
	 lv<-append(lv,a)														# lv contains the lengths of all Lst components
	}
		#Lst list has 5 entries with each entry consisting of a number of
		#array rows representing isomeric QTOF m/z values. Each of the 5
		#Lst entries represent a FT m/z value and these 5 FT m/z values
		#are neighboring peaks. 
	v<-rep(1,length(Lst))
		#v is a counter, indexing a row of each of the 5 Lst entries.
		#Based on v, 5 retention time differences will be fished from the 
		#6th column of the Lst entries, i.e. between each of the 5 FT m/z
		#values and a particular isomeric QTOF m/z value ['diff.vector'].
		#The standard deviation (sd) will be computed. Each sd should then
		#then be added to the particular combination represented by v.
		#This is done by column binding a finally obtained 'sd.vector'
		#with a finally obtained 'index.matrix' in 'fin.matrix'.
		#lv is the end combination that v should reach. 		
	repeat {
		#vector of retention time differences for a particular v:
	 diff.vector<-c()  
	 for (k in 1:length(Lst)) {
		# v[k] represents index, yields an error if A is absent in syn.o
	  a<-as.numeric(Lst[[k]][v[k],6])											
	  diff.vector<-append(diff.vector,a)
	 }
		#compute sd for the particular v combination:
	 sd.comp<-sd(diff.vector)
	 sd.vector<-append(sd.vector,sd.comp)
		#add particular v combination to 'index.matrix'
	 index.matrix<-rbind(index.matrix,v)
		#next 27 lines of code serve to update the v counter for a next
		#combination for which the sd has to be calculated. Of course,
		#if v reached the lv combination, one should break out of the
		#loop. If none of the v entries is equal to the corresponding lv
		#entry, e.g. if all 5 Lst entries contain multiple rows and you re
		#at the first round in this loop with v=[1 1 1 1 1], then the first
		#element of v will be upgraded with 1 (see the at the end of this
		#code section, at code line 25 from here). More often, the second
		#if construction [any(v==lv)] will be entered.
	 if (all(v==lv)) break	
	 if (any(v==lv)) {
	  tel<-which(v==lv)	# which Lst component reached the end
	  tl<-max(tel)
	  i=which(tel==tl,arr.ind=T)
	  repeat{
		#tel[i]>=2: if the second or any later entry of v is equal to the
		#corresponding entry of the lv vector:
	   if (tel[i]>=2){
		#if the v entry that is located just one before the v entry that
		#reached its maximum value (i.e. the v entry that was equal to the
		#corresponding lv entry), did not yet reach its maximum value,
		#select that one via subtracting 1 from 'i':
	    if (v[tel[i]-1]!=lv[tel[i]-1]){ 
	     if (i==1){
	      v[1]=v[1]+1
	      break
	     }
	     i=i-1
	     next
	    }
	   }
		#if all v entries that are located before the v entry that reached
		#its maximum value (i.e. the v entry that was equal to the corres-
		#ponding lv entry), reached their maximum values, increase the v
		#entry that is located immediately after the v entry that reached
		#its maximum value, by 1. Break out of this v counter upgrade loop
		#and perform again a sd computation. 
	   if (all(v[1:tel[i]]==lv[1:tel[i]])){
	    v[tel[i]+1]=v[tel[i]+1]+1
	    v[1:tel[i]]=1
	    break 
	   }
		#if it is the first entry of v that reached its maximum value,
		#break out of the v counter upgrade loop {i is based on the
		#highest v entry that reached its maximum value
		#"i==which(tel==max(tel))"}:
	   if (i==1) break
		#tel[i] gives v entry that reached its maximum value, but perhaps
		#that entry was >=2, but the v entry just before also reached its
		#maximum value (escaping the v[tel[i]-1]!=lv[tel[i]-1] clausule),
		#but none of the earlier v entries 
		#(escaping the all(v[1:tel[i]]==lv[1:tel[i]]) clausule), then go 
		#to this v entry just before the currently selected v entry tel[i].
		#As this v entry should also be in the tel vector, we have to go to
		#the tel[i-1] entry of the v vector and iterate again through this
		#v counter upgrade loop. This v counter upgrade loop goes towards
		#a first and largest v entry that reaches its maximum value and,
		#then, maximizes the v entries towards the left first and finally
		#those towards the right. 
	   i=i-1
	  }	  
	 }else{
	  v[1]=v[1]+1
	 }
	}
	fin.matrix<-data.frame(cbind(index.matrix,sd.vector))
	o<-order(fin.matrix[[dim(fin.matrix)[2]]])
	fin.matrix<-fin.matrix[o,]
		#fin.matrix is sorted following increasing sd
	for (i in 1:length(Lst)){
		#Lst[[1]] corresponds with column 1 in fin.matrix, 
		#Lst[[2]] with column 2 in fin.matrix and so on.
	 res.int<-Lst[[i]][fin.matrix[1,i],]
	 res.F<-rbind(res.F,res.int)
	}
	return(res.F)
}

# example of a fin.matrix result, the first row is the starting v vector and
# the last row is the lv vector:
#   V1 V2 V3 V4 V5 sd.vector
#1   1  1  1  1  1  3.367040
#2   1  2  1  1  1  3.851974
#3   1  3  1  1  1  3.591889
#4   1  4  1  1  1  3.839794
#5   1  5  1  1  1  3.739402
#6   1  1  1  1  2  3.537899
#7   1  2  1  1  2  3.435571
#8   1  3  1  1  2  3.264228
#9   1  4  1  1  2  3.426769
#10  1  5  1  1  2  3.356714
#11  1  1  1  1  3  3.318842
#12  1  2  1  1  3  3.330045
#13  1  3  1  1  3  3.129515
#14  1  4  1  1  3  3.320026
#15  1  5  1  1  3  3.239434
#16  1  1  1  1  4  3.527234
#17  1  2  1  1  4  3.429442
#18  1  3  1  1  4  3.256821
#19  1  4  1  1  4  3.420586
#20  1  5  1  1  4  3.350065
#21  1  1  1  1  5  3.440613
#22  1  2  1  1  5  3.382778
#23  1  3  1  1  5  3.199309
#24  1  4  1  1  5  3.373465
#25  1  5  1  1  5  3.299004
# first 3 rows of fin.matrix following sorting with increasing sd:
#   V1 V2 V3 V4 V5 sd.vector
#13  1  3  1  1  3  3.129515
#23  1  3  1  1  5  3.199309
#15  1  5  1  1  3  3.239434

#Corresponding Lst example:
#Lst[[1]] corresponds with V1, Lst[[2]] with V2 and so on.
#Lst[[1]]:
#      [,1]    [,2]          [,3]          [,4]          [,5]         
#crd.A "52681" "1.667628333" "162.9397946" "1.608483333" "162.9345297"
#      [,6]       [,7]   
#crd.A "0.059145" "24699"
#Lst[[2]]:
#      [,1]    [,2]       [,3]          [,4]          [,5]         
#crd.A "52564" "1.733975" "422.0253727" "1.19055"     "422.0224726"
#crd.A "52564" "1.733975" "422.0253727" "7.682333333" "422.0227954"
#crd.A "52564" "1.733975" "422.0253727" "6.468383333" "422.0241884"
#crd.A "52564" "1.733975" "422.0253727" "7.631066667" "422.0242155"
#crd.A "52564" "1.733975" "422.0253727" "7.191016667" "422.026423" 
#      [,6]           [,7]   
#crd.A "0.543425"     "24624"
#crd.A "-5.948358333" "24417"
#crd.A "-4.734408333" "24737"
#crd.A "-5.897091667" "24797"
#crd.A "-5.457041667" "24559"
#Lst[[3]]:
#      [,1]    [,2]          [,3]          [,4]          [,5]         
#crd.A "52528" "1.814308333" "162.9397946" "1.608483333" "162.9345297"
#      [,6]       [,7]   
#crd.A "0.205825" "24699"
#Lst[[4]]]:
#      [,1]    [,2]      [,3]          [,4]          [,5]         
#crd.A "52534" "2.34093" "133.0141728" "9.281666667" "133.0088589"
#      [,6]           [,7]   
#crd.A "-6.940736667" "24414"
#Lst[[5]]:
#      [,1]    [,2]          [,3]          [,4]          [,5]         
#crd.C "52676" "2.449023333" "422.0253727" "1.19055"     "422.0224726"
#crd.C "52676" "2.449023333" "422.0253727" "7.682333333" "422.0227954"
#crd.C "52676" "2.449023333" "422.0253727" "6.468383333" "422.0241884"
#crd.C "52676" "2.449023333" "422.0253727" "7.631066667" "422.0242155"
#crd.C "52676" "2.449023333" "422.0253727" "7.191016667" "422.026423" 
#      [,6]           [,7]   
#crd.C "1.258473333"  "24624"
#crd.C "-5.23331"     "24417"
#crd.C "-4.01936"     "24737"
#crd.C "-5.182043334" "24797"
#crd.C "-4.741993334" "24559"
#The retention time difference value of Lst[[4]] is suddenly -6.9 which seems
#outlying and responsible for the selection of the third row of Lst[[2]] with
#retention time difference of -4.7 instead of the first row with retention
#time difference of 0.5 which would be more corresponding to the tR differences
#of Lst[[1]] and Lst[[3]] at 0.06 and 0.20 min. A false positive FTm/z,QTOFm/z
#pair likely exists for Lst[[4]] and seems to be created for Lst[[2]]

#corresponding res.F final list:
#        [,1]    [,2]          [,3]          [,4]          [,5]         
#res.int "52681" "1.667628333" "162.9397946" "1.608483333" "162.9345297"
#res.int "52564" "1.733975"    "422.0253727" "6.468383333" "422.0241884"
#res.int "52528" "1.814308333" "162.9397946" "1.608483333" "162.9345297"
#res.int "52534" "2.34093"     "133.0141728" "9.281666667" "133.0088589"
#res.int "52676" "2.449023333" "422.0253727" "6.468383333" "422.0241884"
#        [,6]           [,7]   
#res.int "0.059145"     "24699"
#res.int "-4.734408333" "24737"
#res.int "0.205825"     "24699"
#res.int "-6.940736667" "24414"
#res.int "-4.01936"     "24737"

