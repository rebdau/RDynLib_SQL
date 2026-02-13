folders <- "D:/RDynLib/test"
#folders <- c("D:/These/Ahlam/XCMS/flaxseed/S1", "D:/These/Ahlam/XCMS/flaxseed/S9")
fls <- list.files(path =  folders, pattern = "\\.mzXML$", full.names = TRUE, recursive = TRUE)

library(MSnbase)
library(CompoundDb)
library(RSQLite)

## Create an empty CompDb and store the database in a temporary file
cdb <- emptyCompDb(tempfile())

#Create object of class MSnExp
msdata <- readMSData(files = fls, mode = "onDisk")
#follow workflow to create feature matrix with chromatographic peaks

#For each Feature in the Feature matrix, individually: extract optimal associated Spectra objects and add to database
#Use Filter functions to subset MSnExp object before extracting Spectra objects
##select MS1 spectrum with highest total ion current
##select MS2 spectrum with highest precursor intensity
#Filter functions to be applied on MSnExp object: filterRT(), filterMsLevel(), filterPrecursorMz(), filterFile, ...
##...see for more filters later
msdata.filtered <- 
sps <- msdata %>% filterFile(file=1) %>% filterRt(rt=c(180,181)

#Create dataframe with spectra from filtered MSnExp:
#create list of Spectrum objects from filtered (reduced) MSnExp object 
spc <- spectra(msdata.filtered)
#and set backend to dataframe
setBackend(,MsBackendDataFrame())
#OR extract spectra data as dataframe with spectraData()
spectraData(object, columns = spectraVariables(object))

spc <- spectra(sps)


#Create a compound_id (DYNLIB0000000001) for each feature for which spectra were extracted
#Create compound tibble with compound_id and columns of DynLib compound.csv
#from feature matrix: NODENAME: M...T...;
#for now, set everything except compound_id and nodename to NA

cmp <- data.frame(
compound_id = ,
nodename = ,
name = NA,
formula = NA,
exactmass = NA,
...)

#? can different metadata be added to each experiment in a CompounDb? Maybe not
#Easier will be to just add the metadata corresponding to previous "experiment.csv" file to the cmp dataframe
#columns are EXPID	DATE	USER	MACHINE	COLUMN	MSTYPE	BUFFERA	BUFFERB	GRADIENT_TIME	SOURCE	MODE	SPECIES	TISSUE	CE	META
#maybe CE (= collision energy) can be ommitted as it should be in the spectra data

## add these compounds to the existing (empty) database with insertCompound
cdb <- insertCompound(cdb, compounds = cmp, )#can also add spectra here?

#add compound_id to dataframe with spectra
#add spectra to cdb



