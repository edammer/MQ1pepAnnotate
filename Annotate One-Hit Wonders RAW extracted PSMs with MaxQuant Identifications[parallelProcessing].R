## Single Peptide Identification MS Peptide Spectral Match Annotation Automation Tool
## Eric Dammer - January 18, 2023
#________________________________________________________________________

#Requires:
# R v4.2.0 or newer
# OS with mono-devel (RHEL) installed ( https://www.mono-project.com/ )
# Successful installation of "rawrr" package ( https://www.biorxiv.org/content/10.1101/2020.10.30.362533v1.full OR https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00866 )
# Installation of "protViz", "doparallel", and "stringr" packages

#Inputs:
# Path to Thermo RAW files searched by MaxQuant to obtain the search result text files.
# The above path should have a subfolder named 'txt', containing MaxQuant outputs proteinGroups.txt, peptides.txt, msmsScans.txt, and evidence.txt

#Outputs:
# Table of "one-hit wonder" single peptide (possibly multi-PSM) identifications of protein groups (readable by this script as well)
# PDF of annotated extracted MS/MS PSMs (the top, best scoring, for each identification in the above output table)

#________________________________________________________________________

## One-time setup of Thermo RAW file reading capability in native OS
BiocManager::install("rawrr")
rawrr::installRawFileReaderDLLs()  #Accept license at prompt.
rawrr::installRawrrExe()

#BiocManager::install("ExperimentHub")
BiocManager::install("protViz")

## Configuration section:

rootdir="~/SydneyRAW/IP/"    #Contains RAW files and MaxQuant ./txt/ subfolder in it. The files proteinGroups.txt , msms.txt , and evidence.txt must be present

#####################

## Step 1) Compile "One-hit Wonder" Table for export and also for selection of needed MS/MS spectra to be annotated in step 2

setwd(rootdir)

# Read proteinGroups.txt
proteinGroups<-read.delim(file=paste0(rootdir,"/txt/proteinGroups.txt"),sep="\t",header=TRUE,check.names=FALSE) #,row.names=1 #decoys can have same ID
dim(proteinGroups)
decoyCount=length(which(proteinGroups$Reverse=="+"))
pg.originalRowCount=nrow(proteinGroups)
decoyIndices=which(proteinGroups$Reverse=="+")
FDR=paste0(round(length(decoyIndices)/nrow(proteinGroups)*100,2),"% FDR")
#LFQindices=which(grepl("Intensity.", colnames(proteinGroups)))
cat(paste0("Imported MaxQuant proteinGroups.txt has ",pg.originalRowCount,"x",ncol(proteinGroups)," rows x columns; ",decoyCount," reverse hits, for a net ",FDR,".\n","Summed intensity is available for ",length(LFQindices)," experiment samples.\n"))
proteinGroups<-proteinGroups[which(!proteinGroups$Reverse=="+"),]  #remove decoys
colnames(proteinGroups)
keep.proteinGroupCols.idx<-which(colnames(proteinGroups) %in% c("Majority protein IDs", "Protein names", "Gene names", "Fasta headers", "Peptide IDs", "Best MS/MS"))
masterOHWtable.1<-proteinGroups[which(proteinGroups$Peptides==1), keep.proteinGroupCols.idx]
nrow(masterOHWtable.1)  #count of one-hit wonder peptides in the file

# Read peptides.txt
peptides <-read.delim(file=paste0(rootdir,"/txt/peptides.txt"),sep="\t",header=TRUE, check.names=FALSE)
dim(peptides)
peptides <- peptides[which(peptides$id %in% masterOHWtable.1$'Peptide IDs'),]  #reduce rows to those relevant to one-hit wonders
masterOHWtable.1$Sequence<-peptides$Sequence[match(masterOHWtable.1$'Peptide IDs',peptides$id)]

# Read evidence.txt; fallback to msmsScans.txt if needed
evidence <- read.delim(file=paste0(rootdir,"/txt/evidence.txt"),sep="\t",header=TRUE,check.names=FALSE)
evidence <- evidence[which(evidence$'MS/MS IDs' %in% masterOHWtable.1$'Best MS/MS'), which(colnames(evidence) %in% c("Modified sequence", "Raw file", "Charge", "MS/MS scan number", "m/z", "Mass", "Mass error [ppm]", "Mass error [Da]", "Retention time", "Retention length", "Score", "Delta score", "MS/MS IDs")),]  #reduce rows and columns to those relevant
missingEvidence.idx=which(!masterOHWtable.1$'Best MS/MS' %in% evidence$'MS/MS IDs')
if (length(missingEvidence.idx)>0) {
  msmsScans <- read.delim(file=paste0(rootdir,"/txt/msmsScans.txt"),sep="\t",header=TRUE,check.names=FALSE)
  msmsScans2 <- msmsScans[which(msmsScans$'MS/MS IDs' %in% masterOHWtable.1$'Best MS/MS'[missingEvidence.idx]), c("Modified sequence", "Raw file", "Charge", "Scan number", "m/z", "Mass", "Retention time", "Score", "MS/MS IDs")]  # Read available data from msmsScans.txt missing in evidence.txt
  colnames(msmsScans2) <- c("Modified sequence", "Raw file", "Charge", "MS/MS scan number", "m/z", "Mass", "Retention time", "Score", "MS/MS IDs")
  #msmsScans3 <- msmsScans2[, match(colnames(evidence),colnames(msmsScans2))]
  blankMat<-matrix(NA,nrow=nrow(msmsScans2),ncol=ncol(evidence))
  colnames(blankMat)<-colnames(evidence)
  evidence<-as.data.frame(rbind(evidence, blankMat))
  for (columnName in colnames(msmsScans2)) evidence[c((nrow(evidence)-(nrow(msmsScans2)-1)):nrow(evidence)),columnName] <- msmsScans2[,columnName]
}
masterOHWtable.2 <-  evidence[match(masterOHWtable.1$'Best MS/MS',evidence$'MS/MS IDs'),]

# Finalize OHW Table
masterOHWtable<-cbind(masterOHWtable.1, masterOHWtable.2)
rownames(masterOHWtable)<-NULL
write.table(masterOHWtable,file="Master One-peptide Protein Identification MS2 Spectra Annotations (via protViz)-TOTAL_31_RAW_data.txt",sep="\t",row.names=FALSE)


## Step 2) Build PDFs of PSM annotations for each row in table compiled in previous step. [ PARALLEL, IN CHUNKS ]

library(doParallel)
clusterLocal <- makeCluster(c(rep("localhost",32)),type="SOCK")
registerDoParallel(clusterLocal)
parallel::clusterExport(cl=clusterLocal, varlist=c("rootdir"), env=globalenv())

start_time <- Sys.time()

writePDFforChunk.df <- function(chunk.df) {
	require("protViz")
	require("rawrr")
	pdf(file=paste0("Annotated MS2 spectra of one-peptide protein identifications (Input Sample 32 RAW data fileSearch)[",chunk.df$'Raw file'[1],"].pdf"),width=11,height=8.5)
	par(oma=c(0,1.5,3,2))
	for (tableEntry in 1:nrow(chunk.df)) {
		filenames.wd<-list.files(rootdir)
		rawfile=filenames.wd[which(grepl(paste0(chunk.df[tableEntry,"Raw file"], "\\.[Rr][Aa][Ww]"),filenames.wd))]
		thisScan=chunk.df[tableEntry,"MS/MS scan number"]
		thisPepSeq=gsub("_", "", chunk.df[tableEntry,"Modified sequence"])
		pepSeq.forPlotTitle=gsub("\\(Oxidation \\(M\\)\\)","*",
			gsub("\\(Deamidation \\(NQ\\)\\)","^",
			gsub("\\(Acetyl \\(Protein N-term\\)\\)","(Ac)",thisPepSeq) ))
		primaryPepSeq=chunk.df[tableEntry,"Sequence"]
		varMassShift.vec<-rep(0.0, nchar(primaryPepSeq))
		
		if (unlist(regexpr('\\(Acetyl \\(Protein N-term\\)\\)', thisPepSeq))[1] == 1) {
		  varMassShift.vec[1]=varMassShift.vec[1]+42.010565				#6-digit monoisotopic precision from http://www.unimod.org/
		  thisPepSeq=gsub('\\(Acetyl \\(Protein N-term\\)\\)', '', thisPepSeq)
		}
		
		
		#Handle M(Oxidation) residue mass shifts.
		if (unlist(regexpr('\\(Oxidation \\(M\\)\\)', thisPepSeq))[1] > 0) {
		  oxPTMcount=nrow(stringr::str_locate_all(thisPepSeq, '\\(Oxidation \\(M\\)\\)')[[1]])
		  for (i in 1:oxPTMcount) {
		    thisPosition=unlist(regexpr('\\(Oxidation \\(M\\)\\)', 
			gsub("\\(Deamidation \\(NQ\\)\\)","",thisPepSeq)))[1] -1  #remove other mod annotations which will violate single letter residue notation assumption while counting character positions.
		    varMassShift.vec[thisPosition] <- varMassShift.vec[thisPosition] +15.994915
		    thisPepSeq=sub('\\(Oxidation \\(M\\)\\)', '', thisPepSeq)     #remove first instance match only, so that for loop can continue to find positions of later matches within single letter residue notation
		  }
		}
		
		#Handle NQ(Deamidation residue mass shifts.
		if (unlist(regexpr('\\(Deamidation \\(NQ\\)\\)', thisPepSeq))[1] > 0) {
		  deamPTMcount=nrow(stringr::str_locate_all(thisPepSeq, '\\(Deamidation \\(NQ\\)\\)')[[1]])
		  for (i in 1:deamPTMcount) {
		    thisPosition=unlist(regexpr('\\(Deamidation \\(NQ\\)\\)', 
		        gsub("\\(Oxidation \\(M\\)\\)","",thisPepSeq)))[1] -1     #remove other mod annotations which will violate single letter residue notation assumption while counting character positions.
		    varMassShift.vec[thisPosition] <- varMassShift.vec[thisPosition] +0.984016
		    thisPepSeq=sub('\\(Deamidation \\(NQ\\)\\)', '', thisPepSeq)  #remove first instance match only, so that for loop can continue to find positions of later matches within single letter residue notation
		  }
		}
		
		#Sanity Check (have all PTMs been handled for this peptide?
		if(!thisPepSeq==primaryPepSeq) stop(paste0("Peptide entry ", tableEntry, " stripped of PTM notations is not equivalent to primary sequence for the same entry.\nWere all modifications addressed in the code to get mass shifts by residue position?\n\n"))
		
		# determine the fragment ions with mass shifts
		thisPep.AAweights <-aa2mass(primaryPepSeq)[[1]] + varMassShift.vec
		thisIonSeries = fragmentIon(thisPep.AAweights)[[1]] #[,c("b","y")]  #function would also work on primary pep sequence, but does not directly handle mass shifts(?)
		
	#	PSM[[tableEntry]] <-
		rawrr::readSpectrum(rawfile=rawfile, scan=thisScan) |>
		    lapply(function(x) protViz::peakplot(peptideSequence=thisPepSeq, fi=thisIonSeries, pattern.abc="[ab].*", pattern.xyz="[y].*", x))
	
	#	thisMSMS[[tableEntry]] <- rawrr::readSpectrum(rawfile=rawfile, scan=scan)
	#	protViz::peakplot(peptideSequence = thisPepSeq, fi=thisIonSeries, thisMSMS[[tableEntry]])
	
		mtext(pepSeq.forPlotTitle, side=3,outer=TRUE,adj=0.5,cex=1.75)
	}  # ends for (tableEntry ...
	dev.off()
}  #close function before calling it after %dopar% below.

foreach::foreach(chunk.df=split(masterOHWtable, masterOHWtable$'Raw file')) %dopar% writePDFforChunk.df(chunk.df)

end_time <- Sys.time()
end_time - start_time


#monoisotopic mass definitions
AA
#   letter1 letter3   formula Monoisotopic  Average
#1        A     Ala    C3H5ON     71.03711  71.0788
#2        R     Arg  C6H12ON4    156.10111 156.1875
#3        N     Asn  C4H6O2N2    114.04293 114.1038
#4        D     Asp   C4H5O3N    115.02694 115.0886
#5        C     Cys   C3H5ONS    103.00919 103.1388  #+57 static
#6        E     Glu   C5H7O3N    129.04259 129.1155
#7        Q     Gln  C5H8O2N2    128.05858 128.1307
#8        G     Gly    C2H3ON     57.02146  57.0519
#9        H     His   C6H7ON3    137.05891 137.1411
#10       I     Ile   C6H11ON    113.08406 113.1594
#11       L     Leu   C6H11ON    113.08406 113.1594
#12       K     Lys  C6H12ON2    128.09496 128.1741
#13       M     Met   C5H9ONS    131.04049 131.1926
#14       F     Phe    C9H9ON    147.06841 147.1766
#15       P     Pro    C5H7ON     97.05276  97.1167
#16       S     Ser   C3H5O2N     87.03203  87.0782
#17       T     Thr   C4H7O2N    101.04768 101.1051
#18       W     Trp C11H10ON2    186.07931 186.2132
#19       Y     Tyr   C9H9O2N    163.06333 163.1760
#20       V     Val    C5H9ON     99.06841  99.1326

# aa2mass() and parentIonMass() are useful calculaton functions in the protViz package


