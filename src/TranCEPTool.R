#! /usr/bin/Rscript
# hello world~
###################################################
## Name: TranCEP.R
## input: fasta file containing the unknown/testing protein sequences
## output: The predicted class of each protein sequence, and the classes probabilities in csv format
## Author: Munira Alballa
##################################################
library(argparse)
require(seqinr)
library("Biostrings")
library("stringr")
require(protr)
library(ISLR)
library(e1071)
library(caret)
library(R.utils)

parser <- ArgumentParser()

# positional/ mandatory argument
parser$add_argument("-query",default="",
                    help="The sequence input file in fasta format.",
                    metavar="Query file (fasta file path)")
parser$add_argument("-out", "-o", default='.', 
                    help="Output directory where you want the predicted results, formatted as csv. Default [%(default)].",
                    metavar="DIRECTORY")
parser$add_argument("-trancepdir", default='.', 
                    help="The directory where the base TranCEP files are located. Default [%(default)].",
                    metavar="DIRECTORY")
parser$add_argument("-db", default="./db/",
                    help="The directory where the database is located. Default [%(default)].", metavar="DIRECTORY")
args <- parser$parse_args() #start the parser
if(args$query == ""){
  stop("Please input a query file (using -query=path/to/fasta_file).")
}

query <- args$query
out <- args$out
db <- normalizePath(args$db)
trancepdir <- normalizePath(args$trancepdir)


if (isAbsolutePath(db)){
    dbpath <- db
}else{
    dbpath <- paste0(trancepdir, db)
}

test_fasta <- normalizePath(path.expand(query))
resultspath <- paste0(normalizePath(path.expand(out)),"/")


wd=normalizePath(path.expand(".")) # working directory

# create necessary directory
compostions=paste0(trancepdir,"/Compositions/")
intermediateFiles=paste0(trancepdir,"/output/")
dir.create(compostions, showWarnings = FALSE, recursive = FALSE, mode = "0777")
testname="test"

#dir.create(paste0(compostions,testname,"/"), showWarnings = TRUE, recursive = FALSE, mode = "0777") # intermediate files go here
substates<- c("amino",    "anion" ,   "cation"  , "electron", "other" ,   "protein" ,"sugar" )


#testing data with unknown substrates

# try in trancepdir/src/ if the script is not found in trancepdir 
if(file.exists(paste0(trancepdir,"/TCS_MSA_PAAC.R"))){
    source(paste0(trancepdir,"/TCS_MSA_PAAC.R"))
}else{
    source(paste0(trancepdir,"/src/TCS_MSA_PAAC.R"))
}
# load date
load(paste0(trancepdir,"/tranCEP.rda"))

MSA_TCS_PAAC(testname,test_fasta)

testfeatuers = read.csv(paste0(compostions,testname,"_MSA_TCS_PAAC.csv"),sep=",")[,c(-1,-2)]
svm.predtest <- predict(svm.fit,testfeatuers, probability=T)
substateName <- substates[svm.predtest]
seqs <- readFASTA(test_fasta)
probabilities <- attr(svm.predtest,"probabilities")
colnames(probabilities) <- paste(substates,"probability")
names(seqs) <- sub("\\|.*","",sub(".+?\\|","", names(seqs)) )
print(paste0( "TranCEP output is found at: ", resultspath, "TranCEPout.csv"))

write.csv(cbind(names(seqs),Prediction=substateName, Probabilities=probabilities ),paste0(resultspath,"TranCEPout.csv"))

