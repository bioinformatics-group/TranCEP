#! /usr/bin/Rscript
# hello world~
###################################################
## Name: TranCEP.R
## input: fasta file containing the unknown/testing protein sequences
## output: The predicted class of each protein sequence, and the classes probabilities in csv format
## Author: Munira Alballa
##################################################
library(argparse)

trancepdir <- "."
db<- "./db/"

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

args <- parser$parse_args()

# args <- commandArgs(trailingOnly=TRUE)

# terminate <- FALSE

# out <- "."
# trancepdir <- "."
# db<- "./db/"

# for(i in args){
#    arg = strsplit(i, "=")[[1]];

#    switch(arg[1],
#      "-query"={
#        query <- arg[2]
#      },
#      "-out"={
#        out <- arg[2]
#      },
#      "-trancepdir"={
#        trancepdir <- normalizePath(arg[2])
#      },
#      "-db"={
#          db <- normalizePath(arg[2])
#      },
#      "-help"={
#        cat("TranCEP v1.0 (April 2018)\n")
#        cat("https://doi.org/10.1101/293159\n")
#        cat("\n")
#        cat("Usage: TranCEP -query=<input> [-trancepdir=<trancepdir>] [-out=<outdir>] [-db=<database directory>]\n")
#        cat("\n")
#        cat("\t<input> is your sequence input file in fasta format\n")
#        cat("\t<out> is the output directory where you want the predicted results, formatted as csv\n")
#        cat("\t\t<out> defaults to '",out,"'\n")
#        cat("\t<trancepdir> is the directory where the base TranCEP files are located")
#        cat("\t\t<trancepdir> defaults to '",trancepdir,"'\n")
#        cat("\t\t <db> is the directory where the database is located.")
#        cat("\t\t <db> defaults to ", paste0(trancepdir, "/db/"),"\n")
#        cat("\n")
#        terminate <- TRUE
#        break
#      }
#    )
# }

if(args$query == ""){
  stop("Please input a query file (using -query=path/to/fasta_file).")
}

query <- args$query
out <- args$out
db <- normalizePath(args$db)
trancepdir <- normalizePath(args$trancepdir)

# get directory of the script currently running
script.dir <- dirname(sys.frame(1)$ofile)


test_fasta <- normalizePath(path.expand(query))
resultspath <- paste0(normalizePath(path.expand(out)),"/")

require(seqinr)
library("Biostrings")
library("stringr")
require(protr)
library(ISLR)
library(e1071)
library(caret)
library(R.utils)
wd=normalizePath(path.expand(".")) # change the the tool directory


if (isAbsolutePath(db)){
    dbpath <- db
}else{
    dbpath <- paste0(TooTSCdir, db)
}

compostions=paste0(trancepdir,"/Compositions/")
intermediateFiles=paste0(trancepdir,"/output/")

dir.create(compostions, showWarnings = FALSE, recursive = FALSE, mode = "0777")
testname="test"
#dir.create(paste0(compostions,testname,"/"), showWarnings = TRUE, recursive = FALSE, mode = "0777") # intermediate files go here
substates<- c("amino",    "anion" ,   "cation"  , "electron", "other" ,   "protein" ,"sugar" )


#testing data with unknown substrates

if(file.exists(paste0(trancepdir,"/TCS_MSA_PAAC.R"))){
    source(paste0(trancepdir,"/TCS_MSA_PAAC.R"))
}
else{
    source(paste0(script.dir,"/TCS_MSA_PAAC.R"))
}

load(paste0(trancepdir,"/tranCEP.rda"))

MSA_TCS_PAAC(testname,test_fasta)

testfeatuers = read.csv(paste0(compostions,testname,"_MSA_TCS_PAAC.csv"),sep=",")[,c(-1,-2)]
svm.predtest<-predict(svm.fit,testfeatuers, probability=T)
substateName<- substates[svm.predtest]
seqs<- readFASTA(test_fasta)
probabilities<- attr(svm.predtest,"probabilities")
colnames(probabilities)<- paste(substates,"probability")
names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)) )
print(paste0( "TranCEP output is found at: ", resultspath, "TranCEPout.csv"))
write.csv(cbind(names(seqs),Prediction=substateName, Probabilities=probabilities ),paste0(resultspath,"TranCEPout.csv"))

