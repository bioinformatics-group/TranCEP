###################################################
## Name: TCS_MSA_PAAC.R
## output: TCS_MSA_PAAC feature of given protein sample in  csv format
## requires: TMCOFFEE, BLAST 
## Author: Munira Alballa
##################################################
firstrun=F
blastpSeq<- function(seq, start.pos = 1L, end.pos = nchar(seq), 
                     blastp.path = NULL, makeblastdb.path = NULL, 
                     database.path = NULL, silent = TRUE, 
                     evalue = 10L, output.path=resultspath){
  
  if (Sys.which('makeblastdb') == '' & is.null(makeblastdb.path))
    stop('Please install makeblastdb (included in the NCBI BLAST+) or specify makeblastdb.path.')
  
  if (Sys.which('blastp') == '' & is.null(blastp.path))
    stop('Please install blastp (included in the NCBI BLAST+) or specify blastp.path.')
  
  makeblastdb.path = if (!is.null(makeblastdb.path)) makeblastdb.path else Sys.which('makeblastdb')
  blastp.path = if (!is.null(blastp.path)) blastp.path else Sys.which('blastp')
  
  if (is.null(database.path)) stop('Must specify the database (a FASTA file) path')
  if (is.null(output.path)) stop('Must specify the output path')
  
  N = end.pos - start.pos + 1L
  
  # Prepare data for Blastp
  cmddb = paste0(shQuote(makeblastdb.path), ' -dbtype prot -in ', 
                 shQuote(database.path),' -parse_seqids')        
  if(firstrun)
  {
    print("performing make blastdb")
    # print( cmddb)
    if (silent == TRUE) system(cmddb, ignore.stdout = TRUE) else system(cmddb)
  }
  
  # Basic parameters for Blastp
  tmp = tempfile('Blastp')
  queryFasta = paste0(tmp, '.fasta')
  tabularfile= paste0(tmp, '.txt')
  querySeq = Biostrings::AAStringSet(as.character(seq))
  Biostrings::writeXStringSet(querySeq, queryFasta)          
  # Additional parameters for Blastp
  if (!is.null(evalue)) {
    if (evalue <= 0L) {
      stop('evalue must be > 0')
    }
  }
  # Run Blastp
  
  cmdblastp = paste(
    paste0(shQuote(blastp.path),
           ' -comp_based_stats 1 -db ', shQuote(database.path),
           ' -query ', shQuote(queryFasta),  ' -outfmt 6',' -out ', paste0(shQuote(output.path),"out.txt")))
  
  print("******************************")
  print(cmdblastp)
  if (silent == TRUE) system(cmdblastp, ignore.stdout = F) else system(cmdblastp)      
  #get the hit sequences Id
  if(file.info(paste0(output.path,"out.txt"))$size != 0) # if there some hits are found
  {
    data = read.table(paste0(output.path,"out.txt"))
    HomologousSeqIds= data$V2
    dindex <- which(duplicated(HomologousSeqIds))
    if(length(dindex) !=0 )
    {
      HomologousSeqIds=HomologousSeqIds[-dindex]# remove duplicates if any
    }
    # print(length(HomologousSeqIds))
    
    if(length(HomologousSeqIds)>=120)
      HomologousSeqIds<-HomologousSeqIds[1:120] 
    else
      HomologousSeqIds<- HomologousSeqIds
    
    fileName<-paste0(output.path,"H.txt")
    fileConn<-file(fileName)
    write(as.character(HomologousSeqIds), fileConn)
    close(fileConn)
    
    #get the cossponding Fasta file
    getseqcmd= paste0(shQuote(Sys.which('blastdbcmd')),' -db ',shQuote(database.path), ' -entry_batch ', fileName, ' -out ', paste0(output.path,"/seq.txt"))
    # print(getseqcmd)
    if (silent == TRUE) system(getseqcmd, ignore.stdout = TRUE) else system(getseqcmd)
  }else{
    print (paste0( output.path,"---No hits found"))
  }
}

FilteredMSA= function(path)
{
  #print(path)
  
  setwd(path)
  if(file.exists("seq.txt"))
  {
  tcsScorecmd<-paste0("t_coffee seq.txt -mode psicoffee -blast_server=LOCAL -protein_db ",dbpath,"uniref50-tm.fasta -output tcs_residue_filter3_fasta,clustalw_aln,tcs_column_filter4_fasta,score_html")
  system(tcsScorecmd)
  system('rm *.prf')
  #Removing columns with 80%gaps or more 
  # print("4-Removing columns with 80%gaps or more ")
  removegapscmd= paste0("t_coffee -other_pg seq_reformat -in seq.tcs_column_filter4_fasta", " -output  fasta > filteredSeq.fasta")
  system(removegapscmd)
  }
}

PreparedataforMSAAAC = function(seq, start.pos = 1L, end.pos = nchar(seq), 
                                blastp.path = NULL, makeblastdb.path = NULL, 
                                database.path = NULL, iter = 5, silent = TRUE, 
                                evalue = 10L, output.path=resultspath) {
  #1- run blastp on datafiles
  seqname=names(seq)
  SeqDirectory=paste0(output.path,names(seq),"/")
  if((! dir.exists(SeqDirectory)) ) #first time looking at this sequence
  {
    dir.create(SeqDirectory, showWarnings = FALSE, recursive = FALSE, mode = "0777")   
    blastpSeq(seq, start.pos , end.pos ,  blastp.path, makeblastdb.path , database.path , silent ,evalue, SeqDirectory)
  }
  
  #2- Do Filtered MSA 
  if((! file.exists(paste0(SeqDirectory, "filteredSeq.fasta"))) || (file.info(paste0(SeqDirectory, "filteredSeq.fasta"))$size  ==0))
    {FilteredMSA(SeqDirectory)}
}



MSA_TCS_PAAC<- function(subset,fastafile)
{
  # step#1 perpare data 
  dirName= paste0(intermediateFiles,"/")
  dir.create(dirName, showWarnings = FALSE, recursive = FALSE, mode = "0777") 
  seqs<- readFASTA(fastafile)
  names(seqs)<- sub("\\|.*","",sub(".+?\\|","", names(seqs)))
  for(j in c(1:length(seqs)))
  {
    x<- seqs[j]
    PreparedataforMSAAAC(seq= x,database.path=paste0(dbpath,"/all.fasta"),output.path= dirName)
  }
  
  # at this point, we have N fasta files (where N is the number of blast hits) for each sequence in each set.
  # we have compute the AAC for each of those sequences and then take their mean. 
  # step#2 analyse data ( collect blast hits, do MSA, remove low socore col)
  
  dfMSADC <- data.frame(matrix(ncol = 400+1, nrow =0)) 
   
  setwd(paste0(intermediateFiles))
  AAfiles<-  names(seqs)
  ListofSubsetDetailedMSADCStatistics <- matrix(ncol=400+1,nrow=length(AAfiles))
   rownames( ListofSubsetDetailedMSADCStatistics )<- AAfiles
  print("=====================================================================================")
  print(subset)
  for(i in c(1:length(AAfiles)))
  {
    subdirName<- paste0(getwd(),"/",AAfiles[i],"/")
    print(subdirName)
    
    if((file.exists(paste0(subdirName, "filteredSeq.fasta"))) && (file.info(paste0(subdirName, "filteredSeq.fasta"))$size  >0))
    {
    se<-  readFASTA(paste0(subdirName, "filteredSeq.fasta"))
    seqlist<- unlist(lapply(se,as.character))
    seqlist<-sapply(seqlist,gsub,pattern="[^GPAVLIMCFYWHKRQNEDST]",replacement="") # getting rid of gaps? needs more work
    
    # print("MSADC for one seq") 
    MSADC<- lapply(seqlist, extractDC)
    # colnames( ListofSubsetDetailedMSADCStatistics)<- c(names(MSADC[[1]]))
    #print("MSADC mean") 
    outputDC <- matrix(unlist(MSADC), ncol = 400, byrow = TRUE)
    #colnames(outputDC)<- names(MSADC[[1]])  
    DC<- apply(outputDC, 2, mean)
    ListofSubsetDetailedMSADCStatistics[i,]<-   c("NA",DC)
    }else{ # if the sequence does not have any hits, just compute the PAAC from the single query
      print("filteredSeq.fasta does not exist")
      DC<- extractDC(as.character(seqs[i]))
      ListofSubsetDetailedMSADCStatistics[i,]<-   c("NA",DC)
    }
  }
  dfMSADC<- rbind(dfMSADC,ListofSubsetDetailedMSADCStatistics)
  
  
  write.csv(cbind(UniprotID=AAfiles,dfMSADC), file = paste0(compostions,testname,"_MSA_TCS_PAAC.csv"),row.names = F)
}
