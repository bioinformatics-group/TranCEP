# TranCEP

This tool predicts the substrate class of a given transporter. The class can belong to one of the following seven categories (ordered alphabetically)
```
[1] "amino"    [2] "anion"    [3] "cation"   [4] "electron" [5] "other"    [6] "protein" 
[7] "sugar"
```
 
Input: transporter proteins sequences in Fasta format
Output: the substrate class with highest probability, and the probabilities of the other classes.

## FOLDERS
There are a number of folders that support the running of TranCEP and its outputs.

### Dataset
Contains both the training and the independent testing dataset from the benchmark dataset from TrSSP server available at http://bioinfo.noble.org/TrSSP

### Compositions
Contains the extracted `MSA_TCS_PAAC` features of each protein in the training and independent datasets. The features are saved in csv format, the first column indicates the unitportID, the second column contains the substrate class (encoded as numbers from 1:7 according to their alphabetical order) and the next 400 columns contains the features.

### output
Contains the homology details needed to extract the features. Details of the `MSA`, `TOPCONS scores`, `Blast hits` for each sequence is found here.

This constitutes large amounts of generated data and thus is not included in this repository. You may download the 170M tar.gz file from [here](http://tootsuite.encs.concordia.ca/TranCEP/output.tar.gz) and run `tar -xzf output.tar.gz` from within the output directory to extract these files.

### db
Contains the database to be used.

This constitutes large amounts of static data and thus is not included in this repository. You may download the 290M tar.gz file from [here](http://tootsuite.encs.concordia.ca/TranCEP/db.tar.gz) and run `tar -xzf db.tar.gz` from within the db directory to extract these files.

### src
The scripts needed to use the tool.

## HOW TO USE
 - This tool requires that `TM-COFFEE` and `BLAST` be pre-installed
 - Usage: `TranCEP -query=<input> [-trancepdir=<trancepdir>] [-out=<outdir>]`
  - `<input>` is your sequence input file in fasta format
  - `<out>` is the output directory where you want the predicted 	results, formatted as csv
  - `<trancepdir>` is the directory where the base TranCEP files 	are located
 - `MSA_TCS_PAAC` features of each sequence in the test set is  found under [Compositions/MSA_TCS_PAAC.csv](Compositions/MSA_TCS_PAAC.csv)

## A SIMPLIFIED SETUP
To quickly run this on any environment that supports singularity (and wget):
```
git clone https://github.com/bioinformatics-group/TranCEP.git
cd TranCEP/db
wget http://tootsuite.encs.concordia.ca/TranCEP/db.tar.gz
tar -xzf db.tar.gz
rm -f db.tar.gz
cd ../output
wget http://tootsuite.encs.concordia.ca/TranCEP/output.tar.gz
tar -xzf output.tar.gz
rm -f output.tar.gz
cd ..
cp src/* .
chmod 755 TranCEPTool.R
wget http://tootsuite.encs.concordia.ca/singularity-images/bioinformatics-singularity.simg
singularity -v exec -B .:/TranCEP --pwd /TranCEP bioinformatics-singularity.simg ./TranCEPTool.R -query=testc.fasta -trancepdir=/TranCEP -out=TranCEPoutput
```

The results will be located in `TranCEPoutput/TranCEPout.csv`


