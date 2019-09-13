This tool predicts the substrate class of a given transporter. The class can belong to one of the following seven categories (ordered alphabetically)
[1] "amino"    [2]”anion"    [3]”cation"   [4]”electron" [5]”other"    [6]”protein" 
[7] "sugar" 
Input: transporter proteins sequences in Fasta format
Output: the substrate class


&&&&&&&& FOLDERS &&&&&&&& :
Dataset:
Contains both the training and the independent testing dataset from the benchmark dataset from TrSSP server available at http://bioinfo. noble.org/TrSSP


Compositions:
Contains the extracted MSA_TCS_PAAC features of each protein in the training and independent datasets. The features are saved in csv format, the first column indicates the unitportID, the second column contains the substrate class (encoded as numbers from 1:7 according to their alphabetical order) and the next 400 columns contains the features.


output:
Contains the homology details needed to extract the features. Details of the MSA, TOPCONS scores, Blast hits for each sequence is found here.

src:
The scripts needed to use the tool.

&&&&&&&& HOW TO USE &&&&&&&& :
*This tool requires TM-COFFEE and BLAST pre-installed'
* Usage: TranCEP -query=<input> [-trancepdir=<trancepdir>] [-out=<outdir>]
	<input> is your sequence input file in fasta format
	<out> is the output directory where you want the predicted 	results, formatted as csv
	<trancepdir> is the directory where the base TranCEP files 	are located
*MSA_TCS_PAAC features of each sequence in the test set is  found under TranCEP/Compositions/MSA_TCS_PAAC.csv



