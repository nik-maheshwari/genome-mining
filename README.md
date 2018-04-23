# genome-mining
This repository consists of a series of Perl scripts to mine for lanthipeptide gene clusters from bacterial genome.
Please download as a Zip file and follow the instructions below.

Please read completely.

Dependencies:

The scripts are written in Perl v5.18.2 and BioPerl 1.6.923, tested on Ubuntu 14.04 LTS and Ubuntu 16.04 LTS.
Other Perl libraries used are List::Util, File::chdir and Cwd.

HMMer v3.1b2 installation is required.
Prodigal v2.6.2 executable is required (download the executable at "https://github.com/hyattpd/Prodigal/releases"  and save in the same folder as the scripts and genomic files)
GGSEARCH36 v36.3.7a executable is required (download the executable at "https://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml" and save in the same folder as the scripts and genomic files)

Steps for predicting lanthipeptide gene clusters (LGCs)

1. Unzip and copy all the files to a new folder (say, f1)

2. Move your genomic files, with .fasta extension, to folder f1. 
   The filename should be in the format: strain_file1.fasta
	To allow for multiple fasta file for the same strain, have strain_file1.fasta, strain_file2.fasta, etc. in folder f1.
There should be no _ or spaces in strain name.

Ex: ACCEPTED FILENAMES - StreptococcusPneumoniae_file1.fasta, StreptococcusPneumoniae_file2.fasta, StreptococcusPneumoniae_file1.fasta, nisin_file1.fasta

Ex: NOT ACCEPTED FILENAMES - Streptococcus Pneumoniae_file1.fasta, Streptococcus_Pneumoniae_file1.fasta, Streptococcus.Pneumoniae_file1.fasta

The genomic files should have a NCBI header, with GI number (ex. >gi|27763998|emb|AJ536588.1| Streptomyces cinnamoneus cinnamycin biosynthetic cluster) (Go to NCBI->search AJ536588.1-> Send complete record to file, in fasta format and show GI)

A few test files are included (cinnamycin_dna.fasta and gallidermin_dna.fasta)

3. Do not have any important folder or .txt file in folder f1, other than your .fasta files and scripts from step 1.

4. Open script "addBridge_lanCDeh.pl" in folder f1. At lines 8 and 10, change path to full path of "f1".

5. Open a terminal, change directory to folder f1.

6. run the following command:
	perl run.pl

7. Any predicted clusters will be in the following file in folder f1:	clusters_ggsearch_bridge_lanCDeh.txt.
   Individual files will be in the folder "strain" in folder f1

Thank you for using this pipeline.

Note:

1. To predict the bridging pattern for any peptide sequence, use the following script:
	
	addBridge.pl
