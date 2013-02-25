#-------------------------------------------------------------
#
#	Sequence search in TreeFam database.
#
#
#-------------------------------------------------------------

#First, clone the TreeFam github:
git clone https://github.com/treefam/treefam_tools.git
cd treefam_tools/

#Second, you need to download the TreeFam HMM file:
wget http://www.treefam.org/download/TreeFam.hmm.tar.gz hmm_lib/
cd hmm_lib/
tar -xzvf TreeFam.hmm.tar.gz

#Third you need to construct the binary compressed datafiles for hmmscan.
hmmpress TreeFam.hmm

# Now you are ready to search your sequences (in FASTA format) against the HMMs library:
perl treefam_scan.pl -fasta example_small.fasta -dir hmm_lib/ > output.tab

#Help:
perl treefam_scan.pl -h
