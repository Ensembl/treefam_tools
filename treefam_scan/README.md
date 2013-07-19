#-------------------------------------------------------------
#
#	Sequence search in TreeFam database.
#
#
#-------------------------------------------------------------

# Note that we have tested this script with HMMer version 3.0


#First, clone the TreeFam github:

git clone https://github.com/treefam/treefam_tools.git

cd treefam_tools/


#Second, you need to download the TreeFam HMM file:

mkdir hmm_lib

wget http://dev.treefam.org/static/download/treefam9.hmm3.tar.gz

mv treefam9.hmm3.tar.gz hmm_lib/

cd hmm_lib/

tar -xzvf treefam9.hmm3.tar.gz


#Third you need to construct the binary compressed datafiles for hmmscan.

mv treefam9.hmm3 TreeFam

hmmpress TreeFam


#Example:

perl treefam_scan.pl -fasta example/example_small.fasta -dir example/hmm_lib/ -hmm_file TreeFam.hmm > example/output.tab


#Help:

perl treefam_scan.pl -h
