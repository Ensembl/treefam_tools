
# Sequence search in TreeFam database.

> Note that we have tested this script with HMMer version 3.0. It may not
> work with version 3.1


## First, clone the TreeFam github:

```
git clone https://github.com/treefam/treefam_tools.git
cd treefam_tools/treefam_scan/
```

## Second, you need to download the TreeFam HMM file:

```
mkdir hmm_lib
cd hmm_lib/
wget http://www.treefam.org/static/download/treefam9.hmm3.tar.gz
tar -xzvf treefam9.hmm3.tar.gz
```

## Third you need to construct the binary compressed datafiles for hmmscan.

```
mv treefam9.hmm3 TreeFam9
hmmpress TreeFam9
```

## Run the script on your input Fasta file

Back to `treefam_tools/treefam_scan/`, run:
```
perl treefam_scan.pl -fasta example/example_small.fasta -dir hmm_lib/ -hmm_file TreeFam9 > example/output.tab
```

> Use the same name in `-hmm_file` as in the above `hmmpress` command

## Help:

```
perl treefam_scan.pl -h
```

