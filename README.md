The original GEMF implementation in C was developed by Futing Fan, [Copyright (c) 2016](LICENSE). It has since been updated by [Niema Moshiri](https://niema.net/).

This version is further modified to efficiently and deterministically replicate the analysis of [The molecular epidemiology of multiple zoonotic origins of SARS-CoV-2](https://www.science.org/doi/10.1126/science.abp8337). The only corrections to the single intro model are the starting the primary case as exposed and omitting the pruning of  short lived basal lineage with the stable coalescence. Two intro model draws pairs of single intro epidemics, offsets all events in each by t0 and t1, respectively, combines them with an MRCA and upstream lineages of length t0 and t1, respectively, and then puts the combined epidemic data through the same pipeline as for the single intro model. It's repeated for t0 and t1 in {0, 5, 10, 15, 20, 25, 30} days. 

## Installation

The analysis uses `CoaTran`, which should be downloaded, compiled and installed to your `PATH`:

```bash
git clone https://github.com/niemasd/CoaTran.git
cd CoaTran
make
sudo mv coatran_* /usr/local/bin/
cd ..
```

This repository should then be cloned, and the `GEMF` tool compiled and installed:

```bash
git clone --branch replication https://github.com/nizzaneela/GEMF.git
cd GEMF
sudo mv GEMF /usr/local/bin/
make
```

The analysis can then be reproduced with python script `replicate.py`. It takes one argument: the number of processes to run in parallel, e.g.:

```bash
python3 replicate.py 20
```
Each process needs about 3.5Gb RAM.

The code generates 110,000 simulated epidemics, from which it draws 110,000 pairs of epidemics, and subjects the simulated epidemics and the pairs of simulated epidemics to the same analysis. The results are stored in `one_intro/n` and `two_intro/n`, where n is 0 to 109,999.

Results for any particular simulated epidemic n in the sample 0 to 109,999 can be verified with `verify_one_intro.py $n$`.

Results for any particular pair of simulated epidemics n in the sample 0 to 109,999 can be verified with `verify_two_intro.py $n$`.

The final stochastic product of each simulation and pair of simulations is hashed to `tree_final_hash.txt`, which can be used for verification. E.g., if you want to verify each of the reported CC topologies, search the `one_intro` directory to find the numbers n of the simulations where the content of `CC.txt` is non-zero (there's probably about 300), and then run `verify_one_intro.py n` for each n, and compare the hashes. The idea is to provide verifiable output in a small list of hashes, rather than terabytes of data. Please let me know if you have a better idea for doing this. 


Suggestions for improvements are welcome. nizzaneela@fastmail.com.
