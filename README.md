The original GEMF implementation in C was developed by Futing Fan, [Copyright (c) 2016](LICENSE). It has since been updated by [Niema Moshiri](https://niema.net/).

This version is further modified to efficiently and deterministically reproduce the analysis of [The molecular epidemiology of multiple zoonotic origins of SARS-CoV-2](https://www.science.org/doi/10.1126/science.abp8337). There are no corrections or improvements to the analyses, but the more efficient implementation allows precise measurement of the simulation properties without undue computational cost.

Note that the Bayes factors are a convex function of the binomially distributed topology frequencies, so increasing the sample size will slightly reduce the expected values of the Bayes factors (Jensen's inequality).

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
git clone --branch reproduction https://github.com/nizzaneela/GEMF.git
cd GEMF
sudo mv GEMF /usr/local/bin/
make
```

The analyses can then be reproduced with python script `run.py`. It takes three arguments: the name of a file containing the input parameters for `GEMF`, the number of processes to run in parallel, and the total number of successful epidemics to simulate and analyze, e.g.:

```bash
python3 run.py main.txt 20 1100
```

Parameter files for each of the five simulation cases are included. The analyses are deterministic, seeded with a hash of the parameter files.

Suggestions for improvements are welcome. nizzaneela@fastmail.com.
