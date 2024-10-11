The original GEMF implementation in C was developed by Futing Fan, [Copyright (c) 2016](LICENSE). It has since been updated by [Niema Moshiri](https://niema.net/).

This version will be further modified to efficiently and deterministically reproduce the analysis of [The molecular epidemiology of multiple zoonotic origins of SARS-CoV-2](https://www.science.org/doi/10.1126/science.abp8337), based on the code published in the authors' GitHub [repository](https://github.com/sars-cov-2-origins/multi-introduction) and the methods described in the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf). 

## Installation
The `GEMF` tool is written in C and has no dependencies beyond a modern C++ compiler (and `make` for convenience). You can simply download the latest release tarball (or clone the repo) and compile with `make`:

```bash
git clone https://github.com/niemasd/GEMF.git
cd GEMF
make
sudo mv GEMF /usr/local/bin/ # optional step to install globally
```
