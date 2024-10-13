The original GEMF implementation in C was developed by Futing Fan, [Copyright (c) 2016](LICENSE). It has since been updated by [Niema Moshiri](https://niema.net/).

This version will be further modified to efficiently and deterministically reproduce the analysis of [The molecular epidemiology of multiple zoonotic origins of SARS-CoV-2](https://www.science.org/doi/10.1126/science.abp8337). 

## Installation
The `GEMF` tool is written in C and has no dependencies beyond a modern C++ compiler (and `make` for convenience). You can simply clone the repo from this branch and compile with `make`:

```bash
git clone --branch reproduction https://github.com/nizzaneela/GEMF.git
cd GEMF
make
```

`GEMF` can then be run with the input file `para.txt`.

```bash
./GEMF para.txt
```
