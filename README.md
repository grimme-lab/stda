# *std2* program for computing excited states and response functions via simplified TD-DFT methods (sTDA, sTD-DFT, SF-sTD-DFT, XsTDA, XsTD-DFT, and SF-Xs-TD-DFT)[![DOI](https://zenodo.org/badge/221426808.svg)](https://doi.org/10.5281/zenodo.4022460)

This project provides the `std2` program.

The `std2` program is the rebranded and updated version of the `stda` program. Originally, `stda` was implemented only for the simplified time-dependent density functional theory using the Tamm-Dancoff approximation (sTDA) method. With the implementation of more simplified quantum chemistry (sQC) methods in `stda`, the name was not fitting the application of the program anymore.

## Installation

Two options exist: `make` or `meson`.

### Using `make` (and the intel compiler)

For that option, you need:

+ [`cmake`](https://cmake.org/),
+ the latest intel oneAPI Fortran compiler, `ifx` (**not** `ifort`), with `MKL`.

Then, you need to download the library to compute one- and two-electron integrals : `libcint`:

```bash
# create a `libcint` directory, then download sources in it
mkdir libcint
cd libcint
wget https://github.com/pierre-24/libcint-meson/releases/download/v0.3.0/libcint_v6.1.2.tar.gz -O libcint.tar.gz
tar -xzf libcint.tar.gz
```

The next step depends on your target.
The default (32 bits integers, `LP64`) use a bit less memory but limits the size of the system you can treat.
If you target large system, use the 64 bit integers (`ILP64`) instead.

#### 32 bit integers (default, `LP64`)

First, build `libcint`:

```bash
# use intel
export CC=icx

# In the libcint directory, create a build directory and build `libcint` in it (using cmake)
mkdir build
cd build
cmake ..
cmake --build .
```

Then, compile `std2` itself:

```bash
# go back
cd ../..

# make std2
make
```

#### 64 bit integers (`ILP64`) for larger calculations

First, build `libcint` with the option for 64 bits integers:

```bash
export CC=icx

# In the libcint directory, create a build directory and build `libcint` in it (using cmake)
mkdir build
cd build
cmake .. -DI8=true
cmake --build .
```

Then, compile `std2` itself:

```bash
# go back
cd ../..

# make std2 (using ILP64)
make USEILP64=1
```
Troubleshootings:

In some cases, the path to libraries is a bit different and the Makefile should be adapted as

```
LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5
```

#### Run

You will find a executable named `std2` in this folder.
To make `std2` accessible, do:

```bash
export STD2HOME=/path/to/std2/folder
export LD_LIBRARY_PATH=$STD2HOME/libcint/build:$LD_LIBRARY_PATH
export PATH=$PATH:$STD2HOME
```

in your `.bashrc` or submission scripts.

### Using `meson` (and any compiler)

**Note that for the moment, it is not possible to use Meson to compile the 64 bit version with `ifx` (see [there](https://github.com/mesonbuild/meson/issues/13052))**

If you are not found of `make`, you can use [`meson`](https://mesonbuild.com/) and [`ninja`](https://ninja-build.org/) instead.
Other advantages include: automatic `libcint` import, more flexibility on the linear algebra backend and `gfortran` instead of intel.

First of all, if you want to use other compilers than `gfortran`, use:

```bash
# for latest version of intel compilers
export FC=ifx CC=icx

# for older versions of intel compilers
export FC=ifort CC=icc
```

Then, pick one of the `meson setup` line below:

```bash
# netlib BLAS and LAPACK, 32 bits integers
meson setup _build -Dla_backend=netlib

# openblas and netlib LAPACK, 32 bits integers
meson setup _build -Dla_backend=openblas

# MKL, 32 bits integers
meson setup _build -Dla_backend=mkl

# MKL, 64 bits integers (ILP64)
meson setup _build -Dla_backend=mkl -Dinterface=64
```

You can also:

+ generate a statically linked executable by adding `-Dstatic=true` (only with MKL), or
+ disable OpenMP (not recommended) by adding `-Dopenmp=false`.

And finally, compile everything with

```bash
meson compile -C _build
```

#### Run

You will an executable named `std2` in the `_build` directory.
To make `std2` accessible, export

```bash
export STD2HOME=/path/to/std2/folder
export PATH=$PATH:$STD2HOME/_build/
```

in your `.bashrc` or submission scripts.

## Usage

For parallel usage set the threads for OMP and the MKL linear algebra backend by

```bash
export OMP_NUM_THREADS=<ncores>
```

For larger systems please adjust the stack size accordingly, otherwise
stack overflows *will* occur. Use something along the lines of this:

```bash
ulimit -s unlimited
export OMP_STACKSIZE=4G
```

See the manual on the [release page](https://github.com/grimme-lab/stda/releases/latest).

## Citations

- S. Grimme, A simplified Tamm–Dancoff density functional approach for the electronic excitation spectra of very large molecules, *J. Chem. Phys.*, **2013**, 138, 244104.
  DOI: [10.1063/1.4811331](https://doi.org/10.1063/1.4811331)

- C. Bannwarth, S. Grimme, A simplified time-dependent density functional theory approach for electronic ultraviolet and circular dichroism spectra of very large molecules, *Comput. Theor. Chem.*, **2014**, 1040 – 1041, 45 – 53.
  DOI: [10.1016/j.comptc.2014.02.023](https://doi.org/10.1016/j.comptc.2014.02.023)

- S. Grimme and C. Bannwarth,  Ultra-fast computation of electronic spectra for large systems by tight-binding based simplified Tamm-Dancoff approximation (sTDA-xTB) *J. Chem. Phys.*, **2016**, 145, 054103.
  DOI: [10.1063/1.4959605](https://dx.doi.org/10.1063/1.4959605)

- M. de Wergifosse, S. Grimme, Nonlinear-response properties in a simplified time-dependent density functional theory (sTD-DFT) framework: Evaluation of the first hyperpolarizability, *J. Chem. Phys.*, **2018**, 149 (2), 024108.
  DOI: [10.1063/1.5037665](https://doi.org/10.1063/1.5037665)

- M. de Wergifosse, S. Grimme, Nonlinear-response properties in a simplified time-dependent density functional theory (sTD-DFT) framework: Evaluation of excited-state absorption spectra, *J. Chem. Phys.*, **2019**, 150, 094112.
  DOI: [10.1063/1.5080199](https://doi.org/10.1063/1.5080199)

- M. de Wergifosse, C. Bannwarth, S. Grimme, A simplified spin-flip time-dependent density functional theory (SF-sTD-DFT) approach for the electronic excitation spectra of very large diradicals, *J. Phys. Chem. A*, **2019**, 123 (27), 815–5825.
  DOI: [10.1021/acs.jpca.9b03176](https://doi.org/10.1021/acs.jpca.9b03176)

- M. de Wergifosse, J. Seibert, B. Champagne, and S. Grimme, Are fully conjugated expanded indenofluorenes analogues and diindeno[n]thiophene derivatives diradicals? A simplified (spin-flip) time-dependent density functional theory [(SF-)sTD-DFT] study, *J. Phys. Chem. A*, **2019**, 123 (45), 9828-9839.
  DOI: [DOI: 10.1021/acs.jpca.9b08474](https://doi.org/10.1021/acs.jpca.9b08474)

- M. de Wergifosse, J. Seibert, S. Grimme, Simplified time-dependent density functional theory (sTD-DFT) for molecular optical rotation, *J. Chem. Phys.*, **2020**, 153, 084116.
  DOI: [10.1063/5.0020543](https://doi.org/10.1063/5.0020543)

- M. de Wergifosse, S. Grimme, A unified strategy for the chemically intuitive interpretation of molecular optical response properties, *J. Chem. Theory Comput.*, **2020**, 16 (12), 7709–7720.
  DOI: [10.1021/acs.jctc.0c00990](https://doi.org/10.1021/acs.jctc.0c00990)

- M. de Wergifosse, S. Grimme, Perspective on simplified quantum chemistry methods for excited states and response properties, *J. Phys. Chem. A*, **2021**, *J. Phys. Chem. A*, **2021**, 125 (18) 3841–3851.
  DOI: [10.1021/acs.jpca.1c02362](https://doi.org/10.1021/acs.jpca.1c02362)

- P. Beaujean, B. Champagne, S. Grimme, and M. de Wergifosse, All-atom quantum mechanical calculation of the second-harmonic generation of fluorescent proteins, *J. Phys. Chem. Lett.*, **2021**, 12 (39), 9684-9690.
  DOI: [10.1021/acs.jpclett.1c02911](https://doi.org/10.1021/acs.jpclett.1c02911)

- M. de Wergifosse, P. Beaujean, S. Grimme, Ultrafast evaluation of two-photon absorption with simplified time-dependent density functional theory, *J. Phys. Chem. A*, **2022**, 126 (41) 7534–7547.
  DOI: [10.1021/acs.jpca.2c02395](https://doi.org/10.1021/acs.jpca.2c02395)

- S. Löffelsender, P. Beaujean, M. de Wergifosse. Simplified quantum chemistry methods to evaluate non-linear optical properties of large systems, *WIREs Comput Mol Sci.* **2024**, 14 (1) e1695.
  DOI: [10.1002/wcms.1695](https://doi.org/10.1002/wcms.1695)

- M. de Wergifosse, S. Grimme, The eXact integral simplified time-dependent density functional theory (XsTD-DFT), *J. Chem. Phys.*, **2024**, 160, 204110.
  DOI: [10.1063/5.0206380](https://doi.org/10.1063/5.0206380)

- M. de Wergifosse, Computing excited states of very large systems with range-separated hybrid functionals and the eXact integral simplified time-dependent density functional theory (XsTD-DFT), *J. Phys. Chem. Lett.*, **2024**, 15, (51) 12628–12635.
  DOI: [10.1021/acs.jpclett.4c03193](https://doi.org/10.1021/acs.jpclett.4c03193)

## License

`std2` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

`std2` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.

## Bugs

A bug is a *demonstratable problem* caused by the code in this repository.
Good bug reports are extremely valuable for us - thank you!

Before opening a bug report:

1. Check if the issue has already been reported.
2. Check if it still is an issue or has already been fixed?
   Try to reproduce it with the latest version from the `main` branch.
3. Isolate the problem and create a reduced test case.

A good bug report should not leave others needing to chase you up for more
information. So please try to be as detailed as possible in your report,
answer at least these questions:

1. Which version of `std2` are you using? The current version is always
   a subject to change, so be more specific.
   If possible, also provide the *commit*.
2. What is your environment (your laptop, the cluster of the university)?
3. What steps will reproduce the issue?
   We have to reproduce the issue, so we need all the input files.
4. What would be the expected outcome?
5. What did you see instead?

All these details will help people to fix any potential bugs.
