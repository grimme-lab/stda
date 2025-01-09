# *stda* program for computing excited states and response functions via simplified TD-DFT methods (sTDA, sTD-DFT, and SF-sTD-DFT)[![DOI](https://zenodo.org/badge/221426808.svg)](https://doi.org/10.5281/zenodo.4022460)

This project provides the `stda` program.

## Installation

Statically linked binaries can be found at the projects
[release page](https://github.com/grimme-lab/stda/releases/latest).
To build from source this project uses a `make` based build system and requires
a version of Intel Parallel Studio 17 or newer to be compiled.
To trigger the build run in the root directory

```bash
make
```

You will find a statically linked executable named `stda`.
To make `stda` accessible export

```bash
export STDAHOME=$PWD
export PATH=$PATH:$STDAHOME
```

### Alternatives

If you are not a fan of `make`, you can use [`meson`](https://mesonbuild.com/)
as alternative, but it requires a fairly new version like 0.49 or newer for a
decent Fortran support.
For the default backend [`ninja`](https://ninja-build.org/) version 1.5 or newer
has to be provided.

To perform a build run:

```bash
export FC=ifort
meson setup build_intel
ninja -C build_intel
```

You may also consider to add the following option to `meson`:

+ `-Dinterface=64` to use 64 bit integers,
+ `-Dstatic=false` to use dynamic libraries (usefull if you get error such as `ld: cannot find -lpthread` at the linking stage of `ninja`),
+ `-Dfortran_link_args=-qopenmp`(useful if you use OneAPI 2021 or latter, and meson ends with `Fortran shared or static library 'mkl_intel_thread' not found`).


To install the `stda` binaries to `/usr/local` use (might require `sudo`)

```bash
ninja -C build_intel install
```

For a local installation (or if you want to pack a release), modify the
configuration by using

```bash
meson configure build_intel --prefix=/
DESTDIR=$HOME/.local ninja -C build_intel install
```

The build system will generate binary files in `$DESTDIR/bin` (don't forget to add that to `PATH` if it is not the case).

## Usage

For parallel usage set the threads for OMP and the MKL linear algebra backend by

```bash
export OMP_NUM_THREADS=<ncores> MKL_NUM_THREADS=<ncores>
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

- M. de Wergifosse, C. Bannwarth, S. Grimme, A simplified spin-flip time-dependent density functional theory (SF-sTD-DFT) approach for the electronic excitation spectra of very large diradicals, *J. Phys. Chem. A*, **2019**, 123 (27), 815–5825.
  DOI: [10.1021/acs.jpca.9b03176](https://doi.org/10.1021/acs.jpca.9b03176)

- M. de Wergifosse, S. Grimme, Nonlinear-response properties in a simplified time-dependent density functional theory (sTD-DFT) framework: Evaluation of the first hyperpolarizability, *J. Chem. Phys.*, **2018**, 149 (2), 024108.
  DOI: [10.1063/1.5037665](https://doi.org/10.1063/1.5037665)

- M. de Wergifosse, S. Grimme, Nonlinear-response properties in a simplified time-dependent density functional theory (sTD-DFT) framework: Evaluation of excited-state absorption spectra, *J. Chem. Phys.*, **2019**, 150, 094112.
  DOI: [10.1063/1.5080199](https://doi.org/10.1063/1.5080199)

- M. de Wergifosse, J. Seibert, S. Grimme, Simplified time-dependent density functional theory (sTD-DFT) for molecular optical rotation, *J. Chem. Phys.*, **2020**, 153, 084116.
  DOI: [10.1063/5.0020543](https://doi.org/10.1063/5.0020543)

- M. de Wergifosse, S. Grimme, A unified strategy for the chemically intuitive interpretation of molecular optical response properties, *J. Chem. Theory Comput.*, **2020**, 16 (12), 7709–7720.
  DOI: [10.1021/acs.jctc.0c00990](https://doi.org/10.1021/acs.jctc.0c00990)

- M. de Wergifosse, P. Beaujean, S. Grimme, Ultrafast evaluation of two-photon absorption with simplified time-dependent density functional theory, *J. Phys. Chem. A*, **2022**, XX, XXXX.
  DOI: [10.1021/acs.jpca.2c02395](https://doi.org/10.1021/acs.jpca.2c02395)

## License

`stda` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

`stda` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.

## Bugs

A bug is a *demonstratable problem* caused by the code in this repository.
Good bug reports are extremely valuable for us - thank you!

Before opening a bug report:

1. Check if the issue has already been reported.
2. Check if it still is an issue or has already been fixed?
   Try to reproduce it with the latest version from the `master` branch.
3. Isolate the problem and create a reduced test case.

A good bug report should not leave others needing to chase you up for more
information. So please try to be as detailed as possible in your report,
answer at least these questions:

1. Which version of `stda` are you using? The current version is always
   a subject to change, so be more specific.
2. What is your environment (your laptop, the cluster of the university)?
3. What steps will reproduce the issue?
   We have to reproduce the issue, so we need all the input files.
4. What would be the expected outcome?
5. What did you see instead?

All these details will help people to fix any potential bugs.
