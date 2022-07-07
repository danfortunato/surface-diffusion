# Spectral methods for reaction-diffusion on axisymmetric surfaces

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6762738.svg)](https://doi.org/10.5281/zenodo.6762738)

This repository contains source code accompanying the paper "Forced and spontaneous symmetry breaking in cell polarization" by P. Miller, D. Fortunato, C. Muratov, L. Greengard, and S. Shvartsman.

The code is structured as follows:
* `sor/` contains a driver script and source code for surfaces of revolution, which generalizes the classical double Fourier sphere method to axisymmetric surfaces.
* `sphere/` contains a driver script and source code for the sphere, which represents functions using spherical harmonics.
* `torus` contains a driver script and source code for the torus, which represents functions using bivariate Fourier series.
* `+util/` contains utility functions.
* `@BandedBatch/` contains a class for solving a batch of banded linear systems in parallel.

## Requirements

* [Chebfun](http://www.chebfun.org/) for conveniently computing with smooth functions in MATLAB.
* [This repository](https://github.com/danfortunato/spherical-harmonic-interfaces) for spherical harmonic transforms (optional, only required for the sphere).

For download and installation instructions for these repositories, consult their respective documentations. Once installed, make sure they are added to the MATLAB path.

## Installation

First, clone this repository. In a terminal, run:
```
git clone https://github.com/danfortunato/surface-diffusion.git
```
To build the MEX files for the batched banded linear system solver located in the `@BandedBatch` directory, link `make.inc` to the appropriate template file for your operating system (or modify it accordingly) and run `make`. For example:
```
cd surface-diffusion
ln -s make.inc.mac make.inc
make
```
Finally, add this repository to the MATLAB path. In MATLAB, run:
```
>> addpath surface-diffusion
```

## Driver scripts

Driver scripts are provided for surfaces of revolution (`sor/SORDiffusionDriver.m`) as well as for the sphere (`sphere/SphericalDiffusionDriver.m`) and the torus (`torus/TorusDiffusionDriver.m`). Biological parameters, discretization parameters, and visualization parameters may be modified in the driver scripts. To choose a surface of revolution in `sor/SORDiffusionDriver.m`, change the variable `shape` to one of the geometries listed (e.g., `shape = 'oblate';`). To run the driver script in MATLAB:
```
>> cd surface-diffusion/sor
>> SORDiffusionDriver
```
