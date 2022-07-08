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

Driver scripts are provided for surfaces of revolution (`sor/SORDiffusionDriver.m`) as well as for the sphere (`sphere/SphericalDiffusionDriver.m`) and the torus (`torus/TorusDiffusionDriver.m`). To choose a surface of revolution in `sor/SORDiffusionDriver.m`, change the variable `shape` to one of the geometries listed (e.g., `shape = 'oblate';`). To run the driver script in MATLAB:
```
>> cd surface-diffusion/sor
>> SORDiffusionDriver
```

Geometry parameters, biological parameters, discretization parameters, and visualization parameters may be modified in the driver scripts. A brief summary of the parameters is given below.

### Geometry parameters

<table>
  <tr>
    <td><code>params.rho</code></td>
    <td>Radial spherical coordinate of generating curve</td>
  </tr>
  <tr>
    <td><code>params.theta</code></td>
    <td>Angular spherical coordinate of generating curve</td>
  </tr>
</table>

### Biological parameters

<table>
  <tr>
    <td><code>params.delta</code></td>
    <td>Diffusivity</td>
  </tr>
  <tr>
    <td><code>params.alpha</code></td>
    <td>Coupling parameter</td>
  </tr>
  <tr>
    <td><code>params.beta</code></td>
    <td>Constitutive rate</td>
  </tr>
  <tr>
    <td><code>params.gamma</code></td>
    <td>Membrane recruitment threshold</td>
  </tr>
  <tr>
    <td><code>params.nu</code></td>
    <td>Cooperativity parameter</td>
  </tr>
</table>

### Discretization parameters

<table>
  <tr>
    <td><code>params.nlat</code></td>
    <td>Number of latitudinal coefficients</td>
  </tr>
  <tr>
    <td><code>params.nlon</code></td>
    <td>Number of longitudinal coefficients</td>
  </tr>
  <tr>
    <td><code>params.dt</code></td>
    <td>Time step</td>
  </tr>
  <tr>
    <td><code>params.tend</code></td>
    <td>End time</td>
  </tr>
  <tr>
    <td><code>params.scheme</code></td>
    <td>Time-stepping scheme (<code>'bdf1'</code>, <code>'bdf2'</code>, or <code>'bdf4'</code>)</td>
  </tr>
  <tr>
    <td><code>params.useHeaviside</code></td>
    <td>Replace the nonlinearity with a Heaviside</td>
  </tr>
  <tr>
    <td><code>params.trueConservation</code></td>
    <td>Enforce conservation according to the true volume</td>
  </tr>
</table>

### Visualization parameters

<table>
  <tr>
    <td><code>params.quiet</code></td>
    <td>Flag to keep quiet</td>
  </tr>
  <tr>
    <td><code>params.movie</code></td>
    <td>Flag to plot a movie during simulation</td>
  </tr>
  <tr>
    <td><code>params.colormap</code></td>
    <td>Colormap for plotting</td>
  </tr>
  <tr>
    <td><code>params.movfile</code></td>
    <td>Filename of movie to output</td>
  </tr>
  <tr>
    <td><code>params.keepAll</code></td>
    <td>Flag to return an array of solutions at all time points</td>
  </tr>
</table>

## Replicating results from the paper

Geometries studied in the paper are included as preset options at the beginning of `sor/SORDiffusionDriver.m` and most results from the paper can be easily replicated by choosing the appropriate shape and parameters consistent with those given in the paper.

Altering the PDE requires direct modification of `sor/SORDiffusion.m` (or equivalent files for the spherical or toroidal cases). To change the chemical reactions present (as in the Supplementary Material) one must define a new nonlinear function; see `N_direct` and `N_heaviside` in `sor/SORDiffusion.m` for examples.

For analysis of interface length, we have included a function `+util/interface_length.m` which calculates the arc length of a chosen level set `gamma` of the given coefficient matrix `U` and position matrices `xx`, `yy`, and `zz`, as defined in `sor/SORDiffusionDriver.m`. The analysis in Figure 2 of the paper can be replicated by enabling `params.keepAll` and then calling `util.interface_length` sequentially on `U(:,:,k)` for each time point `k`.
