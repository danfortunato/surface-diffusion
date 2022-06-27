# Spectral methods for reaction-diffusion on axisymmetric surfaces

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6762738.svg)](https://doi.org/10.5281/zenodo.6762738)

Requirements:
* [Chebfun](http://www.chebfun.org/) for conveniently computing with smooth functions in MATLAB.
* [This repository](https://github.com/danfortunato/spherical-harmonic-interfaces) for spherical harmonic transforms (only required for the sphere).

To build the batched banded linear system solver located in the `@BandedBatch` directory, link `make.inc` to the appropriate template for your operating system (or modify it accordingly) and run `make`.
