# Construction of Laguerre pseudospectral differentiation matrices

This repository contains MATLAB code used to generate the figures for the paper
*Construction of Laguerre pseudospectral differentiation matrices*
by E. Nel and N. Hale.

The main contribution is the routine `lagdif_stable.m` which generates first- and second-order differentiation matrices for Laguerre spectral collocation for very large degrees n.

## Dependencies

The code relies on the external library [DMSUITE](https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite) (`lagdif.m`, `poldif.m`, `lagroots.m`). \
The function `lagpts_new.m` in the `src/` directory was adapted from [Chebfun](https://www.chebfun.org/).

## Repository Structure
- `src/`          : core routines
- `figures/`      : scripts to reproduce figures from the paper
- `data/`         : high-precision reference data
- `output/`       : generated figures

## How to Run

Apart from [DMSUITE](https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite), all required functions are located in the `src/` folder. Some scripts rely on precomputed high-precision data stored in the `data/` folder. 

To reproduce a figure, navigate to the `figures/` folder in MATLAB and run, for example: fig1_convergence.
