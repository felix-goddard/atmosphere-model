# Atmosphere model

This codebase implements a parallelised Fortran model of the shallow water equations. The eventual plan is to produce a 3D geophysical fluid simulation -- but we're a long way away from that still!

It was originally built off the code provided as companion to [Modern Fortran: Building Efficient Parallel Applications](https://github.com/modern-fortran/tsunami). At this point, that code forms the basis of the logic for dividing work between cores but the main control and solution routines have been restructured. The method for solving the equations is based on the work of Shian-Jiann Lin and Richard B. Rood, namely Lin & Rood (1996) and Lin & Rood (1997).

# Building

The code has OpenCoarrays and netCDF-Fortran as dependencies, which on my machine were installed via Homebrew so the Makefile reflects that. If your setup is different you will need to edit the Makefile so the linker can find those libraries.

Otherwise, the code can be built simply by running
```bash
make
```
from the base directory of the file tree. The resulting binary, `model`, can then be run with
```shell
$ cafrun -n [cores] ./model
```
where `[cores]` should be replaced with the number of cores to run the script on in parallel, e.g. `1` for serial execution.

# References

Lin, S. J., & Rood, R. B. (1996). Multidimensional flux-form semi-Lagrangian transport schemes. *Monthly weather review, 124*(9), 2046-2070.

Lin, S. J., & Rood, R. B. (1997). An explicit flux‐form semi‐Lagrangian shallow‐water model on the sphere. *Quarterly Journal of the Royal Meteorological Society, 123*(544), 2477-2498.