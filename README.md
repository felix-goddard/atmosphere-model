# Atmosphere model

This codebase implements a simple atmospheric model written in Fortran for execution on multiple cores. The specification of the model can be summed up as:
- a vertically-Lagrangian, finite volume dynamical core solving the compressible hydrostatic primitive equations;
- a two-stream radiation code, using the correlated *k*-distribution method for gas absorption;
- assorted bits and bobs for physics (currently just a simple vertical mixing for dry convection).

The model was originally built off the code provided as companion to [Modern Fortran: Building Efficient Parallel Applications](https://github.com/modern-fortran/tsunami). At this point, that code forms the basis of the logic for dividing work between cores but the main control and solution routines have been pretty fundamentally restructured.

The method for solving the dynamics equations is based on the work of Shian-Jiann Lin and Richard B. Rood, namely Lin & Rood (1996) and Lin & Rood (1997). I also relied heavily upon the technical documentation [and code](https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/tree/main) for the FV3 dynamical core (Harris, 2021), which is based on the same work. In writing the radiation code I relied heavily upon Meador & Weaver (1980), Shonk & Hogan (2008), and Petty (2006). The gas optics data used in the radiation code is derived from [the CKD data provided by the ECMWF](https://confluence.ecmwf.int/display/ECRAD/ECMWF+gas+optics+tool%3A+ecCKD), calculated using the ecCKD tool (Hogan & Matricardi, 2022).

# Building

The code has OpenCoarrays and netCDF-Fortran as dependencies, which on my machine were installed via Homebrew so the Makefile reflects that. If your setup is different you will need to edit the Makefile so the linker can find those libraries.

Otherwise, the code can be built simply by running
```bash
$ make
```
from the base directory of the file tree. The resulting binary, `model`, can then be run with
```shell
$ cafrun -n [cores] ./model
```
where `[cores]` should be replaced with the number of cores to run the script on in parallel, e.g. `1` for serial execution.

# References

Harris, L. (2021). A Scientific Description of the GFDL Finite-Volume Cubed-Sphere Dynamical Core. https://doi.org/10.25923/6NHS-5897

Hogan, R. J., & Matricardi, M. (2022). A Tool for Generating Fast k-Distribution Gas-Optics Models for Weather and Climate Applications. Journal of Advances in Modeling Earth Systems, 14(10), e2022MS003033. https://doi.org/10.1029/2022MS003033

Lin, S., & Rood, R. B. (1997). An explicit flux‐form semi‐lagrangian shallow‐water model on the sphere. Quarterly Journal of the Royal Meteorological Society, 123(544), 2477–2498. https://doi.org/10.1002/qj.49712354416

Lin, S.-J., & Rood, R. B. (1996). Multidimensional Flux-Form Semi-Lagrangian Transport Schemes. Monthly Weather Review, 124(9), 2046–2070. https://doi.org/10.1175/1520-0493(1996)124<2046:MFFSLT>2.0.CO;2

Meador, W. E., & Weaver, W. R. (1980). Two-Stream Approximations to Radiative Transfer in Planetary Atmospheres: A Unified Description of Existing Methods and a New Improvement. https://journals.ametsoc.org/view/journals/atsc/37/3/1520-0469_1980_037_0630_tsatrt_2_0_co_2.xml

Petty, G. W. (2006). A first course in atmospheric radiation (2. ed). Sundog Publ.

Shonk, J. K. P., & Hogan, R. J. (2008). Tripleclouds: An Efficient Method for Representing Horizontal Cloud Inhomogeneity in 1D Radiation Schemes by Using Three Regions at Each Height. https://doi.org/10.1175/2007JCLI1940.1
