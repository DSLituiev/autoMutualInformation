autoMutualInformation
=====================

a MATLAB wrapper script for a mex program to calculate auto-mutual information of a time series. This method is a non-linear alternative to auto-correlation function.

## Files
+ `calcMutInfDelay.m` -- a MATLAB wrapper script
+ `mexCalcMutInfDelay2D.c` -- a MEX/C program 

Calculates mutual information of the delayed time series for each delay. 


## Input
The script takes one to three arguments:

 - `X` -- an N-dimensional array of time series,
 with the time running along the first dimension <br />
          i.e. array size `[ T x S1 x S2 x ... x SN ]`,<br />
          where `T` is the number of time points
- `tauMax` (optional) -- the maximal time delay;
(default: 20)
- `bins` (optional) -- number of bins

## Output
- `Z` -- an array of the mutual information; <br />
the first dimension of `Z` is equal to (`tauMax+1`) <br />
with the corresponding timeline (`0:1:tauMax`) <br />
Uncomputeable time points are replaced with `NaN`
- `t` (optional) -- the lag-timeline: `t = (0:1:tauMax)'`

## References
1. [Andrew M. Fraser and Harry L. Swinney. "Independent coordinates for strange attractors from mutual information", Phys. Rev. A 33 (1986) 1134-1140.]( http://dx.doi.org/10.1103%2fPhysRevA.33.1134 )
2. [Eric Weeks' page](http://www.physics.emory.edu/~weeks/software/minfo.html) with the description of the modified algorithm


===============================
Author: Dmytro S. Lituiev (2013), University of Zurich,
based on the Eric Weeks' C script (1997)

