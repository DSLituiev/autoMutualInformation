function [Z, varargout] = calcMutInfDelay(X, varargin)
% An nD wrapper for the 'mexCalcMutInfDelay2D.c' function.
% Calculates mutual information of the delayed time series for each delay
% Replaces uncomputed time points with NaN
%
% =========  Input =========
% one to three arguments:
% 
% - X                 -- an nD array of time series,
%                        with the time running along the first dimension
%                        i.e. array size [ T x N x M x ... ], 
%                             where T is the number of time points
% - tauMax (optional) -- the maximal time delay; 
%                         (default: 20)
% - bins   (optional)   -- number of bins
%
% ========= Output =========
% - Z                  -- an array of the mutual information;
%                         the first dimension of Z is equal to (tauMax+1)
%                         with the corresponding timeline    (0:1:tauMax)
% - t    (optional)    -- the lag-timeline:
%                          t = (0:1:tauMax)'
%
% ========= References: =========
%  1. Andrew M. Fraser and Harry L. Swinney. "Independent coordinates for strange attractors from mutual information," Phys. Rev. A 33 (1986) 1134-1140.
%                      http://dx.doi.org/10.1103%2fPhysRevA.33.1134
%  2. Eric Weeks' page with the description of the modified algorithm,
%                      http://www.physics.emory.edu/~weeks/software/minfo.html
%
% ===============================
% Author: Dmytro S. Lituiev  (2013), University of Zurich,
%             based on the Eric Weeks' C script (1997)


%% check the second argument, tauMax
if nargin>1 && isscalar(varargin{1})    
    tauMaxIn = varargin{1};
    tauMaxOut = min(varargin{1}, size(X,1)- 2^6);
else
    tauMaxIn = min(20, size(X,1)-5);
    tauMaxOut = tauMaxIn;
end

dim = size(X);

%== squeeze nD into 2D
X = reshape(X, [dim(1), prod(dim(2:end))] );

%% run the MEX file
Z = mexCalcMutInfDelay2D(X, tauMaxOut, varargin{2:end});
%== tauMax is obtained directly from the MEX function output:
tauMaxOut = size(Z,1)-1;

if tauMaxOut < tauMaxIn
    Z( tauMaxOut+2:tauMaxIn+1, :) = NaN;
end

Z = reshape( Z,  [tauMaxIn + 1, dim(2:end)]);

if nargout>1
    varargout{1} = (0:1:tauMaxIn)';
end
