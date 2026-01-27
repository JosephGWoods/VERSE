% [b1v, gv] = mintverse(b1, g, dt, bmax, gmax, smax, dtout)
%
% MEX function version of minimum-time VERSE algorithm.
%
% in:
%      b1    - RF pulse (real or complex, arb units)
%      g     - gradient pulse (arb units)
%      dt    - time steps of each "block" (arb units)
%      bmax  - maximum RF amplitude, in units of b
%      gmax  - maximum Gradient amplitude in units of g
%      smax  - maximum Gradient slew rate in units of g per unit dt
%      emax  - maximum RF energy (units if b1*b1*dt)
%              (optional, no constraint if omitted)
%
% out:
%      b1v - minimum-time VERSE RF pulse.
%      gv  - minimum-time gradient waveform.
%
% To compile the MEX-code:
%   1. cd to mintverse folder
%   2. type "mex -setup C";
%   3. type "mex('mintverse_mex.c','-R2017b');"
%
% Adapted from Brian Hargreaves' mintverse code
% by Joseph G. Woods, June 2022

function [bv, gv] = mintverse(b1, g, dt, bmax, gmax, smax, emax)

% Do some error checking to prevent matlab from crashing
if nargin < 6;    error('b1, g, dt, bmax, gmax, and smax are required inputs.'); end
if isempty(b1);   error('b1 must not be empty.');   end
if isempty(g);    error('g must not be empty.');    end
if isempty(dt);   error('dt must not be empty.');   end
if isempty(bmax); error('bmax must not be empty.'); end
if isempty(gmax); error('gmax must not be empty.'); end
if isempty(smax); error('smax must not be empty.'); end

if ~exist('emax','var') || isempty(emax)
    emax = -1;
end

[bv, gv] = mintverse_mex(b1, g, dt, bmax, gmax, smax, emax);

end
