% [b1v, gv] = minsarverse(b1, g, dt, gmax, smax)
%
% MEX function version of minimum-SAR VERSE algorithm.
%
% in:
%      b1   - RF pulse (real or complex, arb units).
%      g    - gradient pulse (arb units).
%      dt   - time steps of each "block" (arb units).
%      gmax - maximum gradient amplitude in units of g.
%      smax - maximum gradient slew rate in units of g per unit dt.
%
% out:
%      b1v - minimum-SAR VERSE RF pulse.
%      gv  - minimum-SAR gradient waveform.
%
% To compile the MEX-code:
%   1. cd to minsarverse folder
%   2. type "mex -setup C";
%   3. type "mex('minsarverse_mex.c','-R2017b');"
%
% Adapted from Brian Hargreaves' mintverse code
% by Joseph G. Woods, May 2022

function [b1v, gv] = minsarverse(b1, g, dt, gmax, smax)

% Do some error checking to prevent matlab from crashing
if nargin < 5;    error('b1, g, dt, gmax, and smax are required inputs.'); end
if isempty(b1);   error('b1 must not be empty.');    end
if isempty(g);    error('g must not be empty.');    end
if isempty(dt);   error('dt must not be empty.');   end
if isempty(gmax); error('gmax must not be empty.'); end
if isempty(smax); error('smax must not be empty.'); end

[b1v, gv] = minsarverse_mex(b1, g, dt, gmax, smax);
