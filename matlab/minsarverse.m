%% VERSE RF to minimum SAR version then adjust to meet gradient constraints.
%
% [b1v, Gv] = minsarverse(b1, g, dt, gmax, smax, dmax)
%
% in:
%      b1   - complex B1+ waveform (units, array)
%      g    - gradient waveform (units/distance, array)
%      dt   - b1 and gradient waveform raster or resolution (s, scalar)
%      gmax - maximum gradient amplitude permitted (units/distance, scalar)
%      smax - maximum gradient slew rate permitted (units/distance/s, scalar)
%      dmax - maximum time dilation inverse time compression (unitless, scalar)
%
% out:
%      b1v  - minimum-SAR VERSEd complex B1+ waveform (units, array)
%      gv   - minimum-SAR VERSEd gradient waveform (units/distance, array)
%
% notes:
%   This is a matlab version of minsarverse.c. It is simplified to make use
%   of in-built matlab functionality wherever it leads to more clarity.
%
%   Useful information:
%       1. b1 and g must be the same length with identical dt.
%       2. The gradient can be positive, negative or zero.
%       3. b1 and g amplitudes are defined at the centre of each dt.
%       4. If the RF has steep phase ramps (e.g. far off-isocentre pulse or
%          adiabatic), better results can be achieved if the waveform has
%          small dt. The waveforms can then be downsampled to the required
%          dt after VERSEing.
%
% TODO:
%   1. Add option to low-pass filter the gradient waveform and then adjust
%      the RF correspondingly to avoid hardware issues (Hargreaves used 50 kHz).
%
% Written by Joseph G. Woods, May 2022
%
% Following a similar procedure to algorithm described in Conolly et al.
% JMR 1988 (https://doi.org/10.1016/0022-2364(88)90131-X).
% 
% Adapted from Brian Hargreaves' mintverse algorithm, described in
% Hargreaves et al. MRM 2004 (https://doi.org/10.1002/mrm.20168). The code
% can be found here: http://mrsrl.stanford.edu/~brian/mintverse/
%
% My understanding of the VERSE algorithm in general and implementation of
% the time adjustment to maintain a constant pulse duration was greatly
% aided by studying the ss_b1verse.m code from Adam Kerr and Peder Larson.
% This code can be found here:
% https://github.com/LarsonLab/Spectral-Spatial-RF-Pulse-Design/blob/master/ss_b1verse.m

function [b1v, gv] = minsarverse_matlab_simpler(b1, g, dt, gmax, smax, dmax)

if ~exist('dmax','var') || isempty(dmax)
    dmax = inf;
end
dmin = 1/dmax;

n  = length(b1); % Number of B1 points
br = real(b1);   % Get real...
bi = imag(b1);   % ...and imaginary b1
b1abs = abs(b1); % Max abs b1

% Do some error checking
if n<=0; error('n must be greater than zero.'); end
if dt<0; error('dt must be positive.');         end
if gmax<=0; error('gmax must be greater than zero.'); end
if smax<=0; error('smax must be greater than zero.'); end
if any(abs(g)>gmax); error('No element of g should be > gmax!'); end

% Initialise
indfix  = false(n,1);
maxiter = 100;

% Original timing information
dt = dt * ones(n,1);    % Time-widths
t  = cumsum(dt) - dt/2; % Cumulative time array
T  = sum(dt);           % Total pulse duration

% Make gradient positive and revert sign change at the end
gsign = (g./abs(g)); % Save sign of gradient point
gsign(g==0) = 1;     % Fix division by zero
g     = abs(g);      % Positive only gradient

% Generate minimum SAR B1 subject to gmax constraint
b1flat      = sum(b1abs) / sum(b1abs>0); % Calculate non-zero mean b1
dilate      = b1abs / b1flat;            % Calculate time dilation function
ind         = b1abs == 0;                % At points of zero RF...
dilate(ind) = 1;                         % ...don't try to adjust
ind         = (g./dilate) > gmax;         % If gmax will be exceeded...
dilate(ind) = g(ind) / gmax;             % ...simply maximize gradient

% Deal with cases where time is dilated by more dmax
if dmax < inf
    niter = 0;
    dilateold = dilate;
    while (niter < maxiter)

        % Reset fix checker
        indfix(:) = false;

        % Decrease dilation in areas exceeding the max dilation
        indmax = dilate>=dmax;
        indmin = dilate<=dmin;
        dilate(indmax) = dmax;
        dilate(indmin) = dmin;
        indfix(indmax|indmin) = true;

        % Shrink non-adjusted points to maintain a fixed duration
        dilate = adjusttimedilate(dilate, dt.*dilate, T, indfix);

        % Reset b1flat of non-fixed points
        b1flat = mean(b1abs(~indfix)./dilate(~indfix));

        % Recalculate dilation of indmax regions with new b1flat to avoid b1flat being larger than
        % max dilated regions
        ind = b1flat > b1abs(indmax)./dilate(indmax);
        idx = find(indmax);
        dilate(idx(ind)) = b1abs(idx(ind)) / b1flat;

        % Recalculate dilation of indmin regions with new b1flat to avoid b1flat being smaller than
        % min dilated regions
        ind = b1flat < b1abs(indmin)./dilate(indmin);
        idx = find(indmin);
        dilate(idx(ind)) = b1abs(idx(ind)) / b1flat;

        % Check if there were any changes, then stop
        if sum(dilate-dilateold>1e-12) == 0
            break; % Finished
        end
        %fprintf('niter = %s, nfix = %s.\n',num2str(niter),num2str(sum(indfix)));

        dilateold = dilate;
        niter = niter + 1;
    end
end

brv = br ./ dilate;
biv = bi ./ dilate;
gv  = g  ./ dilate;
dtv = dt .* dilate;

% Deal with smax violations and further gmax violations
niter = 0;
while (niter < maxiter)

    % Reset fix checker
    indfix(:) = false;

    % gmax checks
    dilate = ones(n,1);
    ind     = gv > gmax;
    dilate(ind) = gv(ind) / gmax;
    brv     = brv ./ dilate;
    biv     = biv ./ dilate;
    gv      = gv  ./ dilate;
    dtv     = dtv .* dilate;
    indfix(ind) = true;

    % smax checks (forward)
    for ii = 2:n
        SR = (gv(ii)-gv(ii-1)) / (dtv(ii)+dtv(ii-1)) / 0.5; % slew rate
        if (SR > smax)
            dilate    = stretchslew(gv(ii), gv(ii-1), dtv(ii), dtv(ii-1), smax);
            brv(ii)    = brv(ii) / dilate;
            biv(ii)    = biv(ii) / dilate;
            gv(ii)     = gv(ii)  / dilate;
            dtv(ii)    = dtv(ii) * dilate;
            indfix(ii) = true;
        end
    end

    % smax checks (backwards)
    for ii = n-1:-1:1
        SR = (gv(ii)-gv(ii+1)) / (dtv(ii)+dtv(ii+1)) / 0.5; % slew rate
        if (SR > smax)
            dilate    = stretchslew(gv(ii), gv(ii+1), dtv(ii), dtv(ii+1), smax);
            brv(ii)    = brv(ii) / dilate;
            biv(ii)    = biv(ii) / dilate;
            gv(ii)     = gv(ii)  / dilate;
            dtv(ii)    = dtv(ii) * dilate;
            indfix(ii) = true;
        end
    end

    % Check if there were any gmax or slew rate violations
    if sum(indfix) == 0
        break; % Finished
    end
    %fprintf('niter = %s, nfix = %s.\n',num2str(niter),num2str(sum(indfix)));

    % Shrink non-adjusted points to maintain a fixed duration
    [brv, biv, gv, dtv] = adjusttime(brv, biv, gv, dtv,  gmax, T, indfix);

    niter = niter + 1;
end

% Revert gradient sign change
gv = gv .* gsign;

% Resample B1v and Gv to original time-widths
tv  = cumsum(dtv) - dtv/2; % Cumulative time array
brv = interp1(tv, brv, t, 'linear', 'extrap');
biv = interp1(tv, biv, t, 'linear', 'extrap');
gv  = interp1(tv, gv , t, 'linear', 'extrap');

% Re-combine complex waveform
b1v = brv + 1i*biv;

    % Shrink non-adjusted points to maintain a fixed duration
    function dilate = adjusttimedilate(dilate, dt, T, indfix)
        indfix([1,end]) = true;             % Don't adjust first and last values!

        Tp     = sum(dt);                   % Total pulse duration
        Tpfix  = sum(dt(indfix));           % Total fixed pulse duration
        factor = (T-Tpfix) / (Tp-Tpfix);    % Fraction to reduce blocks by

        dilate(~indfix) = dilate(~indfix) * factor; % Only adjust blocks that aren't fixed
    end

    % Compress non-adjusted points to maintain a fixed duration
    function [br, bi, g, dt] = adjusttime(br, bi, g, dt, gmax, T, indfix)
        indfix([1,end]) = true;             % Don't adjust first and last values!
        indfix(sqrt(br(count)^2+bi(count)^2)<1e-12) = true; % Don't adjust points of zero RF!
        indfix(abs(g-gmax)<1e-12) = true;   % Don't adjust points of max gradient!

        Tp     = sum(dt);                   % Total pulse duration
        Tpfix  = sum(dt(indfix));           % Total fixed pulse duration
        factor = (T-Tpfix) / (Tp-Tpfix);    % Fraction to reduce blocks by

        br(~indfix) = br(~indfix) / factor; % Only adjust blocks that aren't fixed
        bi(~indfix) = bi(~indfix) / factor;
         g(~indfix) = g(~indfix)  / factor;
        dt(~indfix) = dt(~indfix) * factor;
    end

    % Adjust gradient amplitude to meet slew rate constraint
    function dilate = stretchslew(gh, gl, dth, dtl, smax)
        % Shrink gh to meet smax but preserve area.
        % Need gh/k = gl+maxSR*(dth*k+dtl)/2
        % or (maxSR*dth)*k^2  + (smax*dtl+2gl)*k  - 2gh = 0,
        % where k is time dilation factor to solve for.
        % This is equivalent to a*k^2 + b*k + c = 0
        % Roots: r = b/2a ± sqrt(b*b-4*a*c)/2a.
        a = smax * dth;
        b = smax * dtl + 2*gl;
        c = -2*gh;
        dilate = (-b + sqrt(b*b-4*a*c))/(2*a);

        % Add small adjustment to slew to deal with floating point issues
        % and ensure the slew rate is not exceeded.
        dilate = dilate + 1e-12;
    end

end
