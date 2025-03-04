%% VERSE RF to a minimum SAR form and adjust gradient to meet constraints.
%
% [b1v, gv] = minsarverse_matlab(b1, g, dt, gmax, smax)
%
% in:
%      b1   - complex B1+ waveform (units, array)
%      g    - gradient waveform (units/distance, array)
%      dt   - b1 and gradient waveform raster or resolution (s, scalar)
%      gmax - maximum gradient amplitude permitted (units/distance, scalar)
%      smax - maximum gradient slew rate permitted (units/distance/s, scalar)
%
% out:
%      b1v  - minimum-SAR VERSEd complex B1+ waveform (units, array)
%      gv   - minimum-SAR VERSEd gradient waveform (units/distance, array)
%
% notes:
%   This is a matlab version of minsarverse.c and avoids in-built matlab
%   functionality to follow the C-code code closely.
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
%   1. Add option to set maximum and minimum scaling of the RF to reduce
%      extreme scalings.
%   2. Add option to low-pass filter the gradient waveform and then adjust
%      the RF correspondingly to avoid hardware issues (Hargreaves used 50 kHz).
%
% Written by Joseph G. Woods, University of Oxford, May 2022
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
% https://github.com/LarsonLab/Spectral-Spatial-RF-Pulse-Design

function [b1v, gv] = minsarverse_matlab(b1, g, dt, gmax, smax)

n  = length(b1); % Number of B1 points
br = real(b1);   % Get real...
bi = imag(b1);   % ...and imaginary b1 to VERSE.

% Do some error checking
if n<=0; error('n must be greater than zero.');    end
if dt<0; error('dt must be positive.');            end
if gmax<=0; error('gmax must be greater than zero.'); end
if smax<=0; error('smax must be greater than zero.'); end
if exceedsgmax(g,gmax,n); error('No element of g should be > gmax!'); end

% Set fix counter
indfix = zeros(n,1);

% Original timing information
dt = dt * ones(n,1);   % Time-widths
t  = cumsumtime(dt,n); % Cumulative time array
T  = sumarray(dt,n);   % Total pulse duration

% Make gradient positive and revert at the end.
[g, gsign] = gradientsign(g, n); % Save sign of the gradient points

% Generate minimum SAR B1 subject to gmax constraint
[brwork, biwork, gwork, dtwork] = minimizesar(br, bi, g, dt, n, gmax);

% Shrink non-adjusted points to maintain a fixed duration
[brwork, biwork, gwork, dtwork] = adjusttime(brwork, biwork, gwork, dtwork, gmax, n, T, indfix);

niter = 0;
maxiter = 100;
while (niter < maxiter)

    % Reset fix checker
    for ii = 1:n
        indfix(ii) = 0;
    end

    % Fix gmax violations
    [brwork, biwork, gwork, dtwork, indfix] = gmaxcheck(brwork, biwork, gwork, dtwork, n, gmax, indfix);

    % Fix slew rate violations
    [brwork, biwork, gwork, dtwork, indfix] = slewcheck(brwork, biwork, gwork, dtwork, n, smax, indfix);

    % Check if there were any gmax or slew rate violations
    if sumarray(indfix,n) == 0 
        break; % Finished
    end
    fprintf('minsarverse(): niter = %s, nfix = %s.\n',num2str(niter),num2str(sumarray(indfix,n)));

    % Shrink non-adjusted points to maintain a fixed duration
    [brwork, biwork, gwork, dtwork] = adjusttime(brwork, biwork, gwork, dtwork, gmax, n, T, indfix);

    niter = niter + 1;
end

% Revert gradient sign change
gwork = multiplyarray(gwork,gsign,n);

% Resample B1v and Gv to original time-widths
tv  = cumsumtime(dtwork,n);
brv = resample(brwork, tv, n, t);
biv = resample(biwork, tv, n, t);
gv  = resample(gwork , tv, n, t);

% Re-combine complex waveform
b1v = brv + 1i*biv;

    % Generate minimum SAR B1 subject to gmax constraint
    function [br, bi, g, dt] = minimizesar(br, bi, g, dt, n, gmax)

        bint  = sumabscomplexarray(br,bi,n); % Calculate magnitude integral of abs(RF)
        bn    = countnonzeroarrays(br,bi,n); % Count non-zero RF points
        bflat = bint / bn;                   % RF magnitude for minimum SAR

        for count = 1:n
            absb1 = sqrt( br(count)*br(count)+bi(count)*bi(count) );
            if ( absb1 > 0 ) % Don't adjust zero RF points

                frac = absb1 / bflat; % Adjust waveform to achieve bflat

                % If gmax will be exceeded, simply maximize gradient
                if ( g(count)/frac) > gmax
                    frac = g(count) / gmax;
                end

                br(count) = br(count) / frac;
                bi(count) = bi(count) / frac;
                g(count)  =  g(count) / frac;
                dt(count) = dt(count) * frac;
            else
                print("minimizesar(): Zero-RF (%d) ignored.\n",count);
            end
        end
    end

    function [br, bi, g, dt, indfix] = gmaxcheck(br, bi, g, dt, n, gmax, indfix)

        for count = 1:n
            if g(count) > gmax % If gmax is exceeded...
                stretch   = g(count) / gmax; % ...simply maximize gradient
                br(count) = br(count) / stretch;
                bi(count) = bi(count) / stretch;
                g(count)  = g(count)  / stretch;
                dt(count) = dt(count) * stretch;
                indfix(count) = 1;
            end
        end
    end

    function [br, bi, g, dt, indfix] = slewcheck(br, bi, g, dt, n, smax, indfix)

        for count = 2:n
            slew = (g(count)-g(count-1)) / (dt(count)+dt(count-1)) / 0.5;
            if abs(slew) > smax
                [br, bi, g, dt, indfix] = adjustslew(br, bi, g, dt, smax, count, 0, indfix);
            end
        end
        
    end

    function [br, bi, g, dt, indfix] = adjustslew(br, bi, g, dt, smax, loc, recursecount, indfix)

        if (loc <= 1)
            error('adjustslew() - loc must be positive.');
        end

        if (g(loc-1) > g(loc)) % Downward slope. Reduce height of higher 
                               % point and work backwards.

            stretch = stretchslew(g(loc-1),g(loc),dt(loc-1),dt(loc),smax);
            br(loc-1) = br(loc-1) / stretch;
            bi(loc-1) = bi(loc-1) / stretch;
            g(loc-1)  = g(loc-1)  / stretch;
            dt(loc-1) = dt(loc-1) * stretch;
            indfix(loc-1) = 1;

            if loc > 2
                slew = (g(loc-1)-g(loc-2)) / (dt(loc-1)+dt(loc-2)) / 0.5;
                if (abs(slew) > smax)
                    [br, bi, g, dt, indfix] = adjustslew(br,bi,g,dt,smax,loc-1,recursecount+1,indfix);
                end
            end

        else % Upward slope. Simply reduce height of higher point.

            stretch = stretchslew(g(loc),g(loc-1),dt(loc),dt(loc-1),smax);
            br(loc) = br(loc) / stretch;
            bi(loc) = bi(loc) / stretch;
            g(loc)  = g(loc)  / stretch;
            dt(loc) = dt(loc) * stretch;
            indfix(loc) = 1;
        end
        
    end

    % Adjust gradient amplitude to meet slew rate constraint
    function stretch = stretchslew(gh, gl, dth, dtl, smax)
        % Shrink gh to meet smax but preserve area.
        % Need gh/k = gl+maxSR*(dth*k+dtl)/2
        % or (maxSR*dth)*k^2  + (smax*dtl+2gl)*k  - 2gh = 0,
        % where k is time stretch factor to solve for.
        % This is equivalent to a*k^2 + b*k + c = 0
        % Roots: r = b/2a ± sqrt(b*b-4*a*c)/2a.
        a = smax*dth;
        b = smax*dtl + 2*gl;
        c = -2*gh;
        stretch = (-b + sqrt(b*b-4*a*c))/(2*a);

        % Add small adjustment to slew to deal with floating point issues
        % and ensure the slew rate is not exceeded.
        stretch = stretch + 1e-12;
    end

    % Shrink non-adjusted points to maintain a fixed duration
    function [br, bi, g, dt] = adjusttime(br, bi, g, dt, gmax, n, T, indfix)

        indfix(1) = 1;                      % don't adjust first value!
        indfix(n) = 1;                      % don't adjust last value!
        for count = 1:n
            if (abs(g(count)-gmax) < 1e-12) 
                indfix(count) = 1;          % don't adjust points of max gradient!
            end
        end

        Tp     = sumarray(dt,n);            % total pulse duration
        Tpfix  = sumarraybool(dt,n,indfix); % total fixed pulse duration
        shrink = (T-Tpfix) / (Tp-Tpfix);    % fraction to adjust blocks by

        % Only adjust blocks that weren't fixed for gmax or slew rate issues
        for count = 1:n
            if (indfix(count) == 0)
                br(count) = br(count) / shrink;
                bi(count) = bi(count) / shrink;
                g(count)  = g(count)  / shrink;
                dt(count) = dt(count) * shrink;
            end
        end
    end

    % Manual linear interpolation with linear extrapolation
    function fr = resample(f, t, n, tr)

        fr = zeros(size(tr));
        count = 2;
        for countr = 1:n

            % Find smallest jj, where t(jj) > tr(ii)
            while (tr(countr) > t(count))
                if (count == n)
                    break;
                end
                count = count + 1;
            end

            % fraction of the countr interval to interpolate into
            xr = (tr(countr)-t(count-1)) / (t(count)-t(count-1));

            % interpolated value
            fr(countr) = (1-xr)*f(count-1) + xr*f(count);
        end

    end

    % Function calculates the sign (1 or -1) of the gradient and updates
    % the gradient to be positive only.
    function [g, gsign] = gradientsign(g, n)

        gsign = zeros(n,1);
        for count = 1:n
            if (g(count) == 0)
                gsign(count) = 1;
            else
                gsign(count) = g(count) / abs(g(count));
                g(count)     = abs(g(count));
            end
        end
    end

    % Function checks the first n points of the array g. If they are all <
    % gmax, it returns 0, otherwise returns 1.
    function val = exceedsgmax(g, gmax, n)
        
        for count = 1:n
            if (abs(g(count)) > gmax)
                val = 1; % An element of g exceeds gmax
                return;
            end
        end
        val =  0;
    end

    % Function returns the cumulative time at the centres of blocks with
    % widths dt.
    function t = cumsumtime(dt, n)
        t = zeros(n,1);
        t(1) = dt(1) * 0.5;
        for count = 2:n
            t(count) = t(count-1) + (dt(count-1)+dt(count))*0.5;
        end
    end

    % Function returns the number of non-zero points in the arrays br and
    % bi.
    function nznum = countnonzeroarrays(br, bi, n)
        nznum = 0;
        for count = 1:n
            if ( br(count)~=0 || bi(count)~=0 )
                nznum = nznum + 1;
            end
        end
    end

    % Function returns the complex magnitude sum of the first n points in
    % the arrays br and bi.
    function bsum = sumabscomplexarray(br, bi, n)
        bsum = 0;
        for count = 1:n
            bsum = bsum + sqrt( br(count)*br(count) + bi(count)*bi(count) );
        end
    end

    % Function returns the conditional sum of the first n points in the
    % array dt. Conditional on indfix[count]==1.
    function dtsum = sumarraybool(dt, n, indfix)
        dtsum = 0;
        for count = 1:n
            if (indfix(count) == 1)
                dtsum = dtsum + dt(count);
            end
        end
    end

    % Function returns the sum of the first n points in the array a.
    function asum = sumarray(a, n)
        asum = 0;
        for count = 1:n
            asum = asum + a(count);
        end
    end

    % Function multiplies the elements of the array g by the elements of
    % the array gsign, updating g with the result.
    function g = multiplyarray(g, gsign, n)
        for count = 1:n
            g(count) = g(count) * gsign(count);
        end
    end

end
