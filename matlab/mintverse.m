% Convert an RF/gradient pair to the minimum-time VERSE equivalent pair.
%
% [b1v, gv] = mintverse(b1, g, dt, bmax, gmax, smax, emax)
%
% 1) Compresses B1/gradient so one is always maximized.
% 2) Removes gradient slew rate violations.
% 3) Resamples waveforms uniformly.
%
% in:
%     b1   - n-point complex RF waveform (arbitrary units)
%     g    - n-point gradient waveform (arbitrary units)
%     dt   - time increment (arbitrary units)
%     n    - number of points in b1, g and dt
%     bmax - max RF value (units of b1)
%     gmax - max gradient value (units of g)
%     smax - max gradient slew rate (units of g per unit dt)
%     emax - max RF energy (units of b1*b1*dt) (-1 to not constrain)
% 
% out:
%     b1v  - VERSEd complex B1 pulse
%     gv   - VERSEd gradient
% 
% Useful information:
%   1. b1 and g must be the same length with identical dt.
%   2. Points of zero RF or g will be left unchanged
%   2. The gradient can be positive, negative or zero.
%   3. b1 and g amplitudes are defined at the centre of each dt.
%   4. If the RF has steep phase ramps (e.g. far off-isocentre pulse or
%      adiabatic), better results can be achieved if the waveform has
%      small dt. The waveforms can then be downsampled to the required
%      dt after VERSEing.
% 
% Adapted from Brian Hargreaves' mintverse.c algorithm 
% by Joseph G. Woods, June 2022
%
% JGW changes:
%   1. Added support for negative gradients.
%   2. Simplified resample() function. Now avoids infinite loop issues.
%   3. Added adjusttime() to correct final pulse duration not matching raster.
%   4. Removed unused support functions.
% 
% Paper: Hargreaves et al. MRM 2004, https://doi.org/10.1002/mrm.20168
% Original code: http://mrsrl.stanford.edu/~brian/mintverse/ 

function [b1v, gv] = mintverse(b1, g, dt, bmax, gmax, smax, emax)

DEBUG = false;
% DEBUG = true;

n  = length(b1); % Number of B1 points
br = real(b1);   % Get real...
bi = imag(b1);   % ...and imaginary b1 to VERSE.

if (n<=0); fprintf("Warning: mintverse(): n must be greater than zero. \n"); end
if (dt<=0); fprintf("Warning: mintverse(): dt must be positive.\n"); end
if (~(isnonnegative(g,n))); fprintf("Warning:  mintverse(): g values must be non-negative.\n"); end
if (bmax<=0); fprintf("Warning: mintverse(): bmax must be greater than zero. \n"); end
if (gmax<=0); fprintf("Warning: mintverse(): gmax must be greater than zero. \n"); end
if (smax<=0); fprintf("Warning: mintverse(): smax must be greater than zero. \n"); end

if ~exist('emax','var') || isempty(emax)
    emax = -1; % Maximum pulse energy
    fprintf("No max B1 energy constraint given.\n");
end

% Set fix counter
indfix = zeros(n,1);

emaxratio = 0.98; % minumum energy/emax for stopping iteration
niter     = 0;    % number of iterations
maxiter   = 1000; % max number of iterations

bmaxc = bmax;
bmaxh = bmax; % maximum bmax for energy OR b-amplitude
bmaxl = 0;

% Make gradient positive and revert at the end.
[g, gsign] = gradientsign(g, n); % Save sign of the gradient points

while (niter < maxiter)

    % copy inputs to the arrays
    brwork = br;
    biwork = bi;
    gwork  = g;
    dtwork = arraycopy(dt, n);

    [brwork,biwork,gwork,dtwork] = compressmax(brwork,biwork,gwork,dtwork,n,bmaxc,gmax);

    [brwork,biwork,gwork,dtwork] = slewcheck(brwork,biwork,gwork,dtwork,n,smax);

    if (emax > 0)

        benergy = calcenergy(brwork, biwork, dtwork, n);

        if (benergy > emax) % reduce bmax to reduce energy
            bmaxh = bmaxc;
        else % increase bmax to increase energy
            if (benergy/emax > emaxratio) % close enough
                break; % finished
            end
            bmaxl = bmaxc;
        end
        bmaxc = bmaxl + (bmaxh - bmaxl)/2.0; % update bmax

        if (bmaxl >= bmaxh)
            fprintf("Warning: mintverse() low-E limit > high-E limit.\n");
            fprintf("Exiting after %d iterations.\n",niter);
            break;
        end
        if (niter >= maxiter)
            fprintf("Warning: mintverse() iteration limit reacheed.\n");
            fprintf("Exiting after %d iterations.\n",niter);
            break;
        end

        if DEBUG
            fprintf("mintverse(): niter = %d, bmax: high/low/next = %g/%g/%g\n",niter,bmaxh,bmaxl,bmaxc);
        end

    else
        break; % no energy constraint, so done!
    end

    niter = niter + 1;
end

% Revert gradient sign change
gwork = multiplyarray(gwork,gsign,n);

% Calculate length of, and allocate output vectors.
nout  = uniformresamplesize(dtwork,n,dt);
dtout = arraycopy(dt, nout);

% Adjust time-points to deal with sum(dtwork) < dt*nout
T = dt * nout;
[brwork, biwork, gwork, dtwork] = adjusttime(brwork, biwork, gwork, dtwork, n, T, indfix);

% Resample with uniform sample step.
twork = cumsumtime(dtwork,n   ); % Calculate VERSEd cumulative time */
tout  = cumsumtime(dtout ,nout); % Calculate output cumulative time */
brv   = resample(brwork,twork,n,tout,nout);
biv   = resample(biwork,twork,n,tout,nout);
gv    = resample(gwork ,twork,n,tout,nout);

% Re-combine complex waveform
b1v = brv + 1i*biv;


function asum = sumarray(a, n)
% Function returns the sum of the first n points in the array a.
    asum = 0;
    for count = 1:n
        asum = asum + a(count);
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


function g = multiplyarray(g, gsign, n)
% Function multiplies the elements of the array g by the elements of
% the array gsign, updating g with the result.
    for count = 1:n
        g(count) = g(count) * gsign(count);
    end
end


function t = cumsumtime(dt, n)
% Function returns the cumulative time at the centres of blocks with
% widths dt.
    t = zeros(n,1);

    % Count time to halfway through blocks
    t(1) = dt(1) * 0.5;

    for count = 2:n
        t(count) = t(count-1) + (dt(count-1)+dt(count))*0.5;
    end
end


function b = isnonnegative(f, n)
% Function checks the first n points of the array f.
% If they are all non-negative, it returns 1, otherwise returns 0.

    for count = 1:n
        if (f(count) < 0)
            b = 0;
            return;
        end
    end
    b = 1;
end


function dtarray = arraycopy(dt, n)
% Function copies dt to the first n points of the array dtwork.

    dtarray = zeros(n,1);
    for count = 1:n
        dtarray(count) = dt;
    end
end


function [g, gsign] = gradientsign(g, n)
% Function calculates the sign (1 or -1) of the gradient and
% updates the gradient to be positive only.

    gsign = zeros(n,1);
    for count = 1:n
        if (g(count) >= 0); gsign(count) =  1;
        else;               gsign(count) = -1;
        end
    end
end


function nfr = uniformresamplesize(dt, n, dtr)
% Function returns the number of points that will are required
% to uniformly resample an array with block-widths in *dt,
% with a uniform time block width dtr.
% 
% in:
%     dt  - n-point array of block widths for input
%     n   - number of points in input
%     dtr - block width for each output block
% 
% out:
%     nfr - number of output points/blocks
        
    tend = sumarray(dt, n);
    nfr = floor((tend+0.99999*dtr) / dtr); % Round up end time

end


function fr = resample(f, t, n, tr, nr)
% Function resamples f onto a different time step.
% The input, f, is linearly interpolated between the centers of
% blocks.
% 
% in:
%     f  - n-point array of heights of block-centers to resample
%     t  - n-point array of block centres for input
%     n  - number of points in input and output
%     fr - nr-point array of resampled heights of block-centers
%     tr - nr-point array of block centres for each output block
%     nr - number of points in output (assumed <=n)
% 
% out:
%     fr - array of heights of output blocks.

    % Manual linear interpolation with linear extrapolation
    fr = zeros(nr,1);
    count = 2;
    for countr = 1:nr

            % Find smallest count, where t(count) > tr(countr)
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


% Shrink non-adjusted points to maintain a fixed duration
function [br, bi, g, dt] = adjusttime(br, bi, g, dt, n, T, indfix)
% Function stretches dt to full time because nout is rounded up.
% 
% in:
%     br - real part of RF waveform (arbitrary units)
%     bi - imag part of RF waveform (arbitrary units)
%     g  - gradient waveform (arbitrary units)
%     dt - time-delta waveform (arbitrary units)
%     nr - number of points in arrays
%     T  - total rounded up time (nout*dtout)
% 
% out:
%     br,bi,g,dt are all updated.

    indfix(1) = 1;                      % don't adjust first value!
    indfix(n) = 1;                      % don't adjust last value!
    
    Tp     = sumarray(dt,n);            % total pulse duration
    Tpfix  = sumarraybool(dt,n,indfix); % total fixed pulse duration
    shrink = (T-Tpfix) / (Tp-Tpfix);    % fraction to adjust blocks by
    if DEBUG
        fprintf("adjusttime(): shrink = %f, T = %f, Tp = %f, Tpfix = %f.\n",shrink,T,Tp,Tpfix);
    end
    
    % Only adjust blocks that aren't fixed
    for count = 1:n
        if (indfix(count) == 0)
            br(count) = br(count) / shrink;
            bi(count) = bi(count) / shrink;
            g(count)  = g(count)  / shrink;
            dt(count) = dt(count) * shrink;
        end
    end
end


function stretch = stretchslew(gh, gl, dth, dtl, smax)
% Function finds the appropriate stretch in time of a gradient block of
% height gh and width dth that is adjacent to a gradient block of height
% gl and width dtl, so that the slew rate smax exactly met, while the area
% of the block is preserved.
% 
% in:
%     gh   - higher gradient (assumed non-negative).
%     gl   - lower gradient (assumed non-negative).
%     dth  - width of higher gradient segment.
%     dtl  - width of lower gradient segment.
%     smax - maximum slew rate.
% 
% Units: Arbitrary, but smax must be in gradient units per time unit.
% 
% out:
%     stretch - factor to stretch the higher gradient block width in TIME!
% 
% notes:
%     Shrink gh to meet smax but preserve area.
%     Need gh/k = gl+maxSR*(dth*k+dtl)/2
%     or (maxSR*dth)*k^2  + (smax*dtl+2gl)*k  - 2gh = 0,
%     where k is time stretch factor to solve for.
%     This is equivalent to a*k^2 + b*k + c = 0
%     Roots: r = b/2a ± sqrt(b*b-4*a*c)/2a.

    % Define a,b,c for quadratic equation a*k^2+b*k+c = 0
    a = smax*dth;
    b = smax*dtl + 2*gl;
    c = -2*gh;
    stretch = (-b + sqrt(b*b-4*a*c))/(2*a);

    % Add small adjustment to slew to deal with floating point issues
    %  and ensure the slew rate is not exceeded.
    stretch = stretch + 1e-12;

end


function [br,bi,g,dt] = adjustslew(br, bi, g, dt, n, smax, loc, recursecount)
% Fix slew rate violations, if they exist, by stretching
%     the block width so that the slew rate is exactly met.
% 
% in:
%     br     - real part of RF waveform (arbitrary units)
%     bi     - imag part of RF waveform (arbitrary units)
%     g      - gradient waveform (arbitrary units)
%     dt     - time-delta waveform (arbitrary units)
%     n      - number of points in b, g and dt
%     smax   - max slew rate in gradient units per time unit
%     loc    - current point to adjust (1 <= loc < n)
% 
% out:
%     br,bi,g,dt are all adjusted.
% 
% note: recursive!

    % Do some error checking
    if (loc <= 1)
        fprintf("Warning: adjustslew() - loc must be positive. \n");
    end

    if (g(loc-1) > g(loc)) % Downward slope Reduce height of higher
                           % point and work backwards.
    
        stretch = stretchslew(g(loc-1),g(loc),dt(loc-1),dt(loc),smax);
        g(loc-1)  = g(loc-1)  / stretch;
        br(loc-1) = br(loc-1) / stretch;
        bi(loc-1) = bi(loc-1) / stretch;
        dt(loc-1) = dt(loc-1) * stretch;

        if (loc > 2)
            slew = (g(loc-1)-g(loc-2)) / (dt(loc-1)+dt(loc-2)) / 0.5;
            if (abs(slew) > smax)
                [br,bi,g,dt] = adjustslew(br,bi,g,dt,n,smax,loc-1,recursecount+1);
            end
        end

    else % Upward slope. Just reduce height of higher point. */
        stretch = stretchslew(g(loc),g(loc-1),dt(loc),dt(loc-1),smax);
        g(loc)  = g(loc)  / stretch;
        br(loc) = br(loc) / stretch;
        bi(loc) = bi(loc) / stretch;
        dt(loc) = dt(loc) * stretch;
    end
end


function [br,bi,g,dt] = slewcheck(br, bi, g, dt, n, smax)
% Adjust pulse so that slew rates are not violated
% 
% in:
%     br     - real part of RF waveform (arbitrary units.)
%     bi     - imag part RF waveform (arbitrary units.)
%     g      - gradient waveform (arbitrary units.)
%     dt     - time-delta waveform (arbitrary units.)
%     n      - number of points in b, g and dt.
%     smax   - max slew rate in gradient units per time unit.
% 
% out:
%     br,bi,g,dt are all adjusted.

    for count = 2:n
        slew = (g(count)-g(count-1)) / (dt(count)+dt(count-1)) / 0.5;
        if (abs(slew) > smax)
            [br,bi,g,dt] = adjustslew(br,bi,g,dt,n,smax,count,0);
        end
    end
end


function benergy = calcenergy(br, bi, dt, n)
% Calculate the energy in a waveform b, with time steps
% in dt and n points.  This is just sum(b_i^2*dt_i).
% 
% in:
%     br - real part of waveform, nominally RF (arb. units)
%     bi - imag part of waveform, nominally RF (arb. units)
%     dt - time-widths of sample blocks in b (arb units)
%     n  - number of points in b and dt
% 
% out:
%     benergy - energy in waveform

    benergy = 0.0;

    for count = 1:n
        benergy = benergy + (br(count)*br(count)+bi(count)*bi(count)) * dt(count);
    end

end


function [br,bi,g,dt] = compressmax(br, bi, g, dt, n, bmax, gmax)
% Stretch or shrink the n block widths (dt) such that at each
% time step either |b| or g reaches its maximum.  Does not
% change the length of |b|, g or dt.
% br, bi and bmax must be the same units, while g and gmax must
%     also be the same units.
% 
% in:
%     br   - real part of RF waveform (arbitrary units)
%     bi   - imag part of RF waveform (arbitrary units)
%     g    - gradient waveform (arbitrary units.)
%     dt   - time-delta waveform (arbitrary units)
%     n    - number of points in b, g and dt
%     bmax - max RF value (units of b)
%     gmax - max gradient value (units of g)
% 
% out:
%     br,bi,g,dt are all updated.

    for count = 1:n

        bfrac = abs( sqrt(br(count)*br(count)+bi(count)*bi(count)) / bmax );
        gfrac = abs( g(count) / gmax );

        if ((bfrac > 0) && (gfrac > 0))
            if (bfrac > gfrac) % RF will be max'd
                br(count) = br(count) / bfrac;
                bi(count) = bi(count) / bfrac;
                g(count)  = g(count)  / bfrac;
                dt(count) = dt(count) * bfrac;

            else                % Gradient will be max'd
                br(count) = br(count) / gfrac;
                bi(count) = bi(count) / gfrac;
                g(count)  = g(count)  / gfrac;
                dt(count) = dt(count) * gfrac;
            end
        else
            fprintf("compressmax(): Zero-RF or zero-gradient (%d), ignored.\n",count);
        end
    end
end

end