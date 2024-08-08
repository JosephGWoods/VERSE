/* ----------------------------------------------------------------------
    minsarverse.c -- Minimum-SAR VERSE functions.

    Useful information:
        1. br, bi, and g must be the same length with identical dt.
        2. The gradient can be positive, negative or zero.
        3. br, bi, and g amplitudes are defined at the centre of each dt.
        4. Points of zero b1 and g magnitude are mostly left alone, since we
           assume they are there for a purpose (e.g. gradient ramps of a
           certain length for coil lead times). However, non-zero b1 is
           currently adjusted by minimizesar() at points of zero gradient.
        5. If the RF has steep phase ramps (e.g. far off-isocentre pulse or
           adiabatic), better results can be achieved if the waveform has
           small dt. The waveforms can then be resampled to the required
           dt after VERSEing.

    TODO:
        1. Add option to set maximum and minimum scaling of the RF to reduce
           extreme scalings.
        2. Add option to low-pass filter the gradient waveform and then adjust
           the RF correspondingly to avoid hardware issues (Hargreaves used 50 kHz).

    Written by Joseph G. Woods, University of Oxford, May 2022

    Following a similar procedure to algorithm described in Conolly et al.
    JMR 1988 (https://doi.org/10.1016/0022-2364(88)90131-X).

    Adapted from Brian Hargreaves' mintverse algorithm, described in
    Hargreaves et al. MRM 2004 (https://doi.org/10.1002/mrm.20168). The code
    can be found here: http://mrsrl.stanford.edu/~brian/mintverse/

    My understanding of the VERSE algorithm in general and implementation of
    the time adjustment to maintain a constant pulse duration was greatly
    aided by studying the ss_b1verse.m code from Adam Kerr and Peder Larson.
    This code can be found here:
    https://github.com/LarsonLab/Spectral-Spatial-RF-Pulse-Design
    
    Please use this software freely with appropriate credit or acknowledgement.
    ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
    Original statement from Brian Hargreaves' mintverse code:

    mintverse.c -- Minimum-time VERSE functions.

    Copyright (c)  Brian A. Hargreaves and Charles H. Cunningham

	The main function of interest is mintverse().  Other support
	functions are included here.

	Please use this software freely.  
	Appropriate credit or acknowledgement is certainly appreciated.
   ---------------------------------------------------------------------- */

/*
#define DEBUG_VERSE
*/

#include <math.h>   /* for M_PI   */
#include <string.h> /* for memcpy */

void multiplyarray(double *g, long *gsign, long n)
/* ------------------------------------------------------------
   Function multiplies the elements of the array g by the elements
   of the array gsign, updating g with the result.
   ------------------------------------------------------------ */
{
    for (long count = 0; count < n; count++)
        g[count] *= gsign[count];
}


double sumarrayd(double *dt, long n)
/* ------------------------------------------------------------
   Function returns the sum of the first n points in the array dt.
   ------------------------------------------------------------ */
{
    double dtsum = 0.0;

    for (long count = 0; count < n; count++)
        dtsum += dt[count];

    return dtsum;
}


long sumarrayl(long *indfix, long n)
 /* ------------------------------------------------------------
    Function returns the sum of the first n points in the array indfix.
    ------------------------------------------------------------ */
{
    long indfixsum = 0;

    for (long count = 0; count < n; count++)
        indfixsum += indfix[count];

    return indfixsum;
}


double sumarraycond(double *dt, long n, long *indfix)
 /* ------------------------------------------------------------
    Function returns the conditional sum of the first n points in the
    array dt. Conditional on indfix[count]==1.
    ------------------------------------------------------------ */
{
    double dtsum = 0.0;

    for (long count = 0; count < n; count++)
        if (indfix[count] == 1)
            dtsum += dt[count];

    return dtsum;
}


double sumabscomplexrfarrays(double *br, double *bi, long n)
    /* ------------------------------------------------------------
    Function returns the complex magnitude sum of the first n points
    in the arrays br and bi.
    ------------------------------------------------------------ */
{
    double bsum = 0;

    for (long count = 0; count < n; count++)
        bsum += sqrt(br[count]*br[count]+bi[count]*bi[count]);

    return bsum;
}


double countnonzerorfarrays(double *br, double *bi, long n)
    /* ------------------------------------------------------------
    Function returns the number of non-zero points in the arrays br and bi.
    ------------------------------------------------------------ */
{
    double nznum = 0;

    for (long count = 0; count < n; count++)
        if ( (br[count]!=0) || (bi[count]!=0) )
            nznum++;

    return nznum;
}


void cumsumtime(double dt, long n, double *t)
 /* ------------------------------------------------------------
    Function returns the cumulative time at the centres of blocks with
    uniform width dt.
    ------------------------------------------------------------ */
{
    /* Count time to halfway through blocks */
    t[0] = dt * 0.5;

    for (long count = 1; count < n; count++)
        t[count] = t[count-1] + dt;
}


void cumsumtimearray(double *dt, long n, double *t)
    /* ------------------------------------------------------------
    Function returns the cumulative time at the centres of blocks with
    widths dt.
    ------------------------------------------------------------ */
{
    /* Count time to halfway through blocks */
    t[0] = dt[0] * 0.5;

    for (long count = 1; count < n; count++)
        t[count] = t[count-1] + (dt[count-1]+dt[count])*0.5;
}


void arraycopy(double dt, double *dtarray, long n)
/* ------------------------------------------------------------
	Function copies dt to the first n points of the array dtwork.
   ------------------------------------------------------------ */
{
    for (long count = 0; count < n; count++)
        dtarray[count] = dt;
}


long exceedsgmax(double *g, double gmax, long n)
    /* ------------------------------------------------------------
    Function checks the first n points of the array g.
    If they are all < gmax, it returns 0, otherwise returns 1.
    ------------------------------------------------------------ */
{
    for (long count = 0; count < n; count++)
        if (fabs(g[count]) > gmax)
            return 1; /* An element of g exceeds gmax */

    return 0;
}


void gradientsign(double *g, long *gsign, long n)
/* ------------------------------------------------------------
    Function calculates the sign (1 or -1) of the gradient and
    updates the gradient to be positive only.
   ------------------------------------------------------------ */
{
    for (long count = 0; count < n; count++) {
        gsign[count] = (g[count] >= 0) ? 1 : -1; /* Assign 1 or -1 based on the sign of g[count] */
        g[count]     = fabs(g[count]);           /* Update g[count] to be positive */
    }
}


void adjusttimefixgmax(double *br, double *bi, double *g, double *dt,
                       long n, double T, long *indfix, double gmax)
    /* ------------------------------------------------------------
    Function finds the appropriate factor to shrink non-adjusted dt
    points by to maintain a fixed pulse duration

    in:
        br     - real part of RF waveform (arbitrary units).
        bi     - imag part of RF waveform (arbitrary units).
        g      - gradient waveform (arbitrary units).
        dt     - time-delta waveform (arbitrary units).
        n      - number of points in b, g and dt.
        T      - pulse duration to achieve (units of dt).
        indfix - array for fixing durations.
        gmax   - max gradient value (units of g).

    out:
        br,bi,g,dt,indfix are all adjusted.
    ------------------------------------------------------------ */
{
    long count;
    double shrink, Tp, Tpfix;

    /* Never adjust first and last values! */
    indfix[0]   = 1;
    indfix[n-1] = 1;

    for (count = 0; count < n; count++) {
        if (sqrt(br[count]*br[count]+bi[count]*bi[count]) < 1e-12)
            indfix[count] = 1;          /* Don't adjust points of zero RF! */
        if (fabs(g[count]-gmax) < 1e-12)
            indfix[count] = 1;          /* Don't adjust points of max gradient! */
    }

    Tp     = sumarrayd(dt,n);           /* Total pulse duration         */
    Tpfix  = sumarraycond(dt,n,indfix); /* Total fixed pulse duration   */
    shrink = (T-Tpfix) / (Tp-Tpfix);    /* Fraction to reduce blocks by */

    /* Only adjust blocks that weren't fixed for gmax or slew rate issues */
    for (count = 0; count < n; count++) {
        if (indfix[count] == 0) {
            br[count] /= shrink;
            bi[count] /= shrink;
            g[count]  /= shrink;
            dt[count] *= shrink;
        }
    }
}


void adjusttime(double *br, double *bi, double *g, double *dt,
                long n, double T, long *indfix)
    /* ------------------------------------------------------------
    Alternate version of adjusttimefixgmax() below to automatically set 
    gmax to a large value - i.e. ignore the gmax contraint below
    ------------------------------------------------------------ */
{
    adjusttimefixgmax(br,bi,g,dt,n,T,indfix,DBL_MAX);
}


double stretchslew(double gh, double gl, double dth, double dtl, double smax)
    /* ------------------------------------------------------------
    Function finds the appropriate stretch in time of a gradient block of
    height gh and width dth that is adjacent to a gradient block of height
    gl and width dtl, so that the slew rate smax exactly met, while the area
    of the block is preserved.

    in:
        gh   - higher gradient (assumed non-negative).
        gl   - lower gradient (assumed non-negative).
        dth  - width of higher gradient segment.
        dtl  - width of lower gradient segment.
        smax - maximum slew rate.

    Units: Arbitrary, but smax must be in gradient units per time unit.

    out:
        stretch - factor to stretch the higher gradient block width in TIME!

    notes:
        Shrink gh to meet smax but preserve area.
        Need gh/k = gl+maxSR*(dth*k+dtl)/2
        or (maxSR*dth)*k^2  + (smax*dtl+2gl)*k  - 2gh = 0,
        where k is time stretch factor to solve for.
        This is equivalent to a*k^2 + b*k + c = 0
        Roots: r = b/2a ± sqrt(b*b-4*a*c)/2a.
    ---------------------------------------------------------------------- */
{
    double stretch;
    double a,b,c;

    /* Define a,b,c for quadratic equation a*k^2+b*k+c = 0 */
    a = smax*dth;
    b = smax*dtl + 2*gl;
    c = -2*gh;
    stretch = (-b + sqrt(b*b-4*a*c))/(2*a);

    /* Add small adjustment to slew to deal with floating point issues
       and ensure the slew rate is not exceeded. */
    stretch = stretch + 1e-12;

    return stretch;
}


void adjustslew(double *br, double *bi, double *g, double *dt, long n,
               double smax, long loc, long recursecount)
/* ------------------------------------------------------------
    Fix slew rate violations, if they exist, by stretching
    the block width so that the slew rate is exactly met.

    in:
        br     - real part of RF waveform (arbitrary units.)
        bi     - imag part of RF waveform (arbitrary units.)
        g      - gradient waveform (arbitrary units.)
        dt     - time-delta waveform (arbitrary units.)
        n      - number of points in b, g and dt.
        smax   - max slew rate in gradient units per time unit.
        loc    - current point to adjust (1 <= loc < n).

    out:
        br,bi,g,dt are all adjusted.

    note: recursive!
   ------------------------------------------------------------ */
{
    double stretch;
    double slew;

    /* Do some error checking. */
    if (loc < 1 || loc >= n) {
        printf("Error: adjustslew() - loc must be >=1 and <n!\n");
        return;
    }

    if (g[loc-1] > g[loc]) { /* Downward slope Reduce height of higher */
                             /* point and work backwards.              */

        stretch = stretchslew(g[loc-1],g[loc],dt[loc-1],dt[loc],smax);
        g[loc-1]  /= stretch;
        br[loc-1] /= stretch;
        bi[loc-1] /= stretch;
        dt[loc-1] *= stretch;

        if (loc > 1) {
            slew = (g[loc-1]-g[loc-2]) / (dt[loc-1]+dt[loc-2]) / 0.5;
            if (fabs(slew) > smax)
                adjustslew(br,bi,g,dt,n,smax,loc-1,recursecount+1);
        }
    } else { /* Upward slope. Just reduce height of higher point. */

        stretch = stretchslew(g[loc],g[loc-1],dt[loc],dt[loc-1],smax);
        g[loc]  /= stretch;
        br[loc] /= stretch;
        bi[loc] /= stretch;
        dt[loc] *= stretch;
    }
}


void adjustslewindfix(double *br, double *bi, double *g, double *dt, long n,
               double smax, long loc, long recursecount, long *indfix)
    /* ------------------------------------------------------------
    Fix slew rate violations, if they exist, by stretching
    the block width so that the slew rate is exactly met.

    in:
        br     - real part of RF waveform (arbitrary units.)
        bi     - imag part of RF waveform (arbitrary units.)
        g      - gradient waveform (arbitrary units.)
        dt     - time-delta waveform (arbitrary units.)
        n      - number of points in b, g and dt.
        smax   - max slew rate in gradient units per time unit.
        loc    - current point to adjust (1 <= loc < n).
        indfix - array for fixing durations

    out:
        br,bi,g,dt,indfix are all adjusted.

    note: recursive!
    ------------------------------------------------------------ */
{
    double stretch;
    double slew;

    /* Do some soft error checking. */
    if (loc < 1 || loc >= n) {
        printf("Error: adjustslewindfix() - loc must be positive.\n");
        return;
    }

    if (g[loc-1] > g[loc]) { /* Downward slope Reduce height of higher */
                             /* point and work backwards.              */
    
        stretch = stretchslew(g[loc-1],g[loc],dt[loc-1],dt[loc],smax);
        g[loc-1]  /= stretch;
        br[loc-1] /= stretch;
        bi[loc-1] /= stretch;
        dt[loc-1] *= stretch;
        indfix[loc-1] = 1;

        if (loc > 1) {
            slew = (g[loc-1]-g[loc-2]) / (dt[loc-1]+dt[loc-2]) / 0.5;
            if (fabs(slew) > smax)
                adjustslewindfix(br,bi,g,dt,n,smax,loc-1,recursecount+1,indfix);
        }
    } else { /* Upward slope. Just reduce height of higher point. */
    
        stretch = stretchslew(g[loc],g[loc-1],dt[loc],dt[loc-1],smax);
        g[loc]  /= stretch;
        br[loc] /= stretch;
        bi[loc] /= stretch;
        dt[loc] *= stretch;
        indfix[loc] = 1;
    }
}


void slewcheck(double *br, double *bi, double *g, double *dt, long n, double smax)
/* ------------------------------------------------------------
    Adjust pulse so that slew rates are not violated

    in:
        br     - real part of RF waveform (arbitrary units.)
        bi     - imag part RF waveform (arbitrary units.)
        g      - gradient waveform (arbitrary units.)
        dt     - time-delta waveform (arbitrary units.)
        n      - number of points in b, g and dt.
        smax   - max slew rate in gradient units per time unit.

    out:
       br,bi,g,dt are all adjusted.
   ------------------------------------------------------------ */
{
    double slew;

    for (long count = 1; count < n; count++) {
        slew = (g[count]-g[count-1]) / (dt[count]+dt[count-1]) / 0.5;
        if (fabs(slew) > smax)
            adjustslew(br,bi,g,dt,n,smax,count,0);
    }
}


void slewcheckindfix(double *br, double *bi, double *g, double *dt, long n, double smax, long *indfix)
/* ------------------------------------------------------------
    Adjust pulse so that slew rates are not violated

    in:
        br     - real part of RF waveform (arbitrary units)
        bi     - imag part RF waveform (arbitrary units)
        g      - gradient waveform (arbitrary units)
        dt     - time-delta waveform (arbitrary units)
        n      - number of points in b, g and dt
        smax   - max slew rate in gradient units per time unit
        indfix - array for fixing durations

    out:
        br,bi,g,dt,indfix are all adjusted.
   ------------------------------------------------------------ */
{
    double slew;

    for (long count = 1; count < n; count++) {
        slew = (g[count]-g[count-1]) / (dt[count]+dt[count-1]) / 0.5;
        if (fabs(slew) > smax)
            adjustslewindfix(br,bi,g,dt,n,smax,count,0,indfix);
    }
}


void gmaxcheck(double *br, double *bi, double *g, double *dt,
               long n, double gmax, long *indfix)
/* ------------------------------------------------------------
    Adjust pulse so that gmax is not violated

    in:
        br     - real part of RF waveform (arbitrary units.)
        bi     - imag part RF waveform (arbitrary units.)
        g      - gradient waveform (arbitrary units.)
        dt     - time-delta waveform (arbitrary units.)
        n      - number of points in b, g and dt.
        gmax   - max gradient value (units of g).
        indfix - array for fixing durations

    out:
        br,bi,g,dt,indfix are all adjusted.
   ------------------------------------------------------------ */
{
    double stretch;

    for (long count = 0; count < n; count++) {
        if (g[count] > gmax) {           /* If gmax is exceeded...      */
            stretch   = g[count] / gmax; /* ...simply maximize gradient */
            br[count] /= stretch;
            bi[count] /= stretch;
             g[count] /= stretch;
            dt[count] *= stretch;
            indfix[count] = 1;
        }
    }
}


void calculateoffcenterphase(double *g, double offset, long nrf, long rfup, long gup, double gamma, double *phase)
    /* ------------------------------------------------------------
    Function calculates off-center slice RF phase waveform, phase, for an
    arbitrary gradient, g.

    Useful for calculating variable phase variation due to off-center
    VERSE pulse.

    in:
        g      - gradient waveform (arbitrary units)
        offset - slice offset (same units of distance as g)
        nrf    - number of points in rf (allows different gradient and rf
                 raster)
        rfup   - rf raster time (µs)
        gup    - gradient raster time (µs)
        gamma  - gyromagnetic ratios (Hz per same unit of field strength as g)

    out:
        phase  - slice offset phase waveform (rad)

    notes: - no error checking is done: gup should be a multiple of rfup.
           - assumes the pulse is symmetrical with repect to setting the
             center phase to 0
           - if gup>rfup, assumes g is constant during each gradient block
           - example units: g (T/m) , offset (m) , rfup (µs), gamma (Hz/T)
             or:            g (G/cm), offset (cm), rfup (µs), gamma (Hz/G)
             or:            g (mT/m), offset (mm), rfup (µs), gamma (MHz/T)
           - the output phase is not phase-wrapped
    ------------------------------------------------------------ */
{
    long M, resd, countrf, countg;
    double off, centerphase;

    off  = 2*M_PI * gamma * offset * (double)rfup*1e-6; /* precalculate (rad*distance per units of g) */
    M    = (long) ceil((double)nrf/2.0) -  1;           /* pulse centre */
    resd = gup / rfup;                                  /* gradient to RF raster difference */

    countg = 0; /* gradient counter */
    for (countrf = 0; countrf < nrf; countrf++) {

        /* Offset phase (rad) (halfway through each rfup) */
        if ((countrf%resd)==0 && countrf>0) { /* for new gradient period */

            countg++; /* countg deals with rf and g having different resolutions */
            phase[countrf] = phase[countrf-1] + off*(g[countg-1]+g[countg])/2;

        } else { /* for a current gradient period */

            if (countrf>0) phase[countrf] = phase[countrf-1] + off*g[countg];
            else           phase[countrf] = off*g[countg]/2;

        }

    }

    /* Subtract phase at pulse centre for accurate pulse phase */
    if (((nrf/resd)%2)==0) centerphase = phase[M] + off*g[M/resd]/2; /* nrf/resd is even */
    else                   centerphase = phase[M];                   /* nrf/resd is odd  */
    for (countrf = 0; countrf < nrf; countrf++) {
        phase[countrf] -= centerphase;
    }
}


long uniformresamplesize(double *dt, long n, double dtr)
/* ------------------------------------------------------------
    Function returns the number of points that will are required
    to uniformly resample an array with block-widths in *dt,
    with a uniform time block width dtr.

    in:
        dt  - n-point array of block widths for input
        n   - number of points in input
        dtr - block width for each output block

    out:
        nfr - number of output points/blocks
   ------------------------------------------------------------ */
{
    long nfr;
    double tend;

    tend = sumarrayd(dt, n);
    nfr = (long) ((tend+0.99999*dtr) / dtr); /* Round up end time */

    return nfr;
}


void resample(double *f, double *t, long n, double *fr, double *tr, long nr)
    /* ------------------------------------------------------------
    Function resamples f onto a different time step.
    The input, f, is linearly interpolated between the centers of blocks.

    in:
        f  - n-point array of heights of block-centers to resample
        t  - n-point array of block centres for input
        n  - number of points in input
        fr - n-point array of resampled heights of block-centers
        tr - n-point array of block centres for each output block
        nr - number of points in output

    out:
        fr - array of heights of output blocks.
   ------------------------------------------------------------ */
{
    double xr;

    /* Manual linear interpolation with linear extrapolation */
    long count = 1;
    for (long countr = 0; countr < nr; countr++) {

            /* Find smallest countr, where t[countr] > tr[count] */
            while (tr[countr] > t[count]) {
                if (count == (n-1))
                    break;
                count++;
            }

            /* fraction of the countr interval to interpolate into */
            xr = (tr[countr]-t[count-1]) / (t[count]-t[count-1]);

            /* interpolated value */
            fr[countr] = (1-xr)*f[count-1] + xr*f[count];
    }
}


void resampleintegrate(double *f, double dt, long n, double *fr, double dtr, long nr)
    /* ------------------------------------------------------------
    Function resamples f onto a different time step.
    The input, f, is linearly interpolated between the centers of blocks.

    Useful for externally resampling waveforms onto different raster times
    after VERSE algorithms have been used. Not used internally by
    minsarverse or mintverse.

    in:
        f  - n-point array of heights of block-centers to resample
        t  - n-point array of block centres for input
        n  - number of points in input
        fr - n-point array of resampled heights of block-centers
        tr - n-point array of block centres for each output block
        nr - number of points in output

    out:
        fr - array of heights of output blocks.
   ------------------------------------------------------------ */
{
    /* Integrate the time values */
    double *t  = (double *) malloc(n *sizeof(double)); /* array for original cumulative time */
    double *tr = (double *) malloc(nr*sizeof(double)); /* array for resampled cumulative time */

    if (!t || !tr) {
        printf("Error: resampleintegrate(): Memory allocation failed.\n");
        free(t); free(tr);
        return;
    }

    cumsumtime(dt ,n ,t );
    cumsumtime(dtr,nr,tr);

    /* Resample f */
    resample(f,t,n,fr,tr,nr);

    /* Free up memory */
    free(t);
    free(tr);
}


double calcenergy(double *br, double *bi, double *dt, long n)
/* ------------------------------------------------------------
    Calculate the energy in a waveform b, with time steps
    in dt and n points.  This is just sum(b_i^2*dt_i).

    in:
        br - real part of waveform, nominally RF (arb. units)
        bi - imag part of waveform, nominally RF (arb. units)
        dt - time-widths of sample blocks in b (arb units)
        n  - number of points in b and dt

	out:
		benergy - energy in waveform
  ------------------------------------------------------------ */
{
    double benergy = 0.0;

    for (long count = 0; count < n; count++)
        benergy += (br[count]*br[count]+bi[count]*bi[count]) * dt[count];

    return benergy;
}


void compressmax(double *br, double *bi, double *g, double *dt,
			    long n, double bmax, double gmax)
/* ------------------------------------------------------------
    Stretch or shrink the n block widths (dt) such that at each
    time step either |b| or g reaches its maximum.  Does not
    change the length of |b|, g or dt.  
    br, bi and bmax must be the same units, while g and gmax must
    also be the same units.

    in:
        br   - real part of RF waveform (arbitrary units)
        bi   - imag part of RF waveform (arbitrary units)
        g    - gradient waveform (arbitrary units.)
        dt   - time-delta waveform (arbitrary units)
        n    - number of points in b, g and dt
        bmax - max RF value (units of b)
        gmax - max gradient value (units of g)
	
	out:
		br,bi,g,dt are all updated.
   ------------------------------------------------------------ */
{
    double bfrac, gfrac;

    for (long count = 0; count < n; count++) {

        bfrac = fabs( sqrt(br[count]*br[count]+bi[count]*bi[count]) / bmax );
        gfrac = fabs( g[count] / gmax );

        if ((bfrac > 0) && (gfrac > 0)) {
            if (bfrac > gfrac)	{ /* RF will be max'd */
                br[count] /= bfrac;
                bi[count] /= bfrac;
                g[count]  /= bfrac;
                dt[count] *= bfrac;
            } else {              /* Gradient will be max'd */
                br[count] /= gfrac;
                bi[count] /= gfrac;
                g[count]  /= gfrac;
                dt[count] *= gfrac;
            }
        } else {
            #ifdef DEBUG_VERSE
            printf("compressmax(): Zero-RF or zero-gradient (%ld), ignored.\n",count);
            #endif
        }
    }
}


void minimizesar(double *br, double *bi, double *g, double *dt, long n, double gmax)
    /* ------------------------------------------------------------
    Stretch or shrink the n block widths (dt) such that at each time step 
    |b| is equal to the mean |b| of the non-zero elements.
    Does not change the length of |b|, g or dt.
    br and bi must be the same units, while g and gmax must also be the 
    same units.

    in:
        br   - real part of RF waveform (arbitrary units)
        bi   - imag part of RF waveform (arbitrary units)
        g    - gradient waveform (arbitrary units)
        dt   - time-delta waveform (arbitrary units)
        n    - number of points in b, g and dt
        gmax - max gradient value (units of g)

    out:
        br,bi,g,dt,indfix are all adjusted.
    ------------------------------------------------------------ */
{
    double bint, bn, bflat, absb1, frac;

    bint  = sumabscomplexrfarrays(br,bi,n); /* Calculate magnitude integral of RF */
    bn    = countnonzerorfarrays(br,bi,n);  /* Count non-zero RF points           */
    bflat = bint / bn;                      /* RF magnitude for minimum SAR       */

    for (long count = 0; count < n; count++) {

        absb1 = sqrt( br[count]*br[count]+bi[count]*bi[count] );
        if ( absb1 > 0 ) { /* Don't adjust zero RF points */

            frac = absb1 / bflat; /* Adjust waveform to achieve bflat */

            /* If gmax will be exceeded, simply maximize gradient */
            if ( (g[count]/frac) > gmax )
                frac = g[count] / gmax;

            br[count] /= frac;
            bi[count] /= frac;
            g[count]  /= frac;
            dt[count] *= frac;
        } else {
            #ifdef DEBUG_VERSE
            printf("minimizesar(): Zero-RF (%ld) ignored.\n",count);
            #endif
        }
    }
}


void minsarverse(double *br, double *bi, double *g, double dt,
                 long n, double gmax, double smax,
                 double *brv, double *biv, double *gv)
/* ------------------------------------------------------------
    Convert an RF/gradient pair to the minimum-SAR VERSE equivalent pair.
    This is the "Do-all" function that:
	
        1) Converts non-zero B1 to have a flat B1 amplitude.
        2) Removes gradient gmax and slew rate violations.
        3) Adjust dts to maintain constant duration.
        4) Resamples waveforms uniformly.
		
    in:
	br   - n-point real part of RF waveform (arbitrary units)
	bi   - n-point imaginary part of RF waveform (arbitrary units)
	g    - n-point gradient waveform (arbitrary units)
	dt   - time increment (arbitrary units)
	n    - number of points in b, g and dt
	gmax - max gradient value (units of g)
	smax - max gradient slew rate (units of g per unit dt)
	
    out:
	brv - output real B1 pulse array
	biv - output imaginary B1 pulse array
	gv  - output gradient array

    notes:
	    - br,bi,g,dt should be length n (not checked)
        - Points of zero RF or g will be left unchanged
   ------------------------------------------------------------ */
{

    /* Do some error checking. */
    if (n<=0) {
        printf("Error: minsarverse() - n must be greater than zero.\n");
        return;
    }
    if (dt<=0) {
        printf("Error: minsarverse() - dt must be positive.\n");
        return;
    }
    if (gmax<=0) {
        printf("Error: minsarverse() - gmax must be greater than zero.\n");
        return;
    }
    if (smax<=0) {
        printf("Error: minsarverse() - smax must be greater than zero.\n");
        return;
    }
    if (exceedsgmax(g,gmax,n)) {
        printf("Error: minsarverse() - No element of g should be > gmax!\n");
        return;
    }

    long count;         /* loop counter             */
    long niter   = 0;   /* number of iterations     */
    long maxiter = 100; /* max number of iterations */

    /*	Allocate working arrays and copy inputs to them. */
    double *brwork = (double *) malloc(n*sizeof(double)); /* working array for real part of B      */
    double *biwork = (double *) malloc(n*sizeof(double)); /* working array for imaginary part of B */
    double *gwork  = (double *) malloc(n*sizeof(double)); /* working array for gradient            */
    double *dtwork = (double *) malloc(n*sizeof(double)); /* working array for dt                  */
    double *t      = (double *) malloc(n*sizeof(double)); /* array for original cumulative time */
    double *twork  = (double *) malloc(n*sizeof(double)); /* array for new cumulative time      */
    long   *gsign  = (long   *) malloc(n*sizeof(long  )); /* array for gradient sign            */
    long   *indfix = (long   *) malloc(n*sizeof(long  )); /* array for fixing durations         */

    if (!brwork || !biwork || !gwork || !dtwork || !t || !twork || !gsign || !indfix) {
        printf("Error: minsarverse(): Memory allocation failed.\n");
        free(brwork); free(biwork); free(gwork); free(dtwork);
        return;
    }

    memcpy(brwork,br,n*sizeof(double));
    memcpy(biwork,bi,n*sizeof(double));
    memcpy(gwork ,g ,n*sizeof(double));
    arraycopy(dt,dtwork,n);

    /* Set fix counter */
    for (count = 0; count < n; count++)
        indfix[count] = 0;

    /* Calculate original timings */
    cumsumtimearray(dtwork,n,t);    /* cumulative time array */
    double T = sumarrayd(dtwork,n); /* total pulse duration  */

    /* Make gradient positive and revert at the end */
    gradientsign(g,gsign,n); /* save sign of gradient points */

    /* Generate minimum SAR B1 subject to gmax constraint */
    #ifdef DEBUG_VERSE
    printf("Calling minimizesar().\n");
    #endif
    minimizesar(brwork,biwork,gwork,dtwork,n,gmax);

    /* Shrink non-adjusted points to maintain a fixed duration */
    #ifdef DEBUG_VERSE
    printf("Calling adjusttime().\n");
    #endif
    adjusttimefixgmax(brwork,biwork,gwork,dtwork,n,T,indfix,gmax);

    while (niter < maxiter) {
        
        /* Reset fix checker */
        for (count = 0; count < n; count++)
            indfix[count] = 0;

        /* Fix gmax violations */
        #ifdef DEBUG_VERSE
        printf("Calling gmaxcheck().\n");
        #endif
        gmaxcheck(brwork,biwork,gwork,dtwork,n,gmax,indfix);

        /* Fix slew rate violations */
        #ifdef DEBUG_VERSE
        printf("Calling slewcheckindfix().\n");
        #endif
        slewcheckindfix(brwork,biwork,gwork,dtwork,n,smax,indfix);

        /* Check if there were any slew rate violations */
        if (sumarrayl(indfix,n) == 0)
            break; /* Finished */
        #ifdef DEBUG_VERSE
        printf("minsarverse(): niter = %ld, nfix = %ld.\n",niter,sumarrayl(indfix,n));
        #endif

        /* Shrink non-adjusted points to maintain a fixed duration */
        #ifdef DEBUG_VERSE
        printf("Calling adjusttime().\n");
        #endif
        adjusttimefixgmax(brwork,biwork,gwork,dtwork,n,T,indfix,gmax);

        niter++;
    }

    /* Revert gradient sign change */
    multiplyarray(gwork,gsign,n);

    /*	Resample with uniform sample step. 	*/
    #ifdef DEBUG_VERSE
    printf("Calling resample().\n");
    #endif
    cumsumtimearray(dtwork,n,twork); /* Calculate VERSEd cumulative time */
    resample(brwork,twork,n,brv,t,n);
    resample(biwork,twork,n,biv,t,n);
    resample(gwork ,twork,n,gv ,t,n);
    
    #ifdef DEBUG_VERSE
    printf("minsarverse() freeing memory.\n");
    #endif
    free(brwork);
    free(biwork);
    free(gwork);
    free(dtwork);
    free(t);
    free(twork);
    free(gsign);
    free(indfix);
}


void mintverse(double *br, double *bi, double *g, double dt, long n, 
			   double bmax, double gmax, double smax, double emax,
               long *nout, double **brv, double **biv, double **gv)
/* ------------------------------------------------------------
    Convert an RF/gradient pair to the minimum-time VERSE
    equivalent pair.  This is the "Do-all" function that:

    1) Compresses B1/gradient so one is always maximized.
    2) Removes gradient slew rate violations.
    3) Resamples waveforms uniformly.

    in:
        br   - n-point real part of RF waveform (arbitrary units)
        bi   - n-point imaginary part of RF waveform (arbitrary units)
        g    - n-point gradient waveform (arbitrary units)
        dt   - time increment (arbitrary units)
        n    - number of points in b, g and dt
        bmax - max RF value (units of b)
        gmax - max gradient value (units of g)
        smax - max gradient slew rate (units of g per unit dt)
        emax - max RF energy (units of b*b*dt) (-1 to not constrain)

    out:
        nout - number of output samples.
        brv  - address of output real B1 pulse array (allocated here)
        biv  - address of output imaginary B1 pulse array (allocated here)
        gv   - address of output gradient array (allocated here)

    notes:
        - Outputs are allocated with malloc.
        - br,bi,g,dt should be length n (not checked)
        - Points of zero RF or g will be left unchanged
   ------------------------------------------------------------ */
{

    if (n<=0) {
        printf("Error:  mintverse(): n must be greater than zero. \n");
        return;
    }
    if (dt<=0) {
        printf("Error:  mintverse(): dt must be greater than zero.\n");
        return;
    }
    if (bmax<=0) {
        printf("Error:  mintverse(): bmax must be greater than zero. \n");
        return;
    }
    if (gmax<=0) {
        printf("Error:  mintverse(): gmax must be greater than zero. \n");
        return;
    }
    if (smax<=0) {
        printf("Error:  mintverse(): smax must be greater than zero. \n");
        return;
    }

    double *brwork = (double *) malloc(n*sizeof(double)); /* working array for real part of B */
    double *biwork = (double *) malloc(n*sizeof(double)); /* working array for imaginary part of B */
    double *gwork  = (double *) malloc(n*sizeof(double)); /* working array for gradient */
    double *dtwork = (double *) malloc(n*sizeof(double)); /* working array for dt */
    double *twork  = (double *) malloc(n*sizeof(double)); /* aray for VERSEd cumulative time */
    long   *gsign  = (long *)   malloc(n*sizeof(long));   /* array for gradient sign */
    long   *indfix = (long *)   malloc(n*sizeof(long));   /* array for fixing durations         */
    double *dtout; /* array for output dt              */
    double *tout;  /* array for output cumulative time */
    double T;      /* rounded up output total pulse duration */
    long count;    /* for loops */

    if (!brwork || !biwork || !gwork || !dtwork || !twork || !gsign || !indfix) {
        printf("Error: mintverse(): Memory allocation failed.\n");
        free(brwork); free(biwork); free(gwork); free(dtwork); free(twork); free(gsign); free(indfix);
        return;
    }

    /* Set fix counter */
    for (count = 0; count < n; count++)
        indfix[count] = 0;

    double emaxratio = 0.98; /* minumum energy/emax for stopping iteration */
    double bmaxh = bmax;     /* bmax that gives energy higher than max. */
    double bmaxl = 0.0;      /* bmax that gives energy lower than max. */
    double bmaxc = bmax;     /* current bmax. */
    double benergy;          /* current energy */
    long niter = 0;		     /* number of iterations 	*/
    long maxiter = 100;      /* max number of iterations */

    /* Make gradient positive and revert at the end */
    gradientsign(g,gsign,n); /* save sign of gradient points */

    while (niter < maxiter) {

        /* copy inputs to the arrays */
        memcpy(brwork,br,n*sizeof(double));
        memcpy(biwork,bi,n*sizeof(double));
        memcpy(gwork ,g ,n*sizeof(double));
        arraycopy(dt,dtwork,n);

        #ifdef DEBUG_VERSE
        printf("Calling compressmax(). \n");
        #endif
        compressmax(brwork,biwork,gwork,dtwork,n,bmaxc,gmax);

        #ifdef DEBUG_VERSE
        printf("Calling slewcheck(). \n");
        #endif
        slewcheck(brwork,biwork,gwork,dtwork,n,smax);

        if (emax > 0) { /* adjust bmaxc to achieve emax */
            
            benergy = calcenergy(brwork, biwork, dtwork, n);

            if (benergy > emax) /* reduce bmax to reduce energy */
                bmaxh = bmaxc;
            else { /* increase bmax to increase energy */
                if (benergy/emax > emaxratio) /* close enough */
                    break; /* finished */
                bmaxl = bmaxc;
            }
            bmaxc = bmaxl + (bmaxh - bmaxl)/2.0;

            if (bmaxl >= bmaxh) {
                printf("Warning: mintverse() low-E limit > high-E limit.\n");
                printf("Exiting after %ld iterations.\n",niter);
                break;
            }
            if (niter >= maxiter) {
                printf("Warning: mintverse() iteration limit reacheed.\n");
                printf("Exiting after %ld iterations.\n",niter);
                break;
            }

            #ifdef DEBUG_VERSE
            printf("mintverse(): niter = %ld, bmax: high/low/next = %g/%g/%g\n",niter,bmaxh,bmaxl,bmaxc);
            #endif

        } else
            break; /* no energy constraint, so done! */

        niter++;
    }

    /* Revert gradient sign change */
    multiplyarray(gwork,gsign,n);

    /* Calculate length of, and allocate output vectors */
    *nout = uniformresamplesize(dtwork,n,dt);
    *brv  = (double *) malloc(*nout*sizeof(double));
    *biv  = (double *) malloc(*nout*sizeof(double));
    *gv   = (double *) malloc(*nout*sizeof(double));
    dtout = (double *) malloc(*nout*sizeof(double));
    tout  = (double *) malloc(*nout*sizeof(double));

    if (!*brv || !*biv || !*gv || !dtout || !tout) {
        printf("Error: mintverse(): Output memory allocation failed.\n");
        free(brwork); free(biwork); free(gwork); free(dtwork); free(twork); free(gsign); free(indfix);
        free(*brv); free(*biv); free(*gv); free(dtout); free(tout);
        return;
    }

    arraycopy(dt,dtout,*nout);

    /* Adjust time-points to deal with sum(dtwork) < dt*nout */
    #ifdef DEBUG_VERSE
    printf("Calling adjusttime(). \n");
    #endif
    T = *nout * dt;
    adjusttime(brwork,biwork,gwork,dtwork,n,T,indfix);

    /* Resample with uniform sample step */
    #ifdef DEBUG_VERSE
    printf("Calling resample(). \n");
    #endif
    cumsumtimearray(dtwork,n    ,twork); /* Calculate VERSEd cumulative time */
    cumsumtimearray(dtout ,*nout,tout ); /* Calculate output cumulative time */
    resample(brwork,twork,n,*brv,tout,*nout);
    resample(biwork,twork,n,*biv,tout,*nout);
    resample(gwork ,twork,n,*gv ,tout,*nout);

    #ifdef DEBUG_VERSE
    printf("mintverse() freeing memory. \n");
    #endif
    free(brwork);
    free(biwork);
    free(gwork);
    free(dtwork);
    free(twork);
    free(gsign);
    free(indfix);
    free(dtout);
    free(tout);
}
