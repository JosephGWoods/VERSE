#ifndef VERSE_H
#define VERSE_H

#ifdef __cplusplus
extern "C" {
#endif

void mintverse(
    double *br, double *bi, double *g, double dt, long n,
    double bmax, double gmax, double smax, double emax,
    long *nout, double **brv, double **biv, double **gv);

void minsarverse(
    double *br, double *bi, double *g, double dt,
    long n, double gmax, double smax,
    double *brv, double *biv, double *gv);

void calculateoffcenterphase(
    double *g, double offset, long nrf,
    long rfup, long gup, double gamma, double *phase);

#ifdef __cplusplus
}
#endif

#endif // VERSE_H