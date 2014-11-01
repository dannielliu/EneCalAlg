#ifndef myfunction_h
#define myfunction_h
#include "TMath.h"

double BreitWigner(double *x, double *par);

double line(double *x, double *par);

double line2(double *x, double *par);

double parabola(double *x, double *par);

double fitfun(double *x, double *par);

double GausLineBack(double *x,double *par);

#endif
