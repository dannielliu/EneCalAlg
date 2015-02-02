#ifndef myfunction_h
#define myfunction_h
#include "TMath.h"

double BreitWigner(double *x, double *par);

double line(double *x, double *par);

double line2(double *x, double *par);

double parabola(double *x, double *par);

double fitfun(double *x, double *par);

double GausLineBack(double *x,double *par);

double CalInvMass(double m1, double px1, double py1, double pz1,
                  double m2, double px2, double py2, double pz2,
                  int n=0, const double *x=0, const double *par=0);

double CalMom(double px, double py, double pz);

double CalEne(double m, double px,double py,double pz,double factor=1.);

#endif
