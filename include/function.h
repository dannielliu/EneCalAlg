#ifndef myfunction_h
#define myfunction_h
#include "TMath.h"

double BreitWigner(double *x, double *par);

double line(double *x, double *par);

double line2(double *x, double *par);

double parabola(double *x, double *par);

double fitfun(double *x, double *par);

double GausLineBack(double *x,double *par);

// gaussian function, ignore the background
double maxlikelihood(double *x,double *par);

// function with background
double maxlikelihood1(double *x,double *par);
double maxlikelihood1_1(double *x,double *par);
double maxlikelihood1_3(double *x,double *par);

// for two diffrent particles end
double maxlikelihood2(double *x,double *par);//pure gaussian
double maxlikelihood2_0(double *x,double *par);//2D, background variable
double maxlikelihood2_1(double *x,double *par);//1D, background fixed

// for two pairs of particles end
double maxlikelihood4(double *x,double *par);
double maxlikelihood4_0(double *x,double *par);//2D function
double maxlikelihood4_1(double *x,double *par);//1D

// for three particle
double maxlikelihood3(double *x,double *par);
double maxlikelihood3_0(double *x,double *par);//two dimension
double maxlikelihood3_1(double *x,double *par);//one dimension


#endif
