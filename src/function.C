#include "function.h"
#include "TMath.h"

double BreitWigner(double *x, double *par)
{
  return par[0]*par[2]/(2*TMath::Pi()*(TMath::Power(x[0]-par[1],2)+TMath::Power(par[2]/2.,2)));
}

double line(double *x, double *par)
{
  return par[0]+par[1]*x[0];
}

double line2(double *x, double *par)
{ 
  return (x[0]-par[0])*par[1];
}

double parabola(double *x, double *par)
{
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

double fitfun(double *x, double *par)
{
  return BreitWigner(x,par)+line(x,&par[3]);
}

double GausLineBack(double *x,double *par)
{
  //TF1 f("gaus","gaus");
  
  return par[0]*TMath::Gaus(*x,par[1],par[2])+line(x,&par[3]);
}
