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

double CalInvMass(double m1, double px1, double py1, double pz1,
                  double m2, double px2, double py2, double pz2,
                  int n, const double *x, const double *par)
{
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  double f1,f2;
  
  p1=TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
  p2=TMath::Sqrt(px2*px2+py2*py2+pz2*pz2);
  if (n==0) {f1=1; f2=1;}
  else if(n==1) {f1=x[0]; f2=x[0];}
  else if(n==-2) {f1=x[0]; f2=x[1];}
  else {
    int tmpindex;
    tmpindex=(int)((p1-par[0])/(par[1]-par[0])*n);
    if (tmpindex<0) tmpindex =0;
    if (tmpindex>n-1) tmpindex=n-1;
    f1 = x[tmpindex];
    tmpindex=(int)((p2-par[0])/(par[1]-par[0])*n);
    if (tmpindex<0) tmpindex =0;
    if (tmpindex>n-1) tmpindex=n-1;
    f2 = x[tmpindex];
  }
  e1=TMath::Sqrt(m1*m1+f1*f1*p1*p1);
  e2=TMath::Sqrt(m2*m2+f2*f2*p2*p2);
  px=f1*px1+f2*px2;
  py=f1*py1+f2*py2;
  pz=f1*pz1+f2*pz2;
  minv = TMath::Sqrt((e1+e2)*(e1+e2)-(px*px+py*py+pz*pz));
  
  return minv;
}

double CalMom(double px, double py, double pz)
{
  return TMath::Sqrt(px*px+py*py+pz*pz);
}

double CalEne(double m, double px,double py,double pz,double factor)
{
  double p=CalMom(px,py,pz);
  return TMath::Sqrt(m*m+factor*factor*p*p);
}

