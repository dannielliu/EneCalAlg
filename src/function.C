#include "function.h"
#include "TMath.h"
#include <vector>
#include <iostream>
extern std::vector<double> px1,py1,pz1,px2,py2,pz2;
extern double m0;
extern double me;

double maxlikelihood(double *x,double *par)
{
  double logl=0;
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(me*me+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(me*me+x[0]*x[0]*p2*p2);
    px=px1.at(i)+px2.at(i);
    py=py1.at(i)+py2.at(i);
    pz=pz1.at(i)+pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-x[0]*x[0]*(px*px+py*py+pz*pz));
    logl += 2*(minv-m0)*(minv-m0);//+2*TMath::Log(2*TMath::Pi());//here sigma is just a const
  }
  return logl;
}

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
