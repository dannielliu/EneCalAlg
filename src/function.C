#include "function.h"
#include "TMath.h"
#include <vector>
#include <iostream>
std::vector<double> px1,py1,pz1,px2,py2,pz2;
double m0;
double mparticle,mparticle2,mparticle3,mpartcle4;
double factor2,factor3,factor4;
double sigma;

double maxlikelihood(double *x,double *par)
{ 
  double logl=0;
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  //double sigma=0.023;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
   for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    px=px1.at(i)+px2.at(i);
    py=py1.at(i)+py2.at(i);
    pz=pz1.at(i)+pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-x[0]*x[0]*(px*px+py*py+pz*pz));
    logl += (minv-m0)*(minv-m0)/(sigma*sigma)+TMath::Log(2.*TMath::Pi())+2.*TMath::Log(sigma);//here sigma is just a const
  }
  return logl;
} 

double weight,width;
double maxlikelihood1(double *x,double *par)
{
  double logl=0;
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  //double sigma=0.023;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
   for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    px=px1.at(i)+px2.at(i);
    py=py1.at(i)+py2.at(i);
    pz=pz1.at(i)+pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-x[0]*x[0]*(px*px+py*py+pz*pz));
    logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
            /(TMath::Sqrt(2*TMath::Pi()))/sigma*x[1] +(1.-x[1])/width);//here sigma is just a const
  }
  return logl;
}
double maxlikelihood1_1(double *x,double *par)
{
  double logl=0;
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  //double sigma=0.023;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    px=px1.at(i)+px2.at(i);
    py=py1.at(i)+py2.at(i);
    pz=pz1.at(i)+pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-x[0]*x[0]*(px*px+py*py+pz*pz));
    logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
            /(TMath::Sqrt(2*TMath::Pi()))/sigma*weight +(1.-weight)/width);//here sigma is just a const
  }
  return logl;
}

double maxlikelihood2(double *x,double *par)
{
  double logl=0;
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  //double sigma=0.0068;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle2*mparticle2+factor2*factor2*p2*p2);
    px=x[0]*px1.at(i)+factor2*px2.at(i);
    py=x[0]*py1.at(i)+factor2*py2.at(i);
    pz=x[0]*pz1.at(i)+factor2*pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-(px*px+py*py+pz*pz));
    logl += (minv-m0)*(minv-m0)/(sigma*sigma)+2*TMath::Log(2*TMath::Pi());//here sigma is just a const
  }
  return logl;
}
double maxlikelihood2_0(double *x,double *par)//two dimension
{
  double logl=0;
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  //double sigma=0.0068;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle2*mparticle2+factor2*factor2*p2*p2);
    px=x[0]*px1.at(i)+factor2*px2.at(i);
    py=x[0]*py1.at(i)+factor2*py2.at(i);
    pz=x[0]*pz1.at(i)+factor2*pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-(px*px+py*py+pz*pz));
    //logl += (minv-m0)*(minv-m0)/(sigma*sigma)+2*TMath::Log(2*TMath::Pi());//here sigma is just a const
    logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
            /(TMath::Sqrt(2*TMath::Pi())*sigma)*x[1] +(1.-x[1])/width);//here sigma is just a const
  }
  return logl;
}
double maxlikelihood2_1(double *x,double *par)//one dimension
{ 
  double logl=0;
  double p1,p2,px,py,pz;
  double e1,e2;
  double minv;
  //double sigma=0.0068;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
   for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle2*mparticle2+factor2*factor2*p2*p2);
    px=x[0]*px1.at(i)+factor2*px2.at(i);
    py=x[0]*py1.at(i)+factor2*py2.at(i);
    pz=x[0]*pz1.at(i)+factor2*pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-(px*px+py*py+pz*pz));
    //logl += (minv-m0)*(minv-m0)/(sigma*sigma)+2*TMath::Log(2*TMath::Pi());//here sigma is just a const
    logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
            /(TMath::Sqrt(2*TMath::Pi())*sigma)*weight +(1.-weight)/width);//here sigma is just a const
  }
  return logl;
}

 std::vector<double> px3,py3,pz3,px4,py4,pz4;
 std::vector<double> le1,le2;
double maxlikelihood4(double *x,double *par)
{
  double logl=0;
  double p1,p2,px,py,pz,p3,p4;
  double e1,e2,e3,e4;
  double minv;
  double mjpsi;
  //double sigma=0.03;// from large statistic, and fit it get sigma
  //double sigma=0.003;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    e3=le1.at(i);
    e4=le2.at(i);
    px=x[0]*(px1.at(i)+px2.at(i))+px3.at(i)+px4.at(i);
    py=x[0]*(py1.at(i)+py2.at(i))+py3.at(i)+py4.at(i);
    pz=x[0]*(pz1.at(i)+pz2.at(i))+pz3.at(i)+pz4.at(i);
    mjpsi=TMath::Sqrt((e3+e4)*(e3+e4)
          -(px3.at(i)+px4.at(i))*(px3.at(i)+px4.at(i))
	  -(py3.at(i)+py4.at(i))*(py3.at(i)+py4.at(i))
	  -(pz3.at(i)+pz4.at(i))*(pz3.at(i)+pz4.at(i)));
    minv =TMath::Sqrt((e1+e2+e3+e4)*(e1+e2+e3+e4)-(px*px+py*py+pz*pz))
          -mjpsi+3.096916;
    logl += (minv-m0)*(minv-m0)/(sigma*sigma)+TMath::Log(2*TMath::Pi())+2*TMath::Log(sigma);
  }
  return logl;
}
double maxlikelihood4_0(double *x,double *par)//2D function
{
  double logl=0;
  double p1,p2,px,py,pz,p3,p4;
  double e1,e2,e3,e4;
  double minv;
  double mjpsi;
  //double sigma=0.03;// from large statistic, and fit it get sigma
  //double sigma=0.003;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    e3=le1.at(i);
    e4=le2.at(i);
    px=x[0]*(px1.at(i)+px2.at(i))+px3.at(i)+px4.at(i);
    py=x[0]*(py1.at(i)+py2.at(i))+py3.at(i)+py4.at(i);
    pz=x[0]*(pz1.at(i)+pz2.at(i))+pz3.at(i)+pz4.at(i);
    mjpsi=TMath::Sqrt((e3+e4)*(e3+e4)
          -(px3.at(i)+px4.at(i))*(px3.at(i)+px4.at(i))
	  -(py3.at(i)+py4.at(i))*(py3.at(i)+py4.at(i))
	  -(pz3.at(i)+pz4.at(i))*(pz3.at(i)+pz4.at(i)));
    minv =TMath::Sqrt((e1+e2+e3+e4)*(e1+e2+e3+e4)-(px*px+py*py+pz*pz))
          -mjpsi+3.096916;
    logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
            /(TMath::Sqrt(2*TMath::Pi()))/sigma*x[1] +(1.-x[1])/width);
  }
  return logl;
}
double maxlikelihood4_1(double *x,double *par)//fix signal weight, 1D function
{
  double logl=0;
  double p1,p2,px,py,pz,p3,p4;
  double e1,e2,e3,e4;
  double minv;
  double mjpsi;
  //double sigma=0.03;// from large statistic, and fit it get sigma
  //double sigma=0.003;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    e3=le1.at(i);
    e4=le2.at(i);
    px=x[0]*(px1.at(i)+px2.at(i))+px3.at(i)+px4.at(i);
    py=x[0]*(py1.at(i)+py2.at(i))+py3.at(i)+py4.at(i);
    pz=x[0]*(pz1.at(i)+pz2.at(i))+pz3.at(i)+pz4.at(i);
    mjpsi=TMath::Sqrt((e3+e4)*(e3+e4)
          -(px3.at(i)+px4.at(i))*(px3.at(i)+px4.at(i))
	  -(py3.at(i)+py4.at(i))*(py3.at(i)+py4.at(i))
	  -(pz3.at(i)+pz4.at(i))*(pz3.at(i)+pz4.at(i)));
    minv =TMath::Sqrt((e1+e2+e3+e4)*(e1+e2+e3+e4)-(px*px+py*py+pz*pz))
          -mjpsi+3.096916;
    //logl += (minv-m0)*(minv-m0)/(sigma*sigma)+TMath::Log(2*TMath::Pi())+2*TMath::Log(sigma);
    logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
            /(TMath::Sqrt(2*TMath::Pi()))/sigma*weight +(1.-weight)/width);
  //std::cout<<mjpsi<<"\t"<<minv<<"\t"<< -(minv-m0)*(minv-m0)/(2.*sigma*sigma)<<"\t"<<logl<<std::endl;
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
