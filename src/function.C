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

double weight,width,slope=0.;
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
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*x[1] +(1.-x[1])*(slope*(minv-m0)+1./width));//here sigma is just a const
    //logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
    //        /(TMath::Sqrt(2*TMath::Pi()))/sigma*x[1] +(1.-x[1])*(slope*(minv-m0)+1./width));//here sigma is just a const
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
  for(int i=0; i<px1.size()-1; i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    px=px1.at(i)+px2.at(i);
    py=py1.at(i)+py2.at(i);
    pz=pz1.at(i)+pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-x[0]*x[0]*(px*px+py*py+pz*pz));
    //logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
    //        /(TMath::Sqrt(2*TMath::Pi()))/sigma*weight +(1.-weight)*(slope*(minv-m0)+1./width));
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*weight +(1.-weight)*(slope*(minv-m0)+1./width));
    //std::cout<<"factor: "<<x[0]<<"\tminv: "<<minv<<"\tsigma:"<<sigma<<"\tlogl:"<<logl<<std::endl;
  }
  return logl;
}
double maxlikelihood1_3(double *x,double *par)
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
            /(TMath::Sqrt(2*TMath::Pi())*sigma)*x[1] +(1.-x[1])*(x[2]*(minv-m0)+1./width));
    //x[0] is the factor, x[1] is signal weight, x[2] is background slope.
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
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*weight +(1.-weight)/width);//here sigma is just a const
    //logl += (minv-m0)*(minv-m0)/(sigma*sigma)+2*TMath::Log(2*TMath::Pi());//here sigma is just a const
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
    e2=TMath::Sqrt(mparticle2*mparticle2+x[0]*x[0]*p2*p2);
    px=x[0]*px1.at(i)+x[0]*px2.at(i);
    py=x[0]*py1.at(i)+x[0]*py2.at(i);
    pz=x[0]*pz1.at(i)+x[0]*pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-(px*px+py*py+pz*pz));
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*x[1] +(1.-x[1])/width);//here sigma is just a const
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
    e2=TMath::Sqrt(mparticle2*mparticle2+x[0]*x[0]*p2*p2);
    px=x[0]*px1.at(i)+x[0]*px2.at(i);
    py=x[0]*py1.at(i)+x[0]*py2.at(i);
    pz=x[0]*pz1.at(i)+x[0]*pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-(px*px+py*py+pz*pz));
    //logl += (minv-m0)*(minv-m0)/(sigma*sigma)+2*TMath::Log(2*TMath::Pi());//here sigma is just a const
    //logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
    //        /(TMath::Sqrt(2*TMath::Pi())*sigma)*weight +(1.-weight)/width);//here sigma is just a const
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*weight +(1.-weight)/width);//here sigma is just a const
  }
  return logl;
}
double maxlikelihood2_2(double *x,double *par)//two factor in different momentum range
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
    double f1=x[0],f2=x[0];
    if (p1>par[0]) f1=x[1];
    if (p2>par[0]) f2=x[1];
    e1=TMath::Sqrt(mparticle*mparticle  +f1*f1*p1*p1);
    e2=TMath::Sqrt(mparticle2*mparticle2+f2*f2*p2*p2);
    px=f1*px1.at(i)+f2*px2.at(i);
    py=f1*py1.at(i)+f2*py2.at(i);
    pz=f1*pz1.at(i)+f2*pz2.at(i);
    minv = TMath::Sqrt((e1+e2)*(e1+e2)-(px*px+py*py+pz*pz));
    //logl += (minv-m0)*(minv-m0)/(sigma*sigma)+2*TMath::Log(2*TMath::Pi());//here sigma is just a const
    //logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
    //        /(TMath::Sqrt(2*TMath::Pi())*sigma)*weight +(1.-weight)/width);//here sigma is just a const
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*weight +(1.-weight)/width);//here sigma is just a const
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
    //logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
    //        /(TMath::Sqrt(2*TMath::Pi()))/sigma*x[1] +(1.-x[1])/width);
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*x[1] +(1.-x[1])/width);
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
    //logl += -2*TMath::Log(TMath::Exp(-(minv-m0)*(minv-m0)/(2.*sigma*sigma))
     //       /(TMath::Sqrt(2*TMath::Pi()))/sigma*weight +(1.-weight)/width);
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*weight +(1.-weight)/width);
  //std::cout<<mjpsi<<"\t"<<minv<<"\t"<< -(minv-m0)*(minv-m0)/(2.*sigma*sigma)<<"\t"<<logl<<std::endl;
  }
  return logl;
}

double maxlikelihood3(double *x,double *par)
{
  double logl=0;
  double p1,p2,p3,px,py,pz;
  double e1,e2,e3;
  double minv;
  //double sigma=0.0068;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    p3=TMath::Sqrt(px3.at(i)*px3.at(i)+py3.at(i)*py3.at(i)+pz3.at(i)*pz3.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    e3=TMath::Sqrt(mparticle2*mparticle2+x[0]*x[0]*p3*p3);
    px=x[0]*px1.at(i)+x[0]*px2.at(i)+x[0]*px3.at(i);
    py=x[0]*py1.at(i)+x[0]*py2.at(i)+x[0]*py3.at(i);
    pz=x[0]*pz1.at(i)+x[0]*pz2.at(i)+x[0]*pz3.at(i);
    minv = TMath::Sqrt((e1+e2+e3)*(e1+e2+e3)-(px*px+py*py+pz*pz));
    logl += (minv-m0)*(minv-m0)/(sigma*sigma)+2*TMath::Log(2*TMath::Pi());//here sigma is just a const
  }
  return logl;
}
double maxlikelihood3_0(double *x,double *par)//two dimension
{
  double logl=0;
  double p1,p2,p3,px,py,pz;
  double e1,e2,e3;
  double minv;
  //double sigma=0.0068;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    p3=TMath::Sqrt(px3.at(i)*px3.at(i)+py3.at(i)*py3.at(i)+pz3.at(i)*pz3.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    e3=TMath::Sqrt(mparticle2*mparticle2+x[0]*x[0]*p3*p3);
    px=x[0]*px1.at(i)+x[0]*px2.at(i)+x[0]*px3.at(i);
    py=x[0]*py1.at(i)+x[0]*py2.at(i)+x[0]*py3.at(i);
    pz=x[0]*pz1.at(i)+x[0]*pz2.at(i)+x[0]*pz3.at(i);
    minv = TMath::Sqrt((e1+e2+e3)*(e1+e2+e3)-(px*px+py*py+pz*pz));
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*x[1] +(1.-x[1])/width);//here sigma is just a const
  }
  return logl;
}
double maxlikelihood3_1(double *x,double *par)//one dimension
{ 
  double logl=0;
  double p1,p2,p3,px,py,pz;
  double e1,e2,e3;
  double minv;
  //double sigma=0.0068;// from large statistic, and fit it get sigma
  //std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
  for(int i=0; i<px1.size(); i++){
    p1=TMath::Sqrt(px1.at(i)*px1.at(i)+py1.at(i)*py1.at(i)+pz1.at(i)*pz1.at(i));
    p2=TMath::Sqrt(px2.at(i)*px2.at(i)+py2.at(i)*py2.at(i)+pz2.at(i)*pz2.at(i));
    p3=TMath::Sqrt(px3.at(i)*px3.at(i)+py3.at(i)*py3.at(i)+pz3.at(i)*pz3.at(i));
    e1=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p1*p1);
    e2=TMath::Sqrt(mparticle*mparticle+x[0]*x[0]*p2*p2);
    e3=TMath::Sqrt(mparticle2*mparticle2+x[0]*x[0]*p3*p3);
    px=x[0]*px1.at(i)+x[0]*px2.at(i)+x[0]*px3.at(i);
    py=x[0]*py1.at(i)+x[0]*py2.at(i)+x[0]*py3.at(i);
    pz=x[0]*pz1.at(i)+x[0]*pz2.at(i)+x[0]*pz3.at(i);
    minv = TMath::Sqrt((e1+e2+e3)*(e1+e2+e3)-(px*px+py*py+pz*pz));
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true)*weight +(1.-weight)/width);//here sigma is just a const
  }
  return logl;
}

std::vector<double> mass1;
double maxlikelihood0(double *x, double *par)
{
  double ret=0;
  for(int i=0; i<mass1.size(); i++){
    ret += -2*TMath::Log(TMath::Gaus(x[0]*mass1.at(i),m0,sigma,true));
  }
  return ret;
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
  
  return par[0]*TMath::Gaus(*x,par[1],par[2],true)+line(x,&par[3]);
}

double GausLineBack1(double *x,double *par)
{
  //TF1 f("gaus","gaus");
  
  return par[0]*TMath::Gaus(*x,par[1],par[2],true)+line(x,&par[3]);
}

double maxlikelihood2_n(int n,const double *x,const double *par)//two factor in different momentum range
{ 
  double logl=0;
  double minv;
  //std::cout<<"par0: "<<par[0]<<", par1: "<<par[1]<<std::endl;
  for(int i=0; i<px1.size(); i++){
    minv = CalInvMass(mparticle,px1.at(i),py1.at(i),pz1.at(i),
                      mparticle,px2.at(i),py2.at(i),pz2.at(i),
		  n,x,par);
    logl += -2*TMath::Log(TMath::Gaus(minv,m0,sigma,true));//*weight +(1.-weight)/width);
    //std::cout<<"inv mass: "<<minv<<" logl value: "<<logl<<std::endl;
  }
  return logl;
}

double MultiDimFunction::DoEvalPar(const double *x, const double *p) const
{
  return maxlikelihood2_n(ndim,x,pars);
}

unsigned int MultiDimFunction::NDim() const
{
  return ndim;
}

unsigned int MultiDimFunction::NPar() const
{
  return npar;
}

double MultiDimFunction::DoParameterDerivative(const double *x, const double *p,unsigned ipar) const
{
  return 0;
}

double MultiDimFunction::GetMinimum(double *xmin, double *xmax)
{
  return fValMin;
}

void MultiDimFunction::GetMinimumX(double *x)
{
  double xx[ndim];
  double xxmin[ndim];
  double xxlowedge[ndim];
  for (int i=0; i<ndim; i++){
    xxmin[i]=fXmin;
    xxlowedge[i]=fXmin;
  }
  double fval; 
  int Npx=2;
  double dx=(fXmax-fXmin)/Npx;
  for (int i=0; i<ndim; i++) xx[i]=fXmin;
  double zzmin=DoEvalPar(xx,pars);
  int minindex=0;
  //std::cout<<"initial min value: "<<zzmin<<std::endl;
  int itern =(int)(TMath::Log(fNpx)/TMath::Log(Npx))+10;
  for (int iteri=0; iteri<itern;iteri++){
    for (int i=0; i<TMath::Power(Npx,ndim); i++){
      std::cout<<"point: (";
      for (int idim=0; idim<ndim; idim++){
        xx[idim] = xxlowedge[idim] + ((i%(int)TMath::Power(Npx,idim+1))/(int)TMath::Power(Npx,idim)+0.5)*dx;// it is dangerous, power may give a wrong int!!
      std::cout<<xx[idim]<<", ";
      }
      fval=DoEvalPar(xx,pars);
      std::cout<<"), fvalue: "<<fval<<std::endl;
      if (fval<zzmin){
        zzmin = fval;
        for (int idim=0; idim<ndim; idim++) xxmin[idim]=xx[idim];
        minindex = i;
      }
    }
    for (int idim=0; idim<ndim; idim++){
      xxlowedge[idim]= xxlowedge[idim] + ((minindex%(int)TMath::Power(Npx,idim+1))/(int)TMath::Power(Npx,idim))*dx;// it is dangerous, power may give a wrong int!!!
    }
    dx = dx/Npx;
  }

  for (int i=0; i<ndim; i++){
    x[i] = TMath::Min(xxmin[i],fXmax);
  }
  fValMin = zzmin;

  return;
}

void MultiDimFunction::ShowParameters()
{
  std::cout<<"par[0]: "<<pars[0]<<" , par[1]: "<<pars[1]<<std::endl;
  return;
}

double BiasCoe(const double *x,double *par)
{ 
  // par[0]: expected mean value, par[1]: sigma, 
  // now it's not used here though
  double bias=0;
  double minv;
  double M3;
  double M2;
  int vsize=px1.size();
  for(int i=0; i<vsize; i++){
    minv = CalInvMass(mparticle,px1.at(i),py1.at(i),pz1.at(i),
                      mparticle,px2.at(i),py2.at(i),pz2.at(i),
		  1,x);
    M3 += TMath::Power(minv-par[0],3);
    M2 += TMath::Power(minv-par[0],2);
  }
  M3 = M3/vsize;
  M2 = M2/vsize;
  bias = fabs(M3/TMath::Power(M2,3./2));
  return bias;
}

//vector<double> p1,p2,theta;
//par[0]: uniform constant
//par[1]: correction factor
//x[0]: momentum of particle 1
//x[1]: momentum of particle 2
//x[2]: angle between 1,2
double distribution(double *x, double *par)
{
  double dis=0;
  double minv,e1,e2;
  double p1=par[1]*x[0];
  double p2=par[1]*x[1];
  double dm=0;
  e1=TMath::Sqrt(p1*p1+mparticle*mparticle);
  e2=TMath::Sqrt(p2*p2+mparticle*mparticle);
  minv = TMath::Sqrt((e1+e2)*(e1+e2)-p1*p1-p2*p2-2*p1*p2*cos(x[2]));
  dm = (-par[1]*par[1]/minv
        +TMath::Power(par[1],4)*x[1]*x[1]/TMath::Power(minv,3)*(e1/e2-cos(x[2]))
        +TMath::Power(par[1],4)*x[0]*x[0]/TMath::Power(minv,3)*(e2/e1-cos(x[2]))
        +(-3./(minv*minv)+1./(e1*e2))*TMath::Power(par[1],6)*TMath::Power(x[0]*x[1],2)/TMath::Power(minv,5)*(e2/e1-cos(x[2])*(e1/e2-cos(x[2])))
       )*sin(x[2]);
  dis = par[0]*TMath::Gaus(minv,m0,sigma,true)*dm;
  return dis;
}
