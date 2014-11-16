#include "TRandom.h"
#include <iostream>
#include <fstream>
#include "TMath.h"

int main(int argc,char** argv)
{
  ofstream sample;
  sample.open("sample.txt");
  double cmspsig,cmspsig1,cmspsig2,cmspbck,cmspbck1,cmspbck2,msig,mbck;
  double labpsig,labpsig1,labpsig2,labpbck,labpbck1,labpbck2;
  double cmstheta1,cmstheta2,cmsphi1,cmsphi2;
  double labtheta1,labtheta2,labphi1,labphi2;
  double betac,gammac;
  double cmsbeta;
  double pi=TMath::Pi();
  double twopi=2*pi;
  double mpeak=1.01946;
  double sigma=0.0025;
  int sigNo=300;
  int backNo=1000;
  double mk=0.493677;
  double c = 299792458;
  gRandom->SetSeed(time(0));
  //signal part,
  for(int i=0;i<sigNo;i++){
    //initial
    msig=gRandom->Gaus(mpeak,sigma);
    cmspsig = gRandom->Exp(0.5);
    betac = cmspsig/TMath::Sqrt(msig*msig+cmspsig*cmspsig);
    gammac = TMath::Sqrt(msig*msig+cmspsig*cmspsig)/msig;
    //split
    //cms
    cmspsig1 = TMath::Sqrt(msig*msig/4.0 - mk*mk);
    cmspsig2 = cmspsig1;
    cmstheta1 = gRandom->Uniform(0,pi);
    cmsphi1   = gRandom->Uniform(0,twopi);
    cmstheta2 = pi - cmstheta1;
    cmsphi2   = cmsphi1 + pi;
//std::cout<<cmstheta1<<" "<<cmstheta2<<" "<<cmsphi1<<" "<<cmsphi2<<"\n";
    // change to lab
    cmsbeta   = cmspsig1 / (msig/2.);
    labtheta1 = atan(sin(cmstheta1)/(gammac*(cos(cmstheta1)+betac/cmsbeta)));
    if(labtheta1<0) labtheta1 += pi;
    labpsig1 = cmspsig1*sin(cmstheta1)/sin(labtheta1);
    labphi1  = cmsphi1;
    labtheta2 = atan(sin(cmstheta2)/(gammac*(cos(cmstheta2)+betac/cmsbeta)));
    if(labtheta1<0) labtheta2 += pi;
    labpsig2 = cmspsig1*sin(cmstheta2)/sin(labtheta2);
    labphi2  = cmsphi2;
    double p1[3],p2[3];
    p1[0] = labpsig1*sin(labtheta1)*cos(labphi1);
    p1[1] = labpsig1*sin(labtheta1)*sin(labphi1);
    p1[2] = labpsig1*cos(labtheta1);
    p2[0] = labpsig2*sin(labtheta2)*cos(labphi2);
    p2[1] = labpsig2*sin(labtheta2)*sin(labphi2);
    p2[2] = labpsig2*cos(labtheta2);
    //psig1 = gRandom->Gaus(1.,0.5);
    //psig2 = TMath::Sqrt(mpeak*mpeak - msig*msig)-psig1;
    sample<<"\n"<<p1[0]<<"\t"<<p1[1]<<"\t"<<p1[2]<<"\t";
    sample<<"\t"<<p2[0]<<"\t"<<p2[1]<<"\t"<<p2[2]<<"\t";
  }
  //background part,
  sample<<"\n";
  for(int i=0;i<backNo;i++){
    //initial
    msig    = gRandom->Uniform(mpeak-10*sigma,mpeak+30*sigma);
    cmspsig = gRandom->Exp(0.5);
    betac   = cmspsig/TMath::Sqrt(msig*msig+cmspsig*cmspsig);
    gammac  = TMath::Sqrt(msig*msig+cmspsig*cmspsig)/msig;
    //split
    //cms
    cmspsig1 = TMath::Sqrt(msig*msig/4.0 - mk*mk);
    cmspsig2 = cmspsig1;
    cmstheta1 = gRandom->Uniform(0,pi);
    cmsphi1   = gRandom->Uniform(0,twopi);
    cmstheta2 = pi - cmstheta1;
    cmsphi2   = cmsphi1 + pi;
    // change to lab
    cmsbeta   = cmspsig1 / (msig/2.);
    if(cos(cmstheta1)+betac/cmsbeta==0){
      labtheta1 = pi/2;
    }
    else
      labtheta1 = atan(sin(cmstheta1)/(gammac*(cos(cmstheta1)+betac/cmsbeta)));
    if(labtheta1<0) labtheta1 += pi;
    labpsig1 = cmspsig1*sin(cmstheta1)/sin(labtheta1);
    labphi1  = cmsphi1;
    if(cos(cmstheta2)+betac/cmsbeta==0){
      labtheta2 = pi/2;
    }else
      labtheta2 = atan(sin(cmstheta2)/(gammac*(cos(cmstheta2)+betac/cmsbeta)));
    if(labtheta2<0) labtheta2 += pi;
    labpsig2 = cmspsig1*sin(cmstheta2)/sin(labtheta2);
    labphi2  = cmsphi2;
    double p1[3],p2[3];
    p1[0] = labpsig1*sin(labtheta1)*cos(labphi1);
    p1[1] = labpsig1*sin(labtheta1)*sin(labphi1);
    p1[2] = labpsig1*cos(labtheta1);
    p2[0] = labpsig2*sin(labtheta2)*cos(labphi2);
    p2[1] = labpsig2*sin(labtheta2)*sin(labphi2);
    p2[2] = labpsig2*cos(labtheta2);
    //psig1 = gRandom->Gaus(1.,0.5);
    //psig2 = TMath::Sqrt(mpeak*mpeak - msig*msig)-psig1;
    sample<<"\n"<<p1[0]<<"\t"<<p1[1]<<"\t"<<p1[2]<<"\t";
    sample<<"\t"<<p2[0]<<"\t"<<p2[1]<<"\t"<<p2[2]<<"\t";
  //  mbck=gRandom->Uniform(mpeak-30*sigma,mpeak-30*sigma);
  //  pbck1 = gRandom->Uniform(0,2.0);
  //  pbck2 = TMath::Sqrt(mbck);
  //  sample<<"\n"<<pbck;
  }
  sample.close();
  return 0;
}
