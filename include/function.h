#ifndef myfunction_h
#define myfunction_h
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include <vector>

double BreitWigner(double *x, double *par);

double line(double *x, double *par);

double line2(double *x, double *par);

double parabola(double *x, double *par);

double fitfun(double *x, double *par);

double GausLineBack(double *x,double *par);

double CalInvMass(double m1, double px1, double py1, double pz1,
                  double m2, double px2, double py2, double pz2,
                  int n=0, const double *x=0, const double *par=0);

double CalMom(double px, double py, double pz=0);

double CalEne(double m, double px,double py,double pz,double factor=1.);


struct Event
{
  double mpi;
  double px1;
  double py1;
  double pz1;
  double px2;
  double py2;
  double pz2;

public:
  Event() {mpi = 0.13957;}
  Event(double mass) { mpi = mass;}
  void SetMass(double mass) { mpi = mass;}
  void SetVal(double x1, double y1, double z1, double x2, double y2, double z2)
  {
    px1=x1;
    py1=y1;
    pz1=z1;
    px2=x2;
    py2=y2;
    pz2=z2;
  }
  double InvMass(double f1=1.0,double f2=1.0)
  {
    double f[2];
	f[0] = f1;
	f[1] = f2;
    return CalInvMass(mpi,px1,py1,pz1,mpi,px2,py2,pz2,-2,f);
  }
  double GetP1()
  {
    return CalMom(px1,py1,pz1);
  }
  double GetP2()
  {
    return CalMom(px2,py2,pz2);
  }
  double GetCostheta()
  {
    return (px1*px2+py1*py2+pz1*pz2)/(GetP1()*GetP2());
  }
  double GetCostheta1()
  {
    return pz1/GetP1();
  }
  double GetCostheta2()
  {
    return pz2/GetP2();
  }
  double GetPhi1()
  {
    if (py1>0) return acos(px1/sqrt(px1*px1+py1*py1));
    else return 2*TMath::Pi()-acos(px1/sqrt(px1*px1+py1*py1)); 
  }
  double GetPhi2()
  {
    if (py2>0) return acos(px2/sqrt(px2*px2+py2*py2));
    else return 2*TMath::Pi()-acos(px2/sqrt(px2*px2+py2*py2)); 
  }
};

struct Psip
{
  double mpi;
  double pipx1;
  double pipy1;
  double pipz1;
  double pipx2;
  double pipy2;
  double pipz2;
  double ml;
  double lpx1;
  double lpy1;
  double lpz1;
  double lpx2;
  double lpy2;
  double lpz2;

public:
  Psip(){
    mpi = 0.13957;
	ml = 0.000511;
  }
  void Setval(double px1,double py1,double pz1,
              double px2,double py2,double pz2,
              double px3,double py3,double pz3,
              double px4,double py4,double pz4 )
  {
    pipx1 = px1;
    pipy1 = py1;
    pipz1 = pz1;
    pipx2 = px2;
    pipy2 = py2;
    pipz2 = pz2;
    lpx1  = px3;
    lpy1  = py3;
    lpz1  = pz3;
    lpx2  = px4;
    lpy2  = py4;
    lpz2  = pz4;
  }
  void SetLeptonM(double m)
  {
    ml = m;
  }
  double InvMass(double f1=1.0,double f2=1.0,double f3=1.00,double f4=1.00)
  {
    double totpx,totpy,totpz,tote;
	double lee[2],pie[2];
	double mass;
    totpx = (f1*pipx1+f2*pipx2)+(f3*lpx1+f4*lpx2);
    totpy = (f1*pipy1+f2*pipy2)+(f3*lpy1+f4*lpy2);
    totpz = (f1*pipz1+f2*pipz2)+(f3*lpz1+f4*lpz2);
    lee[0]=CalEne(ml,f3*lpx1,f3*lpy1,f3*lpz1);
    lee[1]=CalEne(ml,f4*lpx2,f4*lpy2,f4*lpz2);
    //double massjpsi = CalInvMass(ml,lpx1,lpy1,lpz1,ml,lpx2,lpy2,lpz2);
    pie[0]=CalEne(mpi,f1*pipx1,f1*pipy1,f1*pipz1);
    pie[1]=CalEne(mpi,f2*pipx2,f2*pipy2,f2*pipz2);
    tote=lee[0]+lee[1]+pie[0]+pie[1];
    mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
    mass = mass-PsiM(f3,f4)+3.096916;
	return mass;
  } 
  double PsiM(double f3=1.00, double f4 = 1.00)
  {
    double f[2];
	f[0] = f3;
	f[1] = f4;
    return CalInvMass(ml,lpx1,lpy1,lpz1,ml,lpx2,lpy2,lpz2,-2,f);
  }
  double GetP1()
  {
    return CalMom(pipx1,pipy1,pipz1);
  }
  double GetP2()
  {
    return CalMom(pipx2,pipy2,pipz2);
  }
  double GetCostheta1()
  {
    return pipz1/GetP1();
  }
  double GetCostheta2()
  {
    return pipz2/GetP2();
  }
  double GetPhi1()
  {
    if (pipy1>0) return acos(pipx1/sqrt(pipx1*pipx1+pipy1*pipy1));
    else return 2*TMath::Pi()-acos(pipx1/sqrt(pipx1*pipx1+pipy1*pipy1)); 
  }
  double GetPhi2()
  {
    if (pipy2>0) return acos(pipx2/sqrt(pipx2*pipx2+pipy2*pipy2));
    else return 2*TMath::Pi()-acos(pipx2/sqrt(pipx2*pipx2+pipy2*pipy2)); 
  }

};

class EEto6pi
{
  double mpi;
  double px[6];
  double py[6];
  double pz[6];
public:
  EEto6pi() {mpi=0.13957;}
  ~EEto6pi() {}
  void Setval(double *ppx,double *ppy,double *ppz,
              double *mpx,double *mpy,double *mpz)
  {
    px[0] = ppx[0];
    px[1] = ppx[1];
    px[2] = ppx[2];
    px[3] = mpx[0];
    px[4] = mpx[1];
    px[5] = mpx[2];
    py[0] = ppy[0];
    py[1] = ppy[1];
    py[2] = ppy[2];
    py[3] = mpy[0];
    py[4] = mpy[1];
    py[5] = mpy[2];
    pz[0] = ppz[0];
    pz[1] = ppz[1];
    pz[2] = ppz[2];
    pz[3] = mpz[0];
    pz[4] = mpz[1];
    pz[5] = mpz[2];
	return ;
  }
  double GetP(int i=0)
  {
    return CalMom(px[i],py[i],pz[i]);
  }
  double InvMass(double f0,double f1,double f2,double f3,double f4,double f5)
  {
    double totpx,totpy,totpz,tote;
    double pie[6];
    double mass;
    totpx = f0*px[0]+f1*px[1]+f2*px[2]+f3*px[3]+f4*px[4]+f5*px[5];
    totpy = f0*py[0]+f1*py[1]+f2*py[2]+f3*py[3]+f4*py[4]+f5*py[5];
    totpz = f0*pz[0]+f1*pz[1]+f2*pz[2]+f3*pz[3]+f4*pz[4]+f5*pz[5];
    //tote  = pipe[0] +pipe[1] +pipe[2] +pime[0] +pime[1] +pime[2];
    pie[0]=TMath::Sqrt(mpi*mpi+TMath::Power(f0*GetP(0),2) );
    pie[1]=TMath::Sqrt(mpi*mpi+TMath::Power(f1*GetP(1),2) );
    pie[2]=TMath::Sqrt(mpi*mpi+TMath::Power(f2*GetP(2),2) );
    pie[3]=TMath::Sqrt(mpi*mpi+TMath::Power(f3*GetP(3),2) );
    pie[4]=TMath::Sqrt(mpi*mpi+TMath::Power(f4*GetP(4),2) );
    pie[5]=TMath::Sqrt(mpi*mpi+TMath::Power(f5*GetP(5),2) );
    tote  = pie[0] +pie[1] +pie[2] +pie[3] +pie[4] +pie[5];
    mass = TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
    return mass;
  }
  double InvMass(double *f)
  {
    return InvMass(f[0],f[1],f[2],f[3],f[4],f[5]);
  }
  double InvMass()
  {
    return InvMass(1.0,1.0,1.0,1.0,1.0,1.0);
  }

};

struct KKpipi
{
  double mpi;
  double pipx1;
  double pipy1;
  double pipz1;
  double pipx2;
  double pipy2;
  double pipz2;
  double mk;
  double kpx1;
  double kpy1;
  double kpz1;
  double kpx2;
  double kpy2;
  double kpz2;

public:
  KKpipi(){
    mpi = 0.13957;
	mk = 0.493677;
  }
  void SetVal(double px1,double py1,double pz1,
              double px2,double py2,double pz2,
              double px3,double py3,double pz3,
              double px4,double py4,double pz4 )
  {
    pipx1 = px1;
    pipy1 = py1;
    pipz1 = pz1;
    pipx2 = px2;
    pipy2 = py2;
    pipz2 = pz2;
    kpx1  = px3;
    kpy1  = py3;
    kpz1  = pz3;
    kpx2  = px4;
    kpy2  = py4;
    kpz2  = pz4;
  }
  double InvMass(double f1=1.0,double f2=1.0,double f3=1.00,double f4=1.00)
  {
    double totpx,totpy,totpz,tote;
	double ke[2],pie[2];
	double mass;
    totpx = (f1*pipx1+f2*pipx2)+(f3*kpx1+f4*kpx2);
    totpy = (f1*pipy1+f2*pipy2)+(f3*kpy1+f4*kpy2);
    totpz = (f1*pipz1+f2*pipz2)+(f3*kpz1+f4*kpz2);
    ke[0]=CalEne(mk,f3*kpx1,f3*kpy1,f3*kpz1);
    ke[1]=CalEne(mk,f4*kpx2,f4*kpy2,f4*kpz2);
    //double massjpsi = CalInvMass(ml,lpx1,lpy1,lpz1,ml,lpx2,lpy2,lpz2);
    pie[0]=CalEne(mpi,f1*pipx1,f1*pipy1,f1*pipz1);
    pie[1]=CalEne(mpi,f2*pipx2,f2*pipy2,f2*pipz2);
    tote=ke[0]+ke[1]+pie[0]+pie[1];
    mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	return mass;
  } 
  double GetP1pi()
  {
    return CalMom(pipx1,pipy1,pipz1);
  }
  double GetP2pi()
  {
    return CalMom(pipx2,pipy2,pipz2);
  }
  double GetP1k()
  {
    return CalMom(kpx1,kpy1,kpz1);
  }
  double GetP2k()
  {
    return CalMom(kpx2,kpy2,kpz2);
  }

};

struct NPi
{
  int num;
  double mpi;
  std::vector<double> pipx;
  std::vector<double> pipy;
  std::vector<double> pipz;

public:
  NPi(){
    num = 0;
    mpi = 0.13957;
  }
  void SetVal( int n, double* px1,double* py1,double* pz1)
  {
    num = n;
	for (int i=0;i<num;i++){
      pipx.push_back(px1[i]);
      pipy.push_back(py1[i]);
      pipz.push_back(pz1[i]);
	}

  }
  void AddParticle(double px, double py, double pz)
  {
    num ++;
	pipx.push_back(px);
	pipy.push_back(py);
	pipz.push_back(pz);
  }
  bool CheckVolume()
  {
    if (num==pipx.size()) return true;
	return false;
  }
  int Clear(){
    pipx.clear();
    pipy.clear();
    pipz.clear();
	num = pipx.size();
	return num;
  }
  double InvMass(double *f = 0)
  {
    double facs[num];
	double mass;
    double totpx=0,totpy=0,totpz=0,tote=0;
	for (int i=0; i<num;i++){
      f==0 ? facs[i] = 1.0: facs[i] = f[i];
      totpx += facs[i]*pipx.at(i);
      totpy += facs[i]*pipy.at(i);
      totpz += facs[i]*pipz.at(i);
      tote  += CalEne(mpi,facs[i]*pipx.at(i),facs[i]*pipy.at(i),facs[i]*pipz.at(i));
	}
    mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	return mass;
  } 
  double GetP(int idx)
  {
    return CalMom(pipx.at(idx),pipy.at(idx),pipz.at(idx));
  }

};

namespace PeakEstimate{
//private:
  // par[0]: p distribution slope
  // par[1]: gaus mean;
  // par[2]: gaus sigma
  // par[3]: dm(p) E1
  // E1 = par[3];
  // E2 = par[4];
  // p1 = par[5];
  // p2 = par[6];
  // costheta = par[7];
  //double pdfpar[3];
  //double pspar[2];
  //double dmpar[5];

//public:
  //void SetPdfPar(double a, double b, double c)
  //{}
  //void SetPsPar(double mean, double sigma)
  //{
  //  pspar[0] = mean;
	//pspar[1] = sigma;
  //}
  double mresonace = 1.86484;
  double mfinal    = 0.493677;
  double mk = mfinal;
  double pars[8] ;
  double &p1dismean = pars[0] = 1.5;
  double &p1dissigma= pars[1] = 0.15;
  double &p2dismean = pars[2] = 1.5;
  double &p2dissigma= pars[3] = 0.15;
 
  void SetPdisPar(double m1, double s1, double m2, double s2)
  {
    p1dismean  = m1;
	p1dissigma = s1;
	p2dismean  = m2;
	p2dissigma = s2;
  }

  double GetResolution(double p)
  {
    return 0.01*p;
  }
  
  double Gauss2D(double *x, double *par);
  double ppdf(double *x, double *par);
  double psmear(double *x, double *par);
  double deltam(double *x, double *par);
  double multifgdelta(double *x, double *par);
  double multifg(double *x, double *par);
  double dmatp(double *x, double *par);
  double multifdm(double *x, double *par);
  double peakshift(double p1low=0, double p1up=0, double p2low=0, double p2up=0,int idx=0);

};

  double PeakEstimate::Gauss2D(double *x, double *par)
  {
    if (par[2]<=0 || par[4]<=0) return 0;
	double rx = (x[0]-par[1])/par[2];
	double ry = (x[1]-par[3])/par[4];
	return par[0]*TMath::Exp(-(rx*rx+ry*ry)/2.0);
  }
  double PeakEstimate::ppdf(double *x, double *par)
  {
    //return par[0]*x[0];
	if (pars[1]<0 || pars[3]<0) return 0;
	  
	double rx = (x[0]-pars[0])/pars[1];
	double ry = (x[1]-pars[2])/pars[3];

	//double res = TMath::Gaus(x[0],par[0],par[1],true);
	double res = TMath::Exp(-(rx*rx+ry*ry)/2.);
    //std::cout<<"in ppdf, p1 is "<<x[0]<<", mean is "<<par[0]<<", sigma is "<< par[1] <<", value is "<<res<< std::endl;
	return res;
  }
  double PeakEstimate::psmear(double *x, double *par)
  { 
    double sigma = GetResolution(x[0]);
    double res = TMath::Gaus(par[0],x[0],sigma,true);
    //std::cout<<"in psmear, p1 "<<par[0]<<", mean is "<<x[0]<<", sigma is "<< sigma <<", value is "<<res << std::endl;
    return res;
  }
  double PeakEstimate::deltam(double *x, double *par)
  {
    //std::cout<<"in deltam"<<std::endl;
    double m,E1,E2,p1,p2,costheta;
    m = mresonace;
	double mk = mfinal;
    //std::cout<<"in deltam  aaa"<<std::endl;
	p1 = par[0];
    p2 = par[1];
    //std::cout<<"in deltam bbb"<<std::endl;
	double p1r = x[0]; // p1 real position, smear to par[0]
	double p2r = x[1];
	E1 = sqrt(mk*mk + p1r*p1r);
	E2 = sqrt(mk*mk + p2r*p2r);
    costheta = ((E1+E2)*(E1+E2)-m*m-(p1r*p1r + p2r*p2r))/(2*p1r*p2r);
    //std::cout<<"E1: "<<E1<<"\tE2: "<<E2 <<"\tp1: "<<p1<<"\tp2: "<<p2<<"\tcostheta: "<<costheta<< std::endl;
    //std::cout<<"in multifgdelta return"<<std::endl;
	if (fabs(costheta)>1) return -1;

    double delta;
    double dp1 = p1 - x[0];
	double dp2 = p2 - x[1];
    delta = 1/m*((E1+E2)*(p1r*dp1/E1+p2r*dp2/E2)-(p1r*dp1+p2r*dp2+p2r*dp1*costheta+p1r*dp2*costheta));
    //std::cout<<"in delta m, p1 is "<<p1<<", dp1 is "<<dp1<<", delta m is "<< delta << std::endl;

	return delta;
  }
  double PeakEstimate::multifgdelta(double *x, double *par)
  {
    //std::cout<<"in multifgdelta"<<std::endl;
	//std::cout<<"fgd "<<ppdf(x,0)<<" "<<psmear(&x[0],&par[0])*psmear(&x[1],&par[1])<<" "<<deltam(x,par)<<std::endl;
    if (deltam(x,par) <-0.9) return 0;
    //std::cout<<"in multifgdelta"<<std::endl;
    return ppdf(x,0)*psmear(&x[0],&par[0])*psmear(&x[1],&par[1])*deltam(x,par);
  }
  double PeakEstimate::multifg(double *x, double *par)
  {
    double res = ppdf(x,0)*psmear(x,&par[0])*psmear(&x[1],&par[1]);
	//if (res < 1e-20)
	//std::cout<<"in fg ppdf is "<< ppdf(x,0)<<"x is ("<<x[0]<<","<<x[1]<< ") psmear is "<<psmear(x,&par[0])<<" "<<psmear(&x[1],&par[1])<<std::endl;
	return res;
  }
  double PeakEstimate::dmatp(double *x, double *par)
  {
    //std::cout<<"in dmatp"<<std::endl;
    double dm;

	double pars1[2];
	pars1[0] = x[0];// p1 smear to x[0] prob
	pars1[1] = x[1];// p2 smear to x[1]

    double p1 = x[0];
	double p2 = x[1];
    double s1 = GetResolution(x[0]);
	double s2 = GetResolution(x[1]);
    TF2 fgd("fgd",multifgdelta,0,3,0,3,2);
    TF2 fg("fg",multifg,0,3,0,3,2);
    fgd.SetParameters(p1,p2);
	fg.SetParameters(p1,p2);
	//std::cout<<"p1,p2:"<<p1<<" "<<p2<<std::endl;
    
	if ( fgd(x[0],x[1])<-0.9 ) return 0;
	if ( fgd(x[0]-5*s1,x[1]-5*s1)<-0.9 ) return 0;
	if ( fgd(x[0]+5*s2,x[1]+5*s2)<-0.9 ) return 0;

	//std::cout<<"bbb"<<std::endl;
    double integral1 = fgd.Integral(p1-5*s1,p1+5*s1,p2-5*s2,p2+5*s2);
	//std::cout<<"aaaaaa fgd  integral is "<< integral1<< std::endl;
	double integral2 = fg.Integral(p1-5*s1,p1+5*s1,p2-5*s2,p2+5*s2);
	//std::cout<<"aaaaaa fg  integral is "<< integral2<< std::endl;
	//std::cout<<"aaaaaa fgd integral is "<< integral1<<", fg int "<<integral2<<" p1 p2 "<<p1<<" "<<p2<< std::endl;
	if (fabs(integral2-0)<1e-100) return 0;
	return integral1/integral2;
  }
  double PeakEstimate::multifdm(double *x, double *par)
  {
    //std::cout<<"in multifdm "<<ppdf(x,0)<<" "<<dmatp(x,0)<<std::endl;
    return ppdf(x,0)*dmatp(x,0);
  }

  double PeakEstimate::peakshift(double p1low,double p1up,double p2low,double p2up,int idx)
  {
 
    TF2 fdm("fdm",   multifdm,0,2,0);
    TF2 pdis("pdis", ppdf,    0,2,0);
    std::cout<<"peak shift " <<std::endl;
	double fdmint = fdm.Integral(p1low,p1up,p2low,p2up);
	double pdisint= pdis.Integral(p1low,p1up,p2low,p2up);
    double ps = fdmint/pdisint;
	//std::cout<<"int1 int2 : "<<fdmint<<" " <<pdisint<<std::endl;
    std::cout<<"peak shift is "<< ps <<std::endl;
    if (idx ==0 ) return ps+mresonace;
	return ps;
  }



#endif
