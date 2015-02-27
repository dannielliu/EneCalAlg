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
  double InvMass(double f1=1.0,double f2=1.0,double f3=1.001,double f4=1.001)
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

#endif
