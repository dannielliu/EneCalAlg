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


struct Event{
  double px1;
  double py1;
  double pz1;
  double px2;
  double py2;
  double pz2;

  void SetVal(double x1, double y1, double z1, double x2, double y2, double z2){
    px1=x1;
    py1=y1;
    pz1=z1;
    px2=x2;
    py2=y2;
    pz2=z2;
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
  double InvMass(double f1=1.0, double f2=1.0)
  {
    double totpx,totpy,totpz,tote;
	double lee[2],pie[2];
	double mass;
    totpx = (f1*pipx1+f2*pipx2)+(lpx1+lpx2);
    totpy = (f1*pipy1+f2*pipy2)+(lpy1+lpy2);
    totpz = (f1*pipz1+f2*pipz2)+(lpz1+lpz2);
    lee[0]=CalEne(ml,lpx1,lpy1,lpz1);
    lee[1]=CalEne(ml,lpx2,lpy2,lpz2);
    //double massjpsi = CalInvMass(ml,lpx1,lpy1,lpz1,ml,lpx2,lpy2,lpz2);
    pie[0]=CalEne(mpi,f1*pipx1,f1*pipy1,f1*pipz1);
    pie[1]=CalEne(mpi,f2*pipx2,f2*pipy2,f2*pipz2);
    tote=lee[0]+lee[1]+pie[0]+pie[1];
    mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
    mass = mass-PsiM()+3.096916;
	return mass;
  }
  double PsiM()
  {
    return CalInvMass(ml,lpx1,lpy1,lpz1,ml,lpx2,lpy2,lpz2);
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


#endif
