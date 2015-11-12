#ifndef EventClass_h
#define EventClass_h
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
using namespace CLHEP;


class PSIP{
public:
  HepLorentzVector pip;
  HepLorentzVector pim;
  HepLorentzVector ele;
  HepLorentzVector pos;

public:
  PSIP()
  {}
  ~PSIP(){}

  void set(HepLorentzVector fp1, HepLorentzVector fp2, HepLorentzVector fp3, HepLorentzVector fp4)
  {
    pip = fp1;
    pim = fp2;
    ele = fp3;
    pos = fp4;
  }
  double m()
  {
    return (pip+pim+ele+pos).m() - M_jpsi() + 3.096916;
  }
  double M_jpsi()
  {
    return (ele+pos).m();
  }
  inline Hep3Vector GetPpip(){return pip.vect();}
  inline Hep3Vector GetPpim(){return pim.vect();}
  inline Hep3Vector GetPele(){return ele.vect();}
  inline Hep3Vector GetPpos(){return pos.vect();}

  PSIP &operator *= (const double f)
  {
    double mpi = pip.m();
    pip.setVectM(pip.vect()*f,mpi);
    pim.setVectM(pim.vect()*f,mpi);
    return *this;
  }
  void setCorrectionFactors(double f1, double f2)
  {
	double mpi = pip.m();
    pip.setVectM(pip.vect()*f1,mpi);
    pim.setVectM(pim.vect()*f2,mpi);
  }


};

class PHI{
public:
  HepLorentzVector kap;
  HepLorentzVector kam;

public:
  PHI() {}
  ~PHI(){}

  void set(HepLorentzVector fp1, HepLorentzVector fp2)
  {
    kap = fp1;
	kam = fp2;
  }
  double m()
  {
    return (kap+kam).m();
  }
  inline Hep3Vector GetPkap(){return kap.vect();}
  inline Hep3Vector GetPkam(){return kam.vect();}

  PHI &operator *= (const double f)
  {
    double mk = kap.m();
    kap.setVectM(kap.vect()*f,mk);
    kam.setVectM(kam.vect()*f,mk);
    return *this;
  }

};

class D0KPI{
public:
  HepLorentzVector kam;
  HepLorentzVector pip;
public:
  D0KPI() {}
  ~D0KPI() {}
  
  void set(HepLorentzVector fp1, HepLorentzVector fp2)
  {
    kam = fp1;
    pip = fp2;
  }
  double m()
  {
    return (kam+pip).m();
  }
  inline Hep3Vector GetPkam(){return kam.vect();}
  inline Hep3Vector GetPpip(){return pip.vect();}

  D0KPI &operator *= (const double f)
  {
    double mk = kam.m();
    double mpi = pip.m();
    kam.setVectM(kam.vect()*f,mk);
    pip.setVectM(pip.vect()*f,mpi);
    return *this;
  }
  void setCorrectionFactors(double f)
  {
    double mk = kam.m();
    double mpi = pip.m();
    kam.setVectM(kam.vect()*f,mk);
    pip.setVectM(pip.vect()*f,mpi);
  }
  void setCorrectionFactors(double f1, double f2)
  {
    double mk = kam.m();
    double mpi = pip.m();
    kam.setVectM(kam.vect()*f1,mk);
    pip.setVectM(pip.vect()*f2,mpi);
  }
  
  HepLorentzVector Get4P()
  {
    return kam+pip;
  }
};

class D0K3PI{
public:
  HepLorentzVector kam;
  HepLorentzVector pip1;
  HepLorentzVector pip2;
  HepLorentzVector pim;
public:
  D0K3PI() {}
  ~D0K3PI() {}
  
  void set(HepLorentzVector fkm1, HepLorentzVector fp1, HepLorentzVector fp2, HepLorentzVector fpm1)
  {
    kam  = fkm1;
    pip1 = fp1;
    pip2 = fp2;
    pim  = fpm1;
  }
  double m()
  {
    return (kam+pip1+pip2+pim).m();
  }
  inline Hep3Vector GetPkam(){return kam.vect();}
  inline Hep3Vector GetPpip1(){return pip1.vect();}
  inline Hep3Vector GetPpip2(){return pip2.vect();}
  inline Hep3Vector GetPpim(){return pim.vect();}

  D0K3PI &operator *= (const double f)
  {
    double mk = kam.m();
    double mpi = pip1.m();
    kam.setVectM(kam.vect()*f,mk);
    pip1.setVectM(pip1.vect()*f,mpi);
    pip2.setVectM(pip2.vect()*f,mpi);
    pim.setVectM(pim.vect()*f,mpi);
    return *this;
  }
  void setCorrectionFactors(double f)
  {
    double mk = kam.m();
    double mpi = pip1.m();
    kam.setVectM(kam.vect()*f,mk);
    pip1.setVectM(pip1.vect()*f,mpi);
    pip2.setVectM(pip2.vect()*f,mpi);
    pim.setVectM(pim.vect()*f,mpi);
  }
  void setCorrectionFactors(double f1, double f2)
  {
    double mk = kam.m();
    double mpi = pip1.m();
    kam.setVectM(kam.vect()*f1,mk);
    pip1.setVectM(pip1.vect()*f2,mpi);
    pip2.setVectM(pip2.vect()*f2,mpi);
    pim.setVectM(pim.vect()*f2,mpi);
  }
  void setCorrectionFactors(double f1, double f2, double f3, double f4)
  {
    double mk = kam.m();
    double mpi = pip1.m();
    kam.setVectM(kam.vect()*f1,mk);
    pip1.setVectM(pip1.vect()*f2,mpi);
    pip2.setVectM(pip2.vect()*f3,mpi);
    pim.setVectM(pim.vect()*f4,mpi);
  }
  
  HepLorentzVector Get4P()
  {
    return kam+pip1+pip2+pim;
  }
};
  
class D2KPiPi{
public:
  HepLorentzVector kam;
  HepLorentzVector pip;
  HepLorentzVector pim;
public:
  D2KPiPi() {}
  ~D2KPiPi() {}
  
  void set(HepLorentzVector fp1, HepLorentzVector fp2, HepLorentzVector fp3)
  {
    kam = fp1;
	pip = fp2;
	pim = fp3;
  }
  double m()
  {
    return (kam+pip+pim).m();
  }
  void setCorrectionFactors(double fk, double fpi)
  {
    double mk = kam.m();
	double mpi = pip.m();
    kam.setVectM(kam.vect()*fk,mk);
    pip.setVectM(pip.vect()*fpi,mpi);
    pim.setVectM(pim.vect()*fpi,mpi);
  
  }
  void setCorrectionFactors(double fk, double fpi1, double fpi2)
  {
    double mk = kam.m();
    double mpi = pip.m();
    kam.setVectM(kam.vect()*fk,mk);
    pip.setVectM(pip.vect()*fpi1,mpi);
    pim.setVectM(pim.vect()*fpi2,mpi);
  
  }
  D2KPiPi &operator *= (const double f)
  {
    double mk = kam.m();
	double mpi = pip.m();
    kam.setVectM(kam.vect()*f,mk);
    pip.setVectM(pip.vect()*f,mpi);
    pim.setVectM(pim.vect()*f,mpi);
    return *this;
  }

};

class PiPiPP{
public:
  HepLorentzVector pip;
  HepLorentzVector pim;
  HepLorentzVector pp;
  HepLorentzVector pm;
public:
  PiPiPP() {}
  ~PiPiPP() {}
  
  void set(	HepLorentzVector fp1, HepLorentzVector fp2, 
  			HepLorentzVector fp3, HepLorentzVector fp4)
  {
	pip = fp1;
	pim = fp2;
	pp  = fp3;
	pm  = fp4;
  }
  double m()
  {
    return (pip+pim+pp+pm).m();
  }
  void setCorrectionFactors(double fpi, double fp)
  {
    double mpi = pip.m();
    double mp  = pp.m();
    pip.setVectM(pip.vect()*fpi,mpi);
    pim.setVectM(pim.vect()*fpi,mpi);
    pp.setVectM(pp.vect()*fp,mp);
    pm.setVectM(pm.vect()*fp,mp);
  }
  void setCorrectionFactors(double fpi1, double fpi2, double fp1, double fp2)
  {
    double mpi = pip.m();
    double mp  = pp.m();
    pip.setVectM(pip.vect()*fpi1,mpi);
    pim.setVectM(pim.vect()*fpi2,mpi);
    pp.setVectM(pp.vect()*fp1,mp);
    pm.setVectM(pm.vect()*fp2,mp);
  }


  PiPiPP &operator *= (const double f)
  {
	double mpi = pip.m();
    double mp = pp.m();
    pip.setVectM(pip.vect()*f,mpi);
    pim.setVectM(pim.vect()*f,mpi);
    pp.setVectM(pp.vect()*f,mp);
    pm.setVectM(pm.vect()*f,mp);
    return *this;
  }

};

template <class T>
inline T operator * (const T &d0, double f)
{
  T d0tmp = d0;
  d0tmp *= f;
  return d0tmp;
}

class PiPiKK{
public:
  HepLorentzVector pip;
  HepLorentzVector pim;
  HepLorentzVector kap;
  HepLorentzVector kam;

public:
  PiPiKK() {}
  ~PiPiKK() {}

  void set( HepLorentzVector fpip, HepLorentzVector fpim, 
            HepLorentzVector fkap, HepLorentzVector fkam) {
	pip = fpip;
	pim = fpim;
	kap = fkap;
	kam = fkam;
  }
  double m()
  {
    return (pip+pim+kap+kam).m();
  }
  void setCorrectionFactors(double fpi, double fka)
  {
    double mpi = pip.m();
    double mka = kap.m();
    pip.setVectM(pip.vect()*fpi,mpi);
    pim.setVectM(pim.vect()*fpi,mpi);
    kap.setVectM(kap.vect()*fka,mka);
    kam.setVectM(kam.vect()*fka,mka);
  }
  void setCorrectionFactors(double fpi1, double fpi2, double fka1, double fka2)
  {
    double mpi = pip.m();
    double mka = kap.m();
    pip.setVectM(pip.vect()*fpi1,mpi);
    pim.setVectM(pim.vect()*fpi2,mpi);
    kap.setVectM(kap.vect()*fka1,mka);
    kam.setVectM(kam.vect()*fka2,mka);
  }

};

class FourPi{
public:
  HepLorentzVector pip1;
  HepLorentzVector pim1;
  HepLorentzVector pip2;
  HepLorentzVector pim2;

public:
  FourPi() {}
  ~FourPi() {}

  void set( HepLorentzVector fpip, HepLorentzVector fpim, 
            HepLorentzVector fkap, HepLorentzVector fkam) {
	pip1 = fpip;
	pim1 = fpim;
	pip2 = fkap;
	pim2 = fkam;
  }
  double m()
  {
    return (pip1+pim1+pip2+pim2).m();
  }
  void setCorrectionFactors(double fpi)
  {
    double mpi = pip1.m();
    pip1.setVectM(pip1.vect()*fpi,mpi);
    pim1.setVectM(pim1.vect()*fpi,mpi);
    pip2.setVectM(pip2.vect()*fpi,mpi);
    pim2.setVectM(pim2.vect()*fpi,mpi);
  }

};

class APair{
public:
  HepLorentzVector partilep;
  HepLorentzVector partilem;

public:
  APair() {}
  ~APair() {}

  void set( HepLorentzVector fpip, HepLorentzVector fpim) {
	partilep = fpip;
	partilem = fpim;
  }
  double m()
  {
    return (partilep+partilem).m();
  }
  void setCorrectionFactors(double f)
  {
    double mparticle = partilep.m();
    partilep.setVectM(partilep.vect()*f,mparticle);
    partilem.setVectM(partilem.vect()*f,mparticle);
  }
  void setCorrectionFactors(double f1, double f2)
  {
    double mparticle1 = partilep.m();
    double mparticle2 = partilem.m();
    partilep.setVectM(partilep.vect()*f1,mparticle1);
    partilem.setVectM(partilem.vect()*f2,mparticle2);
  }

}; 

class Lambdac{
public:
	HepLorentzVector proton;
	HepLorentzVector kam;
	HepLorentzVector pip;
public:
	Lambdac() {}
	~Lambdac() {}

	void set(HepLorentzVector fp, HepLorentzVector fkm, HepLorentzVector fpip){
		proton = fp;
		kam = fkm;
		pip = fpip;
	}
	double m(){
		return (proton+kam+pip).m();
	}
	void setCorrectionFactors(double f)
	{
		double mp = proton.m();
		double mk = kam.m();
		double mpi = pip.m();
		proton.setVectM(proton.vect()*f, mp);
		kam.setVectM(kam.vect()*f,mk);
		pip.setVectM(pip.vect()*f,mpi);
	}
	void setCorrectionFactors(double fp,double fk, double fpi)
	{
		double mp = proton.m();
		double mk = kam.m();
		double mpi = pip.m();
		proton.setVectM(proton.vect()*fp, mp);
		kam.setVectM(kam.vect()*fk,mk);
		pip.setVectM(pip.vect()*fpi,mpi);
	}

	HepLorentzVector get4P()
	{
	  return proton+kam+pip;
	}

};

#endif
