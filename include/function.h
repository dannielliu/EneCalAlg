#ifndef myfunction_h
#define myfunction_h
#include "TMath.h"

double BreitWigner(double *x, double *par);

double line(double *x, double *par);

double line2(double *x, double *par);

double parabola(double *x, double *par);

double fitfun(double *x, double *par);

double GausLineBack(double *x,double *par);

// gaussian function, ignore the background
double maxlikelihood(double *x,double *par);

// function with background
double maxlikelihood1(double *x,double *par);
double maxlikelihood1_1(double *x,double *par);
double maxlikelihood1_3(double *x,double *par);

// for two diffrent particles end
double maxlikelihood2(double *x,double *par);//pure gaussian
double maxlikelihood2_0(double *x,double *par);//2D, background variable
double maxlikelihood2_1(double *x,double *par);//1D, background fixed
double maxlikelihood2_2(double *x,double *par);//2D, background fixed

// for two pairs of particles end
double maxlikelihood4(double *x,double *par);
double maxlikelihood4_0(double *x,double *par);//2D function
double maxlikelihood4_1(double *x,double *par);//1D

// for three particles
double maxlikelihood3(double *x,double *par);
double maxlikelihood3_0(double *x,double *par);//two dimension
double maxlikelihood3_1(double *x,double *par);//one dimension

double maxlikelihood0(double *x, double *par);

double CalInvMass(double m1, double px1, double py1, double pz1,
                  double m2, double px2, double py2, double pz2,
                  int n=0, const double *x=0, const double *par=0);




#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
class MultiDimFunction: public ROOT::Math::IParametricGradFunctionMultiDim
{
private:
  MultiDimFunction(){}
  ~MultiDimFunction(){}
public:
  MultiDimFunction(int nx, double start=0, double stop=2.0)
  {
    npar=2;
    double *par;
    par = new double[npar];
    SetParameters(par);
    par[0] = start;
    par[1] = stop;
    ndim = nx;
    fNpx=30;
    fXmin=0.99;
    fXmax=1.01;
  }
  double GetMinimum(double *xmin, double *xmax);
  void GetMinimumX(double *x);
  void SetRange(double xmin,double xmax){ fXmin=xmin; fXmax=xmax;}
  void ShowParameters();

private:
  const double *pars;// 2 dimension defined in construct function
  int ndim;
  int npar;
  double pstart;
  double pstop;
  int fNpx;
  double fXmin;
  double fXmax;
  double fValMin;

public:
  double DoEvalPar(const double *x, const double *p) const;
  unsigned int NDim() const;
  ROOT::Math::IParametricGradFunctionMultiDim* Clone() const
  { return new MultiDimFunction();}
  const double* Parameters() const { return pars;}
  void SetParameters(const double *p) { pars=p;}
  unsigned int NPar() const;
  double DoParameterDerivative(const double *x, const double *p,unsigned ipar) const;

};


double BiasCoe(const double *x,double *par);

double distribution(double *x, double *p);

#endif
