#include "function.h"
#include "TMath.h"
#include <fstream>

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

 

double GetEnergy(int runNo)
{
  ifstream inene("curene");
  if (inene.is_open()) {
    double tmpene=0;
    inene >> tmpene;
    if (fabs(tmpene-0)>1e-10) return tmpene;
  }
  //runNo=runNo%10000;
  //int a[2] = {1, 2};
  runNo = abs(runNo);
  if (runNo>=32239 && runNo <= 32864) return 4.23;
  else if(runNo>=29677 && runNo<=30367) return 4.26;
  else if(runNo>=31561 && runNo<=31981) return 4.26;
  else if(runNo>=30616 && runNo<=31279) return 4.36;
  else if(runNo>=31327 && runNo<=31390) return 4.42;
  else if(runNo>=36773 && runNo<=38140) return 4.42;
  else if(runNo>=36245 && runNo<=36393) return 4.47;
  else if(runNo>=36398 && runNo<=36588) return 4.53;
  else if(runNo>=36603 && runNo<=36699) return 4.575;
  else if(runNo>=35227 && runNo<=36213) return 4.60;

  else if(runNo>=33490 && runNo<=33556) return 3.81;
  else if(runNo>=33572 && runNo<=33657) return 3.90;
  else if(runNo>=33659 && runNo<=33719) return 4.09;
  else if(runNo>=30372 && runNo<=30437) return 4.19;
  else if(runNo>=31983 && runNo<=32045) return 4.21;
  else if(runNo>=32046 && runNo<=32140) return 4.22;
  else if(runNo>=30438 && runNo<=30491) return 4.23;
  else if(runNo>=32141 && runNo<=32226) return 4.245;
  else if(runNo>=30492 && runNo<=30557) return 4.31;
  else if(runNo>=31281 && runNo<=31325) return 4.39;
  int runNoLow[108] = {
  34011, 34028, 34037, 34046, 34058, 34069, 34077, 34084, 34091, 34097,
  34105, 34118, 34128, 34135, 34142, 34151, 34161, 34175, 34184, 34191,
  34197, 34203, 34211, 34221, 34231, 34240, 34246, 34253, 34258, 34266,
  34272, 34277, 34282, 34298, 34314, 34321, 34328, 34339, 34346, 34351,
  34359, 34369, 34374, 34382, 34390, 34397, 34404, 34412, 34418, 34428,
  34437, 34447, 34461, 34478, 34486, 34494, 34503, 34512, 34527, 34541,
  34555, 34564, 34574, 34585, 34593, 34603, 34613, 34623, 34634, 34642,
  34652, 34661, 34674, 34685, 34695, 34705, 34719, 34729, 34740, 34754,
  34763, 34777, 34785, 34794, 34804, 34812, 34825, 34837, 34848, 34861,
  34869, 34882, 34891, 34900, 34913, 34926, 34936, 34947, 34958, 34968,
  34982, 35010, 35027, 35041, 35060, 35082, 35099, 35119 };// #the last number is not used!
  double energy[108] = {
   3850,  3890,  3895,  3900,  3905,  3910,  3915,  3920,  3925,  3930,
   3935,  3940,  3945,  3950,  3955,  3960,  3965,  3970,  3975,  3980,
   3985,  3990,  3995,  4000,  4005,  4010,  4012,  4014,  4016,  4018,
   4020,  4025,  4030,  0000,  4035,  4040,  4050,  4055,  4060,  4065,
   4070,  4080,  4090,  4100,  4110,  4120,  4130,  4140,  4145,  4150,
   4160,  4170,  0000,  4180,  4190,  4195,  4200,  4203,  4206,  4210,
   4215,  4220,  4225,  4230,  4235,  4240,  4243,  4245,  4248,  4250,
   4255,  4260,  4265,  4270,  4275,  4280,  4285,  4290,  4300,  4310,
   4320,  4330,  4340,  4350,  4360,  4370,  4380,  4390,  4395,  4400,
   4410,  4420,  4425,  4430,  4440,  4450,  4460,  4480,  4500,  4520,
   4540,  4550,  4560,  4570,  4580,  0000,  4590,  4600 };// # the last one is skipped
  for(int i=0;i<107;i++){
    if(runNo>=runNoLow[i]&&runNo<runNoLow[i+1])
      return energy[i]/1000;
  }
  return -1;
}

  
  
  
  
  double GetWeight(double p, int i)
  {  
    if (i==2){
	  if (p<0.25)	return 2.27593*exp(1.78172-34.8690*p)+0.116624;
	  if (p>=0.25 && p<1.4)	return 9.83740e-2+0.161337*p-0.142985*p*p;
	  if (p>1.4)	return 0.299732*exp(1.37126+46.6799*(p-1.50106))+5.15151e-2;
	}
	if (i==1) return 1.0;
	
	return 0;
  }
  double GetResolution(double p, int i)
  { 
    if (i==1) {
      if (p<=0.2) 			return 3.91343e-3+exp(-2.26966-3.71626e1*p);
      if (p>0.2 && p<=0.6) 	return 4.98844e-3-7.35034e-3*p+9.21407e-3*p*p;
      if (p>0.6)  			return 1.9451e-3+3.00209e-3*p;
    }
    if (i==2) {
      if (p<0.8)			return 1.03546e-2+exp(-7.08236-9.50215*(p-0.393409));
      if (p>=0.8)			return 1.12534e-2+exp(-1.51850+13.3226*(p-1.56431));
    }
    if (i==0) {
      double res1 = GetResolution(p,1);
      double res2 = GetResolution(p,2);
      double w2 = GetWeight(p,2);
      return	sqrt((res1*res1+res2*res2*w2)/(1+w2));
    }
    return 1;
  }
  double IntegralTF2User(TF2 &f2, double *xmin, double *xmax, int Ndf=1000)
  {
    double integral_res = 0;
	double xgap = (xmax[0]-xmin[0])/Ndf;
	double ygap = (xmax[1]-xmin[1])/Ndf;
	double xsqu = xgap*ygap;
	for (int i=0; i<Ndf; i++){
	  for (int j=0; j<Ndf; j++){
	    double xpos = xmin[0]+i*xgap;
		double ypos = xmin[1]+j*ygap;
	    integral_res += xsqu*f2.Eval(xpos+xgap/2,ypos+ygap/2);
	  }
	}
    
	return integral_res;

  }
  double IntegralTF2User(TF2 &f2, double xmin, double xmax, double ymin, double ymax, int Ndf=1000)
  {
    double xlow[] = {xmin, ymin};
	double xup[]  = {xmax, ymax};
	return IntegralTF2User(f2,xlow,xup,Ndf);
  }



  void PeakEstimate::SetPdisPar(double m1, double s1, double m2, double s2)
  {
    PeakEstPar *par = PeakEstPar::GetInstance();
    par->p1dismean  = m1;
    par->p1dissigma = s1;
    par->p2dismean  = m2;
    par->p2dissigma = s2;
  }
  void PeakEstimate::SetMResonace(double m)
  {
    PeakEstPar *par = PeakEstPar::GetInstance();
    par->mresonace = m;
  }
  void PeakEstimate::SetMFinal(double m)
  {
    PeakEstPar *par = PeakEstPar::GetInstance();
    par->mfinal = m;
  }

//double GetResolution(double p)
//{
//  double res;
//  //if (p<=0) return 0;

//  //if (p<=0.15)          res = 2.88453e-02-2.47417e-01*p+6.14261e-01*p*p;
//  //if (p>0.15 && p<=0.3) res = 6.69506e-03-1.61844e-02*p+2.33373e-02*p*p;
//  //if (p>0.3  && p<=0.6) res = 4.82524e-03-4.99686e-03*p+7.00845e-03*p*p;
//  //if (p>0.6 )           res = 2.57764e-03+2.82971e-03*p;
//  
//  // mc resolution from J/psi --> pi+ pi- pi0 pi0
//  if (p<=0.2) res = TMath::Exp(-1.12842-42.2116*p)+4.32780e-03;
//  if (p>0.2 && p<=0.6) res = 5.40457e-03-7.57167e-03*p+9.72436e-03*p*p;
//  if (p>0.6) res = 2.55808e-03+2.83029e-03*p;
//  if (res>0.2) res=0.2;
//  
//  // real data, from ee->pi+ pi- K+ K-, missing pi+ get p of pi+ using p conservation, compare it with reconstruct from detector
//  //res = TMath::Exp(-3.63696-2.15276*p)+9.32902e-03;
//  return res*p;
//}

  double PeakEstimate::Gauss2D(double *x, double *par)
  {
    if (par[2]<=0 || par[4]<=0) return 0;
    double rx = (x[0]-par[1])/par[2];
    double ry = (x[1]-par[3])/par[4];
    return par[0]*TMath::Exp(-(rx*rx+ry*ry)/2.0);
  }
  double PeakEstimate::Maxwell(double *x, double *par)
  {
    if (x[0]>par[1] && par[2]>0 && par[0]>0)
	{
	  return par[0]*(x[0]-par[1])/par[2]*TMath::Exp(-(x[0]-par[1])/par[2]);
	}
	else if (x[0]<=par[1] && par[2]>0 && par[0]>0)
	  return par[0]*(x[0])/par[2]*TMath::Exp(-(x[0]-par[1])/par[2]);
	else return 0;
  }
  double PeakEstimate::ppdf(double *x, double *par)
  {
    //return par[0]*x[0];
    PeakEstPar *ins = PeakEstPar::GetInstance();
    double pars[4];
    pars[0] = ins->p1dismean;
    pars[1] = ins->p1dissigma;
    pars[2] = ins->p2dismean;
    pars[3] = ins->p2dissigma;
    if (pars[1]<0 || pars[3]<0) return 0;
    if (x[0]<0 || x[1]<0) return 0;
          
    //double rx = (x[0]-pars[0])/pars[1];
    if (x[0]<pars[0]) return 0;
    //if (x[1]<pars[2]) return 0;
    double rx = sqrt(x[0]-pars[0])*exp(-(x[0]-pars[0])/(pars[1]));
    double ry = (x[1]-pars[2])*exp(-(x[1]-pars[2])*(x[1]-pars[2])/pars[3]);
    //double ry = (x[1]-pars[2])/pars[3];
    //double ry = TMath::Gaus(x[1],pars[2],pars[3],true);
    //double rx = TMath::Gaus(x[0],pars[0],pars[1],true);
    //double ry = TMath::Gaus(x[1],pars[2],pars[3],true);

    //double res = TMath::Gaus(x[0],par[0],par[1],true);
    //double res = TMath::Exp(-(rx*rx+ry*ry)/2.);
    
    double res = rx*ry;
    //std::cout<<"in ppdf, p1 is "<<x[0]<<", mean is "<<pars[0]<<", sigma is "<< pars[1] <<", value is "<<rx<< std::endl;
    //std::cout<<"in ppdf, p2 is "<<x[1]<<", mean is "<<pars[2]<<", sigma is "<< pars[3] <<", value is "<<ry<< std::endl;
	return res;
  }
  double PeakEstimate::psmear(double *x, double *par)
  { 
    //double sigma = GetResolution(x[0]);
    double w2 = GetWeight(par[0],2);
    double s1 = GetResolution(par[0],1);
    double s2 = GetResolution(par[0],2);
    //std::cout<<w2<<' '<<s1<<' ' <<s2<<std::endl;
    double gaus1 = TMath::Gaus(x[0],par[0],s1*par[0],true);
    double gaus2 = TMath::Gaus(x[0],par[0],s2*par[0],true);
    //double xp = (x[0]-par[0])/par[0];
    //double gaus1 = 1/(sqrt(2*TMath::Pi())*s1*par[0])*TMath::Exp(-xp*xp/(2*s1*s1));
    //double gaus2 = 1/(sqrt(2*TMath::Pi())*s2*par[0])*TMath::Exp(-xp*xp/(2*s2*s2));
    return (gaus1+w2*gaus2)/(1+w2);
    //return (TMath::Gaus(xp,0,s1,true)+w2*TMath::Gaus(xp,0,s2,true))/(1+w2); // need to be normalized

	//double res = TMath::Gaus(par[0],x[0],sigma,true);
    //std::cout<<"in psmear, p1 "<<par[0]<<", mean is "<<x[0]<<", sigma is "<< sigma <<", value is "<<res << std::endl;
    //return res;
  }
  double PeakEstimate::deltam(double *x, double *par)
  {
    //std::cout<<"in deltam"<<std::endl;
    PeakEstPar *ins = PeakEstPar::GetInstance();
    double m,E1,E2,p1,p2,costheta;
    m = ins->mresonace;
    double mk = ins->mfinal;
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
    if (fabs(costheta)>1) return -1.e100;

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
    double fg  =  ppdf(x,0)*psmear(&x[0],&par[0])*psmear(&x[1],&par[1]);
	double res =  fg*deltam(x,par);
	
    if (deltam(x,par) <-1.e99) res = 0;
	//std::cout<<"fgd: x is ("<<x[0]<<","<<x[1]<< ") " <<ppdf(x,0)<<" "<<psmear(&x[0],&par[0])*psmear(&x[1],&par[1])<<" "<<deltam(x,par);
	//std::cout<<"\t fgd "<<res<<" , fg "<<fg<<std::endl;
    //std::cout<<"in multifgdelta"<<std::endl;
    return	res;//ppdf(x,0)*psmear(&x[0],&par[0])*psmear(&x[1],&par[1])*deltam(x,par);
  }
  double PeakEstimate::multifg(double *x, double *par)
  {
    double dm = deltam(x,par);
    if (dm < -1.e99) return 0;
	double pdf = ppdf(x,0);
	double psm = psmear(&x[0],&par[0])*psmear(&x[1],&par[1]);
    double res =  pdf*psm;
	//double fgd = res*deltam(x,par);
	//std::cout<<"int fg, "<<"x is ("<<x[0]<<", "<<x[1]<< ") fg is "<<res<<", fgd is "<<fgd<<std::endl;
	
	//if (res < 1e-20)
	//std::cout<<"in fg ppdf is "<<pdf<<" x is ("<<x[0]<<","<<x[1]<< ") psmear is "<<psm<<" deltam is "<< dm <<" fg is "<<res <<std::endl;
	return res;
  }

//#include "Math/WrappedMultiTF1.h"
//#include "Math/AdaptiveIntegratorMultiDim.h"
//#include "Math/IntegratorMultiDim.h"
////#include "TF2.h"
//#include "Math/AllIntegrationTypes.h"
//#include "Math/GSLMCIntegrator.h"
////#include "TMath.h"
  double PeakEstimate::dmatp(double *x, double *par)
  {
    //std::cout<<"in dmatp"<<std::endl;
    double dm;

////double pars1[2];
////pars1[0] = x[0];// p1 smear to x[0] prob
////pars1[1] = x[1];// p2 smear to x[1]

    double p1 = x[0];
	double p2 = x[1];
    double s1 = GetResolution(x[0],2);
	double s2 = GetResolution(x[1],2);
    TF2 fgd("fgd",multifgdelta,0,3,0,3,2);
    TF2 fg("fg",multifg,0,3,0,3,2);
    fgd.SetParameters(p1,p2);
	fg.SetParameters(p1,p2);

    int n = 5;
    double xmin[]={p1-n*s1,p2-n*s2};
    double xmax[]={p1+n*s1,p2+n*s2};
    double integral1 = IntegralTF2User(fgd,xmin,xmax,100);
    double integral2 = IntegralTF2User(fg,xmin,xmax,100);

////ROOT::Math::WrappedMultiTF1 wfgd(fgd);
////ROOT::Math::WrappedMultiTF1 wfg(fg);
////ROOT::Math::GSLMCIntegrator ig1(ROOT::Math::IntegrationMultiDim::kPLAIN);
////ROOT::Math::GSLMCIntegrator ig2(ROOT::Math::IntegrationMultiDim::kPLAIN);
////ig1.SetFunction(wfgd);
////ig2.SetFunction(wfg);
////ig1.SetRelTolerance(1e-1);
////ig2.SetRelTolerance(1e-1);
////int n = 5;
////double xmin[]={p1-n*s1,p2-n*s2};
////double xmax[]={p1+n*s1,p2+n*s2};
////double integral1 = ig1.Integral(xmin,xmax);
////double integral2 = ig2.Integral(xmin,xmax);
	//std::cout<<"p1,p2:"<<p1<<" "<<p2<<std::endl;
    
	//if ( fgd(x[0],x[1])<-0.9 ) return 0;
	//if ( fgd(x[0]-10*s1,x[1]-10*s1)<-0.9 ) return 0;
	//if ( fgd(x[0]+10*s2,x[1]+10*s2)<-0.9 ) return 0;
	
	//std::cout<<"bbb"<<std::endl;
    //double integral1 = fgd.Integral(p1-n*s1,p1+n*s1,p2-n*s2,p2+n*s2);
	//std::cout<<"aaaaaa fgd  integral is "<< integral1<< std::endl;
	//double integral2 = fg.Integral(p1-n*s1,p1+n*s1,p2-n*s2,p2+n*s2);
	//std::cout<<"aaaaaa fg  integral is "<< integral2<< std::endl;
	//std::cout<<"aaaaaa fgd integral is "<< integral1<<", fg int "<<integral2<<" p1 p2 "<<p1<<" "<<p2<<"\n"<< std::endl;
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
    //return PeakEstPar::GetInstance()->mresonace;
 
    TF2 fdm("fdm",   multifdm,0,2,0);
    TF2 pdis("pdis", ppdf,    0,2,0);
    std::cout<<"peak shift " <<std::endl;
////ROOT::Math::WrappedMultiTF1 wfdm(fdm);
////ROOT::Math::WrappedMultiTF1 wfdis(pdis);
////ROOT::Math::GSLMCIntegrator ig1(ROOT::Math::IntegrationMultiDim::kVEGAS);
////ROOT::Math::GSLMCIntegrator ig2(ROOT::Math::IntegrationMultiDim::kVEGAS);
////ig1.SetFunction(wfdm);
////ig2.SetFunction(wfdis);
////ig1.SetRelTolerance(1e-6);
////ig2.SetRelTolerance(1e-6);
////int n = 5;
////double xmin[]={p1low,p2low};
////double xmax[]={p1up, p2up};
////double fdmint = ig1.Integral(xmin,xmax);
////double pdisint = ig2.Integral(xmin,xmax);

    double fdmint = IntegralTF2User(fdm,p1low,p1up,p2low,p2up,100);
    double pdisint = IntegralTF2User(pdis,p1low,p1up,p2low,p2up,100);
    //double fdmint = fdm.Integral(p1low,p1up,p2low,p2up);
    //std::cout<<"peak shift ..." <<std::endl;
    //double pdisint= pdis.Integral(p1low,p1up,p2low,p2up);
    double ps=0;
    if (fabs(pdisint)>1e-100) ps = fdmint/pdisint;
    std::cout<<"int1 int2 : "<<fdmint<<" " <<pdisint<<std::endl;
    std::cout<<"peak shift is "<< ps <<std::endl;
    if (idx ==0 ) return ps+PeakEstPar::GetInstance()->mresonace;
    return ps;
  }
