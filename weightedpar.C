#include <fstream>
#include <vector>
void weightedpar()
{
  cout<<"start"<<endl;
  //TH1D *h1=new TH1D("h1","factor e",50,0.98,1.02);
  //TH1D *h2=new TH1D("h2","factor mu",30,0.95,1.05);
  //TH1D *h3=new TH1D("h3","factor pi(ee)",30,0.95,1.05);
  //TH1D *h4=new TH1D("h4","factor pi(mumu)",30,0.95,1.05);
  TH1D *h1=new TH1D("h","correction factor for pi",30,0.95,1.05);
  ifstream in1("parpipill");
  ifstream in2("parkpi");
  ifstream in3("parkpipi");
  double factore;
  double factormu;
  double factorpi1;
  double factorpi2;
  double factorpi3;
  double factoreerr;
  double factormuerr;
  double factorpi1err;
  double factorpi2err;
  double factorpi3err;
  double energy=0;
  vector<double> x,y1,y2,y3,xe,y1e,y2e,y3e;//y1elow,y1eup;

  x.clear();
  xe.clear();
  y1.clear();
  y2.clear();
  y3.clear();
  //y1elow.clear();
  //y1eup.clear();
  y1e.clear();
  y2e.clear();
  y3e.clear();
  char tmpstr[1000];
  getline(in1,tmpstr);
  while(!in1.eof()){
    in1>>energy;
    //energy+=1;
    //getline(in,tmpstr);
    //getline(in,tmpstr);
    in1>>factore>>factoreerr;
    in1>>factormu>>factormuerr;
    in1>>factorpi1>>factorpi1err;
    //getline(in,tmpstr);
    getline(in1,tmpstr);
    //h1->Fill(factore);
    //h2->Fill(factormu);
    //h3->Fill(factorpi);
    //if(fabs(factorpierr)>1) continue;
    x.push_back(energy);
    xe.push_back(0.0);
    y1.push_back(factorpi1);
    y1e.push_back(factorpi1err);
    //y1elow.push_back(factorelow);
    //y1eup.push_back(factoreup);
    //y2.push_back(factormu);
    //y2e.push_back(factormuerr);
    //y3.push_back(factorpi);
    //y3e.push_back(factorpierr);
  }
  getline(in2,tmpstr);
  while(!in2.eof()){
    in2>>energy;
    //energy+=1;
    in2>>factorpi2>>factorpi2err;
    getline(in2,tmpstr);
    //h1->Fill(factore);
    //h2->Fill(factormu);
    //h3->Fill(factorpi);
    //if(fabs(factorpierr)>1) continue;
    //x.push_back(energy);
    //xe.push_back(0.0);
    y2.push_back(factorpi2);
    y2e.push_back(factorpi2err);
    //y1elow.push_back(factorelow);
    //y1eup.push_back(factoreup);
    //y2.push_back(factormu);
    //y2e.push_back(factormuerr);
    //y3.push_back(factorpi);
    //y3e.push_back(factorpierr);
  }
  getline(in3,tmpstr);
  while(!in3.eof()){
    in3>>energy;
    //energy+=1;
    in3>>factorpi3>>factorpi3err;
    //getline(in,tmpstr);
    getline(in3,tmpstr);
    //h1->Fill(factore);
    //h2->Fill(factormu);
    //h3->Fill(factorpi);
    //if(fabs(factorpierr)>1) continue;
    //x.push_back(energy);
    //xe.push_back(0.0);
    y3.push_back(factorpi3);
    y3e.push_back(factorpi3err);
    //y1elow.push_back(factorelow);
    //y1eup.push_back(factoreup);
    //y2.push_back(factormu);
    //y2e.push_back(factormuerr);
    //y3.push_back(factorpi);
    //y3e.push_back(factorpierr);
  }

  const double n=x.size()-1;
  //const double n=104;
  double ene[n],enee[n];
  double fpi1[n],fpi1e[n],fpi2[n],fpi2e[n],fpi3[n],fpi3e[n];
  for(int i=0;i<n;i++){
    ene[i]=x.at(i);
    enee[i]=0;
    fpi1[i] =y1.at(i);
    fpi1e[i]=y1e.at(i);
    fpi2[i] =y2.at(i);
    fpi2e[i]=y2e.at(i);
    fpi3[i] =y3.at(i);
    fpi3e[i]=y3e.at(i);

  }

  ofstream of("pion.par");
  double factor;
  for(int i=0;i<n;i++){
    factor=(fpi1[i]/fpi1e[i]+fpi2[i]/fpi2e[i]+fpi3[i]/fpi3e[i])
           /(1./fpi1e[i]+1./fpi2e[i]+1./fpi3e[i]);
    of<<ene[i]<<"\t"<<factor<<"\n";
  }
  
   exit();
}
