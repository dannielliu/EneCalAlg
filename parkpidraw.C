#include <fstream>
#include <vector>
void parkpidraw()
{
  cout<<"start"<<endl;
  TH1D *h1=new TH1D("h1","factor K",50,0.996,1.004);
  //TH1D *h2=new TH1D("h2","factor mu",30,0.95,1.05);
  //TH1D *h3=new TH1D("h3","factor pi(ee)",30,0.95,1.05);
  //TH1D *h4=new TH1D("h4","factor pi(mumu)",30,0.95,1.05);
  ifstream in("parkpi");
  double factork;
  double factorerr;
  //double factormu;
  //double factorpi;
  //double factorelow;
  //double factoreup;
  //double factormuerr;
  //double factorpierr;
  double energy=0;
  double run;
  vector<double> x,y1,xe,y1e;

  x.clear();
  y1.clear();
  //y2.clear();
  //y3.clear();
  xe.clear();
  y1e.clear();
  //y3e.clear();
  char tmpstr[1000];
  getline(in,tmpstr);
  while(!in.eof()){
    in>>energy;
    in>>factork>>factorerr;
    //getline(in,tmpstr);
    //if(fabs(factorpi-1)>0.05){
    //  cout<<"May not fitted at energy: "<<energy<<endl;
    //}
    //if(fabs(factoreerr)>1){ 
    //  cout<<"error is too large at energy:  "<<energy<<endl;
    //  continue;
    //}
    if(fabs(factork )>5){
      cout<<"factor is too large/small at energy:  "<<energy<<endl;
      continue;
    }
    h1->Fill(factork);
    //h2->Fill(factormu);
    //h3->Fill(factorpi);
    //if(fabs(factorpierr)>1) continue;
    x.push_back(energy);
    xe.push_back(0.0);
    y1.push_back(factork);
    y1e.push_back(factorerr);
    //y2.push_back(factormu);
    //y2e.push_back(factormuerr);
    //y3.push_back(factorpi);
    //y3e.push_back(factorpierr);
    getline(in,tmpstr);
  }
  h1->Fit("gaus");
  //h2->Fit("gaus");
  //h3->Fit("gaus");

  const double n=x.size()-1;
  //const double n=104;
  double ene[n],enee[n];
  double fe[n],fee[n];
  //double fmu[n],fmue[n];
  //double fpi[n],fpie[n];
  for(int i=0;i<n;i++){
    ene[i]=x.at(i);
	enee[i]=0;
	fe[i]=y1.at(i);
	fee[i]=y1e.at(i);
	//feeup[i]=y1eup.at(i)-fe[i];
	//fmu[i]=y2.at(i);
	//fmue[i]=y2e.at(i);
	//fpi[i]=y3.at(i);
	//fpie[i]=y3e.at(i);
  }
  
  gStyle->SetOptFit(1111);
  //TF1 *f=new TF1("f","[0]",3800,4600);
  //f->SetParameter(0,1.0);
  cout<<"create graph,n point "<<n<<endl;
  TGraphErrors *graph1=new TGraphErrors(n,ene,fe,enee,fee);
  graph1->SetMarkerStyle(5);
  //graph1->Fit(f);
  graph1->SetMinimum(0.995);
  graph1->SetMaximum(1.005);
  //TGraphErrors *graph2=new TGraphErrors(n,ene,fmu,enee,fmue);
  //graph2->SetMarkerStyle(5);
  //graph2->Fit(f);
  //graph2->SetMinimum(0.95);
  //graph2->SetMaximum(1.05);
  //TGraphErrors *graph3=new TGraphErrors(n,ene,fpi,enee,fpie);
  //graph3->SetMarkerStyle(5);
  //graph3->SetMinimum(0.95);
  //graph3->SetMaximum(1.05);
  //graph3->Fit(f);
  TCanvas *c1=new TCanvas();
  //c1->Divide(2,2);
  //c1->cd(1);
  graph1->Draw("AP");
  graph1->SetTitle("D^{0} #rightarrow #pi K");
  c1->Print("factorsK.eps");
  //c1->cd(2);
  //graph2->Draw("AP");
  //c1->Print("factorsmu.eps");
  //c1->cd(3);
  //graph3->Draw("AP");
  //c1->Print("factorspi.eps");
  //exit();
}
