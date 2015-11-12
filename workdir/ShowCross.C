using namespace std;
int ShowCross()
{
  ifstream infile("CrossSection.txt");
  std::cout<<infile.is_open()<<std::endl;
  istringstream iss;
  double x[10][100];
  double xe[10][100];
  double y[10][100];
  double ye[10][100];
  int setid=-1;
  int setNo[10] = {0,0,0,0,0,0,0,0,0,0};
  char title[10][100];
  
  char line[1000];
  while (!infile.eof())
  {
    iss.clear();
    infile.getline(line,1000);
    string str = line;
    iss.str(str.c_str());
    std::cout<<str<<std::endl;
    if (str.length() == 0) continue;
    if (line[0] == ' ') continue;
    //std::cout<<iss.peek()<<std::endl;
    if (iss.peek()=='#')
    {
      setid ++;
      std::cout<<"new set "<<setid<<" "<<str<<std::endl;
      sprintf(title[setid],"%s", &line[1]);
      continue;
    }

    iss >> x[setid][setNo[setid]] >> y[setid][setNo[setid]];
    iss >> line ;
    iss >> ye[setid][setNo[setid]];
    xe[setid][setNo[setid]] = 0;
    setNo[setid] ++;
  }
  std::cout<<"setNo "<<setid<<std::endl;

  TCanvas *c1 = new TCanvas();
  TLegend *legend = new TLegend(0.7,0.7,0.89,0.89);
  legend->SetLineColor(0);
  for (int i=0; i<=setid; i++)
  {
    std::cout<<"set id "<< i <<" set size "<<setNo[i]<<std::endl;
    TGraphErrors *graph = new TGraphErrors(setNo[i],x[i],y[i],xe[i],ye[i]);
    graph->SetMarkerStyle(24+i);
    graph->SetMarkerColor(1+i);
    graph->SetFillColor(0);
    legend->AddEntry(graph,title[i]);
    if (i==0) {
      graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
      graph->GetXaxis()->SetLabelSize(0.04); 
      graph->GetXaxis()->SetTickLength(0.04);
      graph->GetXaxis()->SetTitleSize(0.04); 
      graph->GetXaxis()->SetTitleOffset(1.2);
      
      graph->GetYaxis()->SetTitle("Cross section (nb)");
      graph->GetYaxis()->SetLabelSize(0.04);
      graph->GetYaxis()->SetTickLength(0.04);
      graph->GetYaxis()->SetTitleSize(0.04);
      graph->GetYaxis()->SetTitleOffset(1.2);
      graph->GetYaxis()->SetRangeUser(0,0.3);
      
      graph->SetTitle("multiple hadrons cross section");
      graph->Draw("AP");
      continue;
    }
    graph->Draw("P");
  }
  legend->Draw();
}
