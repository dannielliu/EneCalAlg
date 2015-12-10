//#include <TAxis.h>
//#include <TF1.h>
//#include <TMath.h>
//#include <TCanvas.h>
//#include <iostream>
//#include <iomanip>
//using namespace std;

//int main()
{
  // D0 and D+
  double a1 = 0.742502	;	/*	0.740901;	*/	double a2 = 0.384245;
  double b1 = 0.849109	;	/*	0.847603;	*/	double b2 = 1.058387;
  double c1 = -.000812	;	/*	-0.0008415;	*/	double c2 = -0.000912;

  double a1e = 	0.004775;	/*	0.004073;	*/	double a2e = 0.005133;
  double b1e = 	0.004775;	/*	0.004073;	*/	double b2e = 0.005007;
  double c1e = 	0.000061;	/*	-0.000052;	*/	double c2e = -0.000067;
  
  // D0bar and D-
//double a1 = 0.741767	;		double a2 = 0.395000;
//double b1 = 0.847200	;		double b2 = 1.053754;
//double c1 = -.000775	;		double c2 = -0.000863;

//double a1e = 	0.003585;		double a2e = 0.003457;
//double b1e = 	0.003973;		double b2e = 0.003042;
//double c1e = 	0.000051;		double c2e = -0.000044;

// D0D0bar and D+D-
// 4230
//double a1 = 0.739945	;		double a2 = 0.390051;
//double b1 = 0.847342	;		double b2 = 1.052556;
//double c1 = -0.000921	;		double c2 = -0.000781;

//double a1e = 	0.003574;		double a2e = 0.003294;
//double b1e = 	0.003574;		double b2e = 0.002890;
//double c1e = 	0.000042;		double c2e = 0.000042;
// D0D0bar and D+D-
//double a1 = 0.740879	;		double a2 = 0.390981;
//double b1 = 0.846639	;		double b2 = 1.053499;
//double c1 = -.000790	;		double c2 = -0.000893;

//double a1e = 	0.003240;		double a2e = 0.002609;
//double b1e = 	0.003236;		double b2e = 0.002769;
//double c1e = 	0.000042;		double c2e = -0.000033;

// at 4575 MeV
// D0D0bar and D+D-
//double a1 = 0.729156	;		double a2 = 0.370591;
//double b1 = 0.842974	;		double b2 = 1.044851;
//double c1 = -.000480	;		double c2 = -0.000973;

//double a1e = 	0.019158;		double a2e = 0.018644;
//double b1e = 	0.019293;		double b2e = 0.024609;
//double c1e = 	0.000247;		double c2e = -0.000210;


std::cout<<a1*b2-a2*b1<<std::endl;
//dE = 0.847586*(x-1) + -0.000842      dE = 0.384245*(x-1) + -0.000912
//dE = 0.740922*(x-1) + -0.000841   dE = 1.058387*(x-1) + -0.000912
//dE = 0.740901*(x-1) + -0.000841
//dE = 0.847603*(x-1) + -0.000842
  
  double fpi =  (a1*c2-a2*c1)/(a2*b1-a1*b2);
  double fk  =  (b1*c2-b2*c1)/(b2*a1-b1*a2);
  double fpie = sqrt( 
    			TMath::Power(a2*fk*a1e/(a2*b1-a1*b2),2)
		+	TMath::Power(a1*fk*a2e/(a2*b1-a1*b2),2)
		+	TMath::Power(a2*fpi*b1e/(a2*b1-a1*b2),2)
		+	TMath::Power(a1*fpi*b2e/(a2*b1-a1*b2),2)
		+	TMath::Power(a2*c1e/(a2*b1-a1*b2),2)
		+	TMath::Power(a1*c2e/(a2*b1-a1*b2),2)
			);
  double fke  = sqrt( 
    			TMath::Power(a1*fpi*a2e/(a2*b1-a1*b2),2)
		+	TMath::Power(a2*fpi*a1e/(a2*b1-a1*b2),2)
		+	TMath::Power(a1*fk*b2e/(a2*b1-a1*b2),2)
		+	TMath::Power(a2*fk*b1e/(a2*b1-a1*b2),2)
		+	TMath::Power(a1*c2e/(a2*b1-a1*b2),2)
		+	TMath::Power(a2*c1e/(a2*b1-a1*b2),2)
			);
  std::cout<<setiosflags(ios::fixed)<<setprecision(6)
  	<<"f pi is "<< 1+fpi <<" +/- " << fpie
	<<", f k is "<< 1+fk <<" +/- " << fke <<std::endl;

  std::cout
    <<	a2*fk*a1e/(a2*b1-a1*b2) << " "
    <<	a1*fk*a2e/(a2*b1-a1*b2) << " "
    <<	a2*fpi*b1e/(a2*b1-a1*b2)<< " "
    <<  a1*fpi*b2e/(a2*b1-a1*b2)<< " "
    <<	a2*c1e/(a2*b1-a1*b2)    << " "
    <<  a1*c2e/(a2*b1-a1*b2)    << std::endl;

  std::cout
    <<	a1*fpi*a2e/(a2*b1-a1*b2)<< " "
    <<	a2*fpi*a1e/(a2*b1-a1*b2)<< " "
    <<	a1*fk*b2e/(a2*b1-a1*b2) << " "
    <<  a2*fk*b1e/(a2*b1-a1*b2) << " "
    <<	a1*c2e/(a2*b1-a1*b2)    << " "
    <<  a2*c1e/(a2*b1-a1*b2)    <<std::endl;
	
  TF1 *f1 = new TF1("fk_kpi","[0]*(x-1)+[1]+1",0.995,1.005);
  f1->SetParameters(-b1/a1, -c1/a1);
  f1->SetLineColor(2);
  f1->GetXaxis()->SetTitle("f_{#pi}");
  f1->GetXaxis()->SetTitleFont(62);
  f1->GetYaxis()->SetTitle("f_{K}");
  f1->GetYaxis()->SetTitleFont(62);

  TF1 *f2 = new TF1("fk_kpipi","[0]*(x-1)+[1]+1",0.99,1.01);
  f2->SetParameters(-b2/a2, -c2/a2);
  f2->SetLineColor(3);
  double deltafk1 = sqrt(
			TMath::Power(c1e/a1,2)
		+	TMath::Power(a1e*(c1)/(a1*a1),2)
  );
  double deltafk2 = sqrt(
			TMath::Power(c2e/a2,2)
		+	TMath::Power(a2e*(c2)/(a2*a2),2)
  );
  std::cout<<"fk1 error "<< deltafk1<<std::endl;
  std::cout<<"fk2 error "<< deltafk2<<std::endl;
  TF1 *f1_l = new TF1("fk_kpi_l","[0]*(x-1)+[1]+1-[2]",0.99,1.01);
  f1_l->SetParameters(-b1/a1, -c1/a1, deltafk1);
  f1_l->SetLineColor(2);
  TF1 *f1_u = new TF1("fk_kpi_u","[0]*(x-1)+[1]+1+[2]",0.99,1.01);
  f1_u->SetParameters(-b1/a1, -c1/a1, deltafk1);
  f1_u->SetLineColor(2);
  TF1 *f2_l = new TF1("fk_kpipi_l","[0]*(x-1)+[1]+1-[2]",0.99,1.01);
  f2_l->SetParameters(-b2/a2, -c2/a2,deltafk2);
  f2_l->SetLineColor(3);
  TF1 *f2_u = new TF1("fk_kpipi_u","[0]*(x-1)+[1]+1+[2]",0.99,1.01);
  f2_u->SetParameters(-b2/a2, -c2/a2,deltafk2);
  f2_u->SetLineColor(3);
 
  
  TCanvas *canvas1 = new TCanvas("canvas1","c1");
  f1->Draw();
  f1_l->Draw("same");
  f1_u->Draw("same");
  f2->Draw("same");
  f2_l->Draw("same");
  f2_u->Draw("same");

  // modify draw style
  int n = canvas1->GetListOfPrimitives()->GetSize();
  std::string histname;
  for (int i=0; i<n; i++){
 	if (!strncmp(canvas1->GetListOfPrimitives()->At(i)->ClassName(), "TF1",3)){
		histname = (canvas1->GetListOfPrimitives()->At(i))->GetName();
		break;
	}
  } 
  canvas1->SetMargin(0.15,0.15,0.15,0.15);

  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->SetTitle("Get Factor");
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetXaxis()->SetRangeUser(0.995,1.005);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetXaxis()->SetNdivisions(505);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetXaxis()->SetLabelSize(0.045);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTickLength(0.045);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitle("f_{#pi}");
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitleSize(0.045);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitleOffset(1.3);

  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetYaxis()->SetRangeUser(0.995,1.005);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetYaxis()->SetLabelSize(0.045);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTickLength(0.045);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitle("f_{K}");
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitleSize(0.045);
  ((TF1*)canvas1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitleOffset(1.3);

  TEllipse *circle = new TEllipse(fpi+1,fk+1,0.0005);
  circle->SetLineWidth(2);
  circle->SetFillStyle(0);
  circle->Draw();
  // add arrows point to intersection
//TArrow *arr1 = new TArrow(fpi+1,0.996,fpi+1,fk+1);
//TArrow *arr2 = new TArrow(0.996,fk+1,fpi+1,fk+1);
//arr1->SetArrowSize(0.035);
//arr2->SetArrowSize(0.035);
//arr1->SetLineWidth(2);
//arr2->SetLineWidth(2);
//arr1->Draw();
//arr2->Draw();
  return 0;
}
