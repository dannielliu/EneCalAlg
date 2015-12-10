#include "TMatrixD.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TAxis.h"
#include <iostream>
using namespace std;

int SolveLinear()
{

  // R value
  double a[5] = {0.742425, 0.389758, 0.305, /*0.305040,*/ 0, 0};
  double b[5] = {0.848947, 1.055110, 1.022850, 0, 0};
  double c[5] = {-0.0008515, -0.000880, -0.001050, 1.000902, 1.00071 /*0.994773*/};
  double pk[5] = {0.928916, 0.52321, 0.524507, 0, 0.851903 /*0.492212*/};
  double ppi[5] = {0.892368, 0.474128, 0.404327, 0.250183, 0};
  double num[5] = {22458, 22142, 65432, 8570, 60412};

  // uncertainty
  double e_a[5] = {4.654e-3, 4.061e-3, 2.398e-3, 0, 0};
  double e_b[5] = {4.424e-3, 4.25e-3,  2.451e-3, 0, 0};
  double e_c[5] = {5.8e-5, 5.3e-5, 3.1e-5, 7.7e-5, 8.7e-5};
  double e_pka[5] = {0.2056, 0.231066, 0.197752, 0, 0.217029};
  double e_ppi[5] = {0.171203, 0.193948, 0.167805, 0.078192,0};

  // 4230
//double a[4] = {0.740950, 0.392329, 0.306608, /*0.305040,*/ 0};
//double b[4] = {0.858381, 1.051814, 1.02682, 0};
//double c[4] = {-0.000914, -0.000785, -0.000997, 1.000902};
//double pk[4] = {0.941364, 0.55598, 0.543015, 0};
//double ppi[4] = {0.899468, 0.470103, 0.403764, 0.250245};
  
  // 4260
//double a[4] = {0.741456, 0.382856, 0.304960, 0};
//double b[4] = {0.849607, 1.055529, 1.019443, 0};
//double c[4] = {-0.000827, -0.000739, -0.000895, 1.00079};
//double pk[4] = {0.952286, 0.562599, 0.550774, 0};
//double ppi[4] = {0.90419, 0.472172, 0.403715, 0.250578};

  // 4360
//double a[4] = {0.738415, 0.390205, 0.306430, 0};
//double b[4] = {0.846814, 1.063350, 1.022483, 0};
//double c[4] = {-0.000858, -0.000777, -0.000941, 1.00077};
//double pk[4] = {0.957453, 0.552201, 0.555413, 0};
//double ppi[4] = {0.910932, 0.474918, 0.408327, 0.251365};

  // 4420
//double a[4] = {0.739081, 0.389891, 0.305391, 0};
//double b[4] = {0.847329, 1.054495, 1.020434, 0};
//double c[4] = {-0.001132, -0.001109, -0.001216, 1.00125};
//double pk[4] = {0.965522, 0.552528, 0.555959, 0};
//double ppi[4] = {0.908237, 0.47599, 0.409435, 0.251763};

  // 4600
 // double a[4] = {0.738413, 0.383921, 0.301900, 0};
 // double b[4] = {0.846875, 1.057320, 1.024816, 0};
 // double c[4] = {-0.001067, -0.001152, -0.001182, 1.00132};
 // double pk[4] = {0.984861, 0.579023, 0.577792, 0};
 // double ppi[4] = {0.918968, 0.480159, 0.411049, 0.252376};

  
  
  double aa[5][5];
  double xx[5];
  double yy[5];
  for (int i=0;i<3;i++){
    yy[i] = a[i] + b[i] - c[i];
    aa[i][0] = a[i];
    aa[i][1] = a[i]*pk[i];
    aa[i][2] = b[i];
    aa[i][3] = b[i]*ppi[i];
    aa[i][4] = b[i]*pow(ppi[i],2);
  }
  yy[3] = c[4];
  aa[3][0] = 1.0;
  aa[3][1] = pk[4];
  aa[3][2] = 0;
  aa[3][3] = 0;
  aa[3][4] = 0;


  yy[4] = c[3];//1.000902;
  aa[4][0] = 0;
  aa[4][1] = 0;
  aa[4][2] = 1;
  aa[4][3] = ppi[3];
  aa[4][4] = pow(ppi[3],2);

  cout << "########################################"<<endl;
  cout << "# D->Kpi D->K2pi D->KK psip->pipiJ/psi #"<<endl;
  cout << "########################################"<<endl;
  TMatrixD M_a(5,5);
  TVectorD v_x(5);
  TVectorD v_y(5);
  M_a.SetMatrixArray(*aa);
  v_x.SetElements(xx);
  v_y.SetElements(yy);

  cout<< "Print M_a" << endl;
  M_a.Print();
  M_a.Invert();
  cout<< "Print M_a after invert" << endl;
  M_a.Print();
  cout<< "Print v_y" << endl;
  v_y.Print();
  v_x = M_a*v_y;
  cout<< "aa "<< xx[0] << endl;
  v_x.Print();
  
  TF1* fpi = new TF1("fpi","[0]+[1]*x+[2]*x*x",0,2);
  fpi->SetParameter(0,v_x[2]);
  fpi->SetParameter(1,v_x[3]);
  fpi->SetParameter(2,v_x[4]);
  fpi->GetXaxis()->SetTitle("p_{#pi} (GeV/c)");
  fpi->GetYaxis()->SetTitle("f_{#pi}");
  fpi->Draw();

  //return 0;

  cout << "########################"<<endl;
  cout << "#D->Kpi D->K2pi D->K3pi#"<<endl;
  cout << "########################"<<endl;
  // fk as a constant
  double aaa[3][3];
  double xxx[3];
  double yyy[3];
  for (int i=0;i<3;i++){
    yyy[i] = a[i]+b[i] - c[i];
    aaa[i][0] = a[i];
    aaa[i][1] = b[i];
    aaa[i][2] = b[i]*ppi[i];
  }
  TMatrixD M_a3(3,3);
  TVectorD v_x3(3);
  TVectorD v_y3(3);
  M_a3.SetMatrixArray(*aaa);
  v_y3.SetElements(yyy);

  cout<<"Print Matrix 3x3:"<<endl;
  M_a3.Print();
  M_a3.Invert();
  cout<<"after invert"<<endl;
  M_a3.Print();
  cout<<"Print v_y3"<<endl;
  v_y3.Print();
  v_x3 = M_a3*v_y3;
  cout<<"Print xxx"<<endl;
  v_x3.Print();

  cout << "################################"<<endl;
  cout << "#D->Kpi D->K2pi psip->pipiJ/psi#"<<endl;
  cout << "################################"<<endl;
  // another solution 3
  aaa[2][0] = 0; aaa[2][1] = 1; aaa[2][2] = ppi[3];
  yyy[2] = c[3];
  M_a3.SetMatrixArray(*aaa);
  v_y3.SetElements(yyy);

  cout<<"Print Matrix 3x3:"<<endl;
  M_a3.Print();
  M_a3.Invert();
  cout<<"after invert"<<endl;
  M_a3.Print();
  cout<<"Print v_y3"<<endl;
  v_y3.Print();
  v_x3 = M_a3*v_y3;
  cout<<"Print xxx"<<endl;
  v_x3.Print();
  
  double a2[2][2];
  double y2[2];
  double* pka = &pk[0];
  double e_a2[2][2];
  double e_y2[2];
  
  TMatrixD M_a2(2,2);
  TVectorD v_x2(2);
  TVectorD v_y2(2);
  TMatrixD dM_a2(2,2);
  TVectorD dv_x2(2);
  TVectorD dv_y2(2);
  double pkaave;
  double ppiave;
  
  cout << "#####################"<<endl;
  cout << "#  D->Kpi D->K2pi  ##"<<endl;
  cout << "#####################"<<endl;
    a2[0][0] = a[0];                e_a2[0][0] = e_a[0];
    a2[0][1] = b[0];                e_a2[0][1] = e_b[0];
    a2[1][0] = a[1];                e_a2[1][0] = e_a[1];
    a2[1][1] = b[1];                e_a2[1][1] = e_b[1];
    y2[0]    =/* a[0] + b[0] */- c[0];  e_y2[0] = sqrt(/*pow(e_a[0],2) + pow(e_b[0],2) +*/ pow(e_c[0],2));
    y2[1]    =/* a[1] + b[1] */- c[1];  e_y2[1] = sqrt(/*pow(e_a[1],2) + pow(e_b[1],2) +*/ pow(e_c[1],2));
  M_a2.SetMatrixArray(*a2);
  v_y2.SetElements(y2);
  dM_a2.SetMatrixArray(*e_a2);
  dv_y2.SetElements(e_y2);
  pkaave = (pka[0]*num[0]+pka[1]*num[1])/(num[0]+num[1]);
  ppiave = (ppi[0]*num[0]+ppi[1]*num[1])/(num[0]+num[1]);
  cout << "average pk = "<<pkaave << "\t ppi = "<<ppiave << endl; 

  // dx = inv(A)*( db -dA*x)
  cout<<"Print Matrix 2x2:"<<endl;
  M_a2.Print();
  M_a2.Invert();
 // cout<<"after invert"<<endl;
 // M_a2.Print();
  cout<<"Print v_y2"<<endl;
  v_y2.Print();
  v_x2 = M_a2*v_y2;
  cout<<"Print xx"<<endl;
  v_x2.Print();
  
  cout<<"Print error Matrix 2x2:"<<endl;
  dM_a2.Print();
  cout<<"Print dv_y2"<<endl;
  dv_y2.Print();
  dv_x2 = M_a2*(dv_y2-dM_a2*v_x2);
  cout<<"Print dv_x2"<<endl;
  dv_x2.Print();
  return 0;

  cout << "#####################"<<endl;
  cout << "#  D->Kpi D->K3pi  ##"<<endl;
  cout << "#####################"<<endl;
    a2[0][0] = a[0];
    a2[0][1] = b[0];
    y2[0]    = a[0] + b[0] - c[0];
    a2[1][0] = a[2];
    a2[1][1] = b[2];
    y2[1]    = a[2] + b[2] - c[2];
  M_a2.SetMatrixArray(*a2);
  v_y2.SetElements(y2);
  pkaave = (pka[0]*num[0]+pka[2]*num[2])/(num[0]+num[2]);
  ppiave = (ppi[0]*num[0]+ppi[2]*num[2])/(num[0]+num[2]);
  cout << "average pk = "<<pkaave << "\t ppi = "<<ppiave << endl; 

  cout<<"Print Matrix 2x2:"<<endl;
  M_a2.Print();
  M_a2.Invert();
  cout<<"after invert"<<endl;
  M_a2.Print();
  cout<<"Print v_y2"<<endl;
  v_y2.Print();
  v_x2 = M_a2*v_y2;
  cout<<"Print xx"<<endl;
  v_x2.Print();

  cout << "#####################"<<endl;
  cout << "#  D->K2pi D->K3pi  ##"<<endl;
  cout << "#####################"<<endl;
    a2[0][0] = a[1];
    a2[0][1] = b[1];
    y2[0]    = a[1] + b[1] - c[1];
    a2[1][0] = a[2];
    a2[1][1] = b[2];
    y2[1]    = a[2] + b[2] - c[2];
  M_a2.SetMatrixArray(*a2);
  v_y2.SetElements(y2);
  pkaave = (pka[2]*num[2]+pka[1]*num[1])/(num[2]+num[1]);
  ppiave = (ppi[2]*num[2]+ppi[1]*num[1])/(num[2]+num[1]);
  cout << "average pk = "<<pkaave << "\t ppi = "<<ppiave << endl; 

  cout<<"Print Matrix 2x2:"<<endl;
  M_a2.Print();
  M_a2.Invert();
  cout<<"after invert"<<endl;
  M_a2.Print();
  cout<<"Print v_y2"<<endl;
  v_y2.Print();
  v_x2 = M_a2*v_y2;
  cout<<"Print xx"<<endl;
  v_x2.Print();
 
  return 0;
}
