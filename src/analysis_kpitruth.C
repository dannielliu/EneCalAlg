#include <iostream>
#include "gepep_kpi.h"
#include "mctruth.h"
#include <string>
using namespace std;
string outputdir = ".";
mctruth *truthalg;

int main(int argc, char** argv)
{
  if (argc == 1 ) {
    cout<<"Please input a data file"<<endl;
    return -1;
  }
  TFile *file = new TFile(argv[1]);
  if (!file) return -2;
  TTree *tree = (TTree*)file->Get("gepep_kpi");
  TTree *treetr = (TTree*)file->Get("mctruth");
  if (!tree) return -3;
  if (!treetr) return -4;
  gepep_kpi obj(tree);
  truthalg = new mctruth(treetr);
  obj.Loop();

  cout<<"aaaaaaa"<<endl;
  delete truthalg;
  cout<<"aaaaaaa"<<endl;
  truthalg=NULL;
  //delete tree;
  cout<<"aaaaaaa"<<endl;
  //file->Close();
  cout<<"aaaaaaa"<<endl;
  //delete tree;
  //tree=NULL;

  return 0;
}

