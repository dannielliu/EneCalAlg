#include <iostream>
#include "gepep_kpi.h"
#include <string>
using namespace std;
string outputdir = ".";

int main(int argc, char** argv)
{
  if (argc == 1 ) {
    cout<<"Please input a data file"<<endl;
    return -1;
  }
  TFile *file = new TFile(argv[1]);
  if (!file) return -2;
  TTree *tree = (TTree*)file->Get("gepep_kpi");
  if (!tree) return -3;
  gepep_kpi obj(tree);
  obj.Loop();

  delete tree;
  tree=NULL;

  return 0;
}

