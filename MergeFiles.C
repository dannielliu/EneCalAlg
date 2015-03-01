#include "TFile.h"
#include <iostream>
#include "TFileMerger.h"
#include "TTree.h"

int main(int argc,char** argv)
{
  TFileMerger merger;
  char infilename[1000][1000];
  TFile *file[1000];
  char outfilename[1000];
  int fileNo = argc-2;
  if (argc<3) {
    std::cout<<"Too few input arguements"<<std::endl;
	return -1;
  }
  sprintf(outfilename,"%s",argv[1]);
  merger.OutputFile(outfilename);
  for (int i=0;i<fileNo;i++){
    sprintf(infilename[i],"%s",argv[i+2]);
	file[i] = new TFile(infilename[i]);
	merger.AddFile(file[i]);
  }
  merger.Merge();
  merger.PrintFiles("");
  return 0;
}
