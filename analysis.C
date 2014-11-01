#include <iostream>
#include "gepep_fastpipill.h"
#include "gepep_fast4pi.h"
#include "gepep_fast6pi.h"
//#include "TFile.h"
#include <string>

template<typename T> 
bool analysis(T *type, std::string &filename, std::string treename);
std::string outputdir=".";

int main(int argc,char **argv)
{
  std::string filename;
  //std::cout<<"arg number is "<<argc<<std::endl;
  if(argc>=2){ 
    filename=argv[1];
	std::cout<<"input file name is : "<<filename<<std::endl;
  }
  if(argc >=3){
    outputdir=argv[2];
	std::cout<<"output dir is "<<outputdir<<std::endl;
  }
  if (argc==1){
	//std::cout<<"Please specify a input root file ! "<<std::endl;
    filename = "data_Rvalue_fpipill_e3850.root";
  }

  bool usepipill=false;
  bool use4pi=false;
  bool use6pi=false;
  usepipill = ( filename.find("fpipill")    != std::string::npos
             || filename.find("fastpipill") != std::string::npos);
  use4pi    = ( filename.find("f4pi")    != std::string::npos
             || filename.find("fast4pi") != std::string::npos);
  use6pi    = ( filename.find("f6pi")    != std::string::npos
             || filename.find("fast6pi") != std::string::npos);

  if( usepipill ){
    gepep_fastpipill *a;
	analysis(a,filename, "gepep_fastpipill");
  }
  else if( use4pi ){
    gepep_fast4pi *a;
	analysis(a,filename, "gepep_fast4pi");
  }
  else if( use6pi){
    gepep_fast6pi *a;
	analysis(a,filename, "gepep_fast6pi");
  }
  else{
    std::cout<<"Do not know how to deal with it!"<<std::endl;
	return 1;
  }

  std::cout<<"done"<<std::endl;
  return 0;
}

template<typename T> 
bool analysis(T *type,std::string &filename, std::string treename)
{
  TFile *f=new TFile(filename.c_str());
  if(!f){
    std::cout<<"can not open file "<<filename<<std::endl;
	return false;
  }
  TTree *tree = new TTree();;
  f->GetObject(treename.c_str(),tree);
  if(!tree){
    std::cout<<"can not find tree "<<treename<<std::endl;
	return false;
  }
  //~~~~~~~got :-)~~~~~~~~

  T events(tree);
  events.Loop();
  //delete f;
  delete tree;
  //std::cout<<"after delete"<<std::endl;
  return true;
}

