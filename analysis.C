#include <iostream>
#include "gepep_fastpipill.h"
#include "gepep_fast4pi.h"
#include "gepep_fast6pi.h"
#include "gepep_4k.h"
#include "gepep_kk.h"
#include "gepep_kpi.h"
#include "gepep_kpi2.h"
#include "gepep_fkkpipi.h"
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
	std::cout<<"Please specify an input root file ! "<<std::endl;
    return 0;
	//filename = "data_Rvalue_fpipill_e3850.root";
  }

  bool usefpipill=false;
  bool usef4pi=false;
  bool usef6pi=false;
  bool use4k=false;
  bool usekk=false;
  bool usekpi=false;
  bool usekpi2=false;
  bool usefkkpipi=false;
  usefpipill = ( filename.find("fpipill")   != std::string::npos
              || filename.find("fastpipill")!= std::string::npos);
  usef4pi    = ( filename.find("f4pi")      != std::string::npos
              || filename.find("fast4pi")   != std::string::npos);
  usef6pi    = ( filename.find("f6pi")      != std::string::npos
              || filename.find("fast6pi")   != std::string::npos);
  use4k      = ( filename.find("4k")        != std::string::npos);
  usekk      = ( filename.find("_kk_")      != std::string::npos);
  usekpi     = ( filename.find("_kpi_")     != std::string::npos);
  usekpi2    = ( filename.find("_kpi2_")    != std::string::npos);
  usefkkpipi = ( filename.find("fkkpipi")   != std::string::npos
              || filename.find("fastkkpipi")!= std::string::npos);

  if( usefpipill ){
    gepep_fastpipill *a;
	analysis(a,filename, "gepep_fastpipill");
  }
  else if( usef4pi ){
    gepep_fast4pi *a;
	analysis(a,filename, "gepep_fast4pi");
  }
  else if( usef6pi){
    gepep_fast6pi *a;
	analysis(a,filename, "gepep_fast6pi");
  }
  else if( use4k){
    gepep_4k *a;
	analysis(a,filename, "gepep_4k");
  }
  else if( usekk){
    gepep_kk *a;
	analysis(a,filename, "gepep_kk");
  }
  else if( usekpi){
    gepep_kpi *a;
	analysis(a,filename, "gepep_kpi");
  }
  else if( usekpi2){
    gepep_kpi2 *a;
	analysis(a,filename, "gepep_kpi");
  }
  else if( usefkkpipi){
    gepep_fkkpipi *a;
	analysis(a,filename, "gepep_fastkkpipi");
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
  TFile *file=new TFile(filename.c_str());
  if(!file){
    std::cout<<"can not open file "<<filename<<std::endl;
	return false;
  }
  TTree *tree = new TTree();;
  file->GetObject(treename.c_str(),tree);
  if(!tree){
    std::cout<<"can not find tree "<<treename<<std::endl;
	return false;
  }
  //~~~~~~~got :-)~~~~~~~~

  T events(tree);
  events.Loop();
  //delete f;
  delete tree;
  std::cout<<"after delete"<<std::endl;
  return true;
}

