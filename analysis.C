#include <iostream>
#include "gepep_fastpipill.h"
#include "gepep_fast4pi.h"
//#include "gepep_fast6pi.h"
#include "gepep_4k.h"
#include "gepep_kk.h"
#include "gepep_kpi.h"
#include "gepep_ppi.h"
#include "gepep_kpi2.h"
#include "gepep_kpipi.h"
#include "gepep_fkkpipi.h"
#include "gepep_pipipp.h"
//#include "KsAlg.h"
#include "mumu.h"
#include "Ks0Alg.h"
#include "TKey.h"
#include <string>

template<typename T> 
bool analysis(T *type, std::string &filename, std::string treename="");
std::string outputdir=".";

int main(int argc,char **argv)
{ 
  std::string filename;
  //std::cout<<"arg number is "<<argc<<std::endl;
//if(argc>=2){ 
//  filename=argv[1];
// std::cout<<"input file name is : "<<filename<<std::endl;
//}
//if(argc >=3){
//  outputdir=argv[2];
// std::cout<<"output dir is "<<outputdir<<std::endl;
//}
  if (argc==1){
    std::cout<<"Please specify an input root file ! "<<std::endl;
    return 0;
   //filename = "data_Rvalue_fpipill_e3850.root";
  }

  for (int i=1; i<argc; i++){
    filename = argv[i];
	if (filename.find(".root") ==std::string::npos) {
	  outputdir = filename;
	}
  }
  for (int i=1; i<argc; i++){
    filename = argv[i];
    if (filename.find(".root") == std::string::npos) continue;
    std::string namenopa = filename;
    while (namenopa.find("/") != std::string::npos) namenopa = &namenopa[namenopa.find("/")+1];
    std::cout<<"pure file name: "<< namenopa << std::endl;
    bool usefpipill=false;
    bool usemumu=false;
    bool usef4pi=false;
    //bool usef6pi=false;
    bool use4k=false;
    bool usekk=false;
    bool usekpi=false;
    bool useppi=false;
    bool usekpi2=false;
    bool usekpipi=false;
    bool usefkkpipi=false;
    bool usepipipp=false;
    bool useKs0Alg=false;
    usefpipill = ( namenopa.find("fpipill")   != std::string::npos
                || namenopa.find("fastpipill")!= std::string::npos);
    usemumu    = ( namenopa.find("mumu")      != std::string::npos);
    usef4pi    = ( namenopa.find("f4pi")      != std::string::npos
                || namenopa.find("fast4pi")   != std::string::npos);
    use4k      = ( namenopa.find("4k")        != std::string::npos);
    usekk      = ( namenopa.find("_kk_")      != std::string::npos
                || namenopa.find("kk_")       != std::string::npos
                || namenopa.find("KK.")       != std::string::npos);
    usekpi     = ( namenopa.find("_kpi_")     != std::string::npos
                || namenopa.find("Kpi.")     != std::string::npos);
    useppi     = ( namenopa.find("_ppi_")     != std::string::npos);
    usekpi2    = ( namenopa.find("_kpi2_")    != std::string::npos);
    usekpipi   = ( namenopa.find("_kpipi_")   != std::string::npos);
    usefkkpipi = ( namenopa.find("fkkpipi")   != std::string::npos
                || namenopa.find("fastkkpipi")!= std::string::npos);
    usepipipp = ( namenopa.find("pipipp")    != std::string::npos
                || namenopa.find("fpipipp")   != std::string::npos);
    useKs0Alg  = ( namenopa.find("Ksto2pi")   != std::string::npos
                || namenopa.find("Ks")        != std::string::npos);
 
 
    if( usefpipill ){
      gepep_fastpipill *a;
      analysis(a,filename);
    }
    else if( usemumu ){
      mumu *a;
      analysis(a,filename);
    }
    else if( usef4pi ){
      gepep_fast4pi *a;
      analysis(a,filename);
    }
    else if( use4k){
      gepep_4k *a;
      analysis(a,filename);
    }
    else if( usekk){
      gepep_kk *a;
      analysis(a,filename);
    }
    else if( usekpi){
      gepep_kpi *a;
      analysis(a,filename,"gepep_kpi");
    }
    else if( useppi){
      gepep_ppi *a;
      analysis(a,filename);
    }
    else if( usekpi2){
      gepep_kpi2 *a;
      analysis(a,filename);
    }
    else if( usekpipi){
      gepep_kpipi *a;
      analysis(a,filename,"gepep_kpi");
    }
    else if( usefkkpipi){
      gepep_fkkpipi *a;
      analysis(a,filename);
    }
    else if( usepipipp){
      gepep_pipipp *a;
      analysis(a,filename);
    }
    else if( useKs0Alg){
      Ks0Alg *a;
      analysis(a,filename, "Ks_info");
    }
    else{
      std::cout<<"Do not know how to deal with it!"<<std::endl;
      return 1;
    }
  
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

  // get tree
  if (treename == "") { 
  	int n = file->GetListOfKeys()->GetSize();
	for (int i=0; i<n; i++){
		if (!strncmp(((TKey*)file->GetListOfKeys()->At(i))->GetClassName(), "TTree",5)){
			treename = ((TKey*)file->GetListOfKeys()->At(i))->GetName();
			break;
		}
	}
  }
  std::cout<<"Tree name is "<<treename<<std::endl;
  
  TTree *tree = new TTree();
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
  tree = NULL; 
  std::cout<<"after delete"<<std::endl;
  return true;
}

