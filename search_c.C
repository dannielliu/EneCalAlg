#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "sys/io.h"
#include <string>
#include "dir.h"
#include <fstream>

using namespace std;

int main(int argc,char** argv)
{
  char folder[500]={"."};
  if (argc>1) sprintf(folder,"%s",argv[1]);
  char filename[500];
  sprintf(filename,"%s/*.txt",folder);

  _finddata_t fileinfo;
  long handle;
  handle = _findfirst(filename,&fileinfo);
  if (handle==-1) {
    std::cout<<"fail to find file."<<std::endl;
    return -1;
  }
  while (_findnext(handle,&fileinfo)==0){
    cout<<"file: "<<fileinfo.name<<std::endl;
  }
  _findclose(handle);
  return 0;
}
