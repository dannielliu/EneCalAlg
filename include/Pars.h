#include <map>
#include <fstream>

struct pars
{
  double par[4];
};

class ParMap
{
public:
  ParMap();
  ParMap(std::string parfile);
  ~ParMap();
  void init();
  double GetPar(double ene);
  double GetPar(int runNo);
  double GetEnergy(int runNo);

private:
  std::map<int,pars> fParMap;
  std::string fParFileName;
  bool InitStat;
};

