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

private:
  std::map<int,pars> fParMap;
  std::string fParFileName;
  bool InitStat;
};

ParMap::ParMap()
{
  fParFileName="";
  InitStat = false;
}

ParMap::ParMap(std::string parfile)
{
  fParFileName = parfile;
  init();
}

ParMap::~ParMap()
{
  fParMap.clear();
}

void ParMap::init()
{
  ifstream in;
  int ene,runNo;
  pars fact;
  in.open(fParFileName.c_str());
  while(!in.eof())
  {
    in>>ene>>runNo;
    in>>fact.par[0]>>fact.par[1]>>fact.par[2];
    fParMap.insert(std::make_pair(ene,fact));
  }
  InitStat = true;
}

double ParMap::GetPar(double ene)
{
  int ene2=(int)ene;
  if(InitStat) return fParMap[ene2].par[0];
  else return -1;
}

double ParMap::GetPar(int runNo)
{
  runNo=runNo%10000;
  //int a[2] = {1, 2};
  int runNoLow[108] = {
   4011, 4028, 4037, 4046, 4058, 4069, 4077, 4084, 4091, 4097,
   4105, 4118, 4128, 4135, 4142, 4151, 4161, 4175, 4184, 4191,
   4197, 4203, 4211, 4221, 4231, 4240, 4246, 4253, 4258, 4266,
   4272, 4277, 4282, 4298, 4314, 4321, 4328, 4339, 4346, 4351,
   4359, 4369, 4374, 4382, 4390, 4397, 4404, 4412, 4418, 4428,
   4437, 4447, 4461, 4478, 4486, 4494, 4503, 4512, 4527, 4541,
   4555, 4564, 4574, 4585, 4593, 4603, 4613, 4623, 4634, 4642,
   4652, 4661, 4674, 4685, 4695, 4705, 4719, 4729, 4740, 4754,
   4763, 4777, 4785, 4794, 4804, 4812, 4825, 4837, 4848, 4861,
   4869, 4882, 4891, 4900, 4913, 4926, 4936, 4947, 4958, 4968,
   4982, 5010, 5027, 5041, 5060, 5082, 5099, 5119 };// #the last number is not used!
  double energy[108] = {
   3850, 3890, 3895, 3900, 3905, 3910, 3915, 3920, 3925, 3930,
   3935, 3940, 3945, 3950, 3955, 3960, 3965, 3970, 3975, 3980,
   3985, 3990, 3995, 4000, 4005, 4010, 4012, 4014, 4016, 4018,
   4020, 4025, 4030, 0000, 4035, 4040, 4050, 4055, 4060, 4065,
   4070, 4080, 4090, 4100, 4110, 4120, 4130, 4140, 4145, 4150,
   4160, 4170, 0000, 4180, 4190, 4195, 4200, 4203, 4206, 4210,
   4215, 4220, 4225, 4230, 4235, 4240, 4243, 4245, 4248, 4250,
   4255, 4260, 4265, 4270, 4275, 4280, 4285, 4290, 4300, 4310,
   4320, 4330, 4340, 4350, 4360, 4370, 4380, 4390, 4395, 4400,
   4410, 4420, 4425, 4430, 4440, 4450, 4460, 4480, 4500, 4520,
   4540, 4550, 4560, 4570, 4580, 0000, 4590, 4600, };// # the last one is skipped
  for(int i=0;i<107;i++){
    if(runNo>=runNoLow[i]&&runNo<runNoLow[i+1])
      return GetPar(energy[i]);
  }
  return -2;
}
