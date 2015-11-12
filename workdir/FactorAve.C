#include <vector>
#include <utility>
#include <iostream>
#include "TMath.h"
int AveragePar(std::vector<std::pair<double, double> > factor)
{
    //std::vector<double> facs;
	//facs.push_back(1.000902);
    ;

	//factor.push_back(std::make_pair(1.000902, 0.000077));
	//factor.push_back(std::make_pair(1.000714, 0.000086));
	//factor.push_back(std::make_pair(1.000540, 0.000033));
	//factor.push_back(std::make_pair(1.000630, 0.000047));
	//factor.push_back(std::make_pair(4.571834, 0.004976));
	//factor.push_back(std::make_pair(4.570046, 0.000744));
	double facsum=0, errsum=0;
	for (int i=0; i<factor.size(); i++){
		facsum += factor.at(i).first/TMath::Power(factor.at(i).second,2);
		errsum += 1./TMath::Power(factor.at(i).second,2);
	}
	double ave_fac = facsum/errsum;
	double ave_err = 1./sqrt(errsum);
	std::cout << "Average factor is " << ave_fac <<", average error is "<< ave_err<<std::endl;
	return 0;
}

int main()
{
	std::vector<std::pair<double, double> > factor;
	
	std::cout<<"with correction factor : 1.0000"<<std::endl;
	factor.clear();
	factor.push_back(std::make_pair(4.571834, 0.004976));
	factor.push_back(std::make_pair(4.570046, 0.000744));
	AveragePar(factor);
	
	std::cout<<"with correction factor : 1.000815"<<std::endl;
	factor.clear();
	factor.push_back(std::make_pair(4.577786, 0.002674));
	factor.push_back(std::make_pair(4.573563, 0.000719));
	AveragePar(factor);
	
	std::cout<<"with correction factor : 1.000610"<<std::endl;
	factor.clear();
	factor.push_back(std::make_pair(4.576954, 0.002645));
	factor.push_back(std::make_pair(4.572780, 0.000683));
	AveragePar(factor);
}
