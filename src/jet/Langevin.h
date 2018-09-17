#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>
#include "lorentz.h"
#include "random.h"
#include "JetEnergyLossModule.h"
#include "JetScapeConstants.h"

using namespace Jetscape; 

class Langevin : public JetEnergyLossModule<Langevin> //, public std::enable_shared_from_this<Langevin>
{
 private: 
	double alpha_s;
  	double g;
  	double pcut;        // below this scale, no further Eloss
  	double M; 
  	double mu; 
  	const double nc = 3.; 
  	const double CR = 3.; 

 public:  
  	Langevin();
  	virtual ~Langevin();

  	//main//
  	void Init();
  	void DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
 
	double qperp(double E, double T);
	double qpara(double E, double T);
	void Ito_update( double dt, double T, std::vector<double> v, const fourvec & pIn, fourvec & pOut);
}; 
#endif

