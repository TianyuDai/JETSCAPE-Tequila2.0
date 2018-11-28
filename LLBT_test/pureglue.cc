#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapePartonReader.h"
#include "JetScapeBanner.h"
#include "fjcore.hh"

#include <GTL/dfs.h>

using namespace std;
using namespace fjcore;

using namespace Jetscape;

// You could overload here and then simply use ofstream << p;
// ostream & operator<<(ostream & ostr, const fjcore::PseudoJet & jet);


// -------------------------------------

int main(int argc, char** argv)
{
  double scale_list[3] = {0.5, 1., 2.}; 
  double alpha_list[1] = {0.3}; 
  double omegacut_list[1] = {1.}; 

  for (int i = 0; i < 3; i++)
  {
  	  
	  double scale = scale_list[i]; 
	  for (int j = 0; j < 1; j++)
	  {
	  	double alpha_s = alpha_list[j]; 
		for (int k = 0; k < 1; k++)
		{
		  double omegacut = omegacut_list[k]; 
		  auto reader=make_shared<JetScapeReaderAscii>("./realhydro_20GeV_gluon_muscale"+std::to_string(scale)+"alpha"+std::to_string(alpha_s)+"omega"+std::to_string(omegacut)+"_medium_inel.dat"); 
		  std::ofstream jet_output (("../../Result/Tequila/elas/realhydro_20GeV_gluon_muscale"+std::to_string(scale)+"alpha"+std::to_string(alpha_s)+"omega"+std::to_string(omegacut)+"_medium_inel.txt").c_str());
		  const int nEvents = 10000; 
		  double Emean = 0.; 
		  double jetpTMin = 4., jetRadius = 0.4; 
		  int nParticles = 0; 
		  fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
		  vector <fjcore::PseudoJet> fjInputs; 
		  
		  vector <double> EBin(101); 
		  for (unsigned int j = 0; j < EBin.size(); j++)
		  	EBin[j] = 2.+(double) j * 20. / 100.; 
		  vector <double> cs(EBin.size()-1, 0.), err(EBin.size()-1, 0.); 
		  vector <int> sqSum(EBin.size()-1, 0); 
		  
		  while (!reader->Finished())
			{
			  reader->Next();
			  vector <int> ct(EBin.size()-1, 0); 
			  vector <shared_ptr <Parton> > partons; 
			  partons = reader->GetPartons(); 
			  nParticles += partons.size(); 
			  for (unsigned int i = 0; i < partons.size(); i++)
			  {
			  		double energy = partons[i]->e(); 
			  		Emean += energy; 
			  		for (unsigned int j = 0; j < EBin.size() - 1; j++)
			  			if (energy>=EBin[j] && energy<EBin[j+1])
						{
							ct[j]++; 
							break; 
						}
    			}
			  for (unsigned int j = 0; j < EBin.size() - 1; j++)
			   {
		  	   cs[j] += ct[j]; 
			   sqSum[j] += pow(ct[j], 2);
			   }
			 }
		  Emean /= nParticles; 
		  // cout<<"The number of particles is "<<nParticles<<endl;
		  cout<<"The average energy of the final jet is "<<Emean<<endl; 
		  cout<<"Finished!"<<endl; 
		  for (unsigned int j=0; j<EBin.size()-1; j++)
		  {
			cs[j] /= nEvents; 
			err[j] = sqrt(((double)sqSum[j]/nEvents - pow(cs[j], 2)) / nEvents); 
			jet_output << (EBin[j] + EBin[j+1]) / 2 << " " << cs[j] / (EBin[j+1] - EBin[j]) << " " << err[j] / (EBin[j+1] - EBin[j]) << endl; 
		  }
		  reader->Close(); 
		}
	}
  }
}
