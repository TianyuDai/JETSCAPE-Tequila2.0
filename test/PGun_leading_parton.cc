/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

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
  
  auto reader=make_shared<JetScapeReaderAscii>("PGun_16GeV_1fm_elas_omega001_noangle.dat");
  std::ofstream jet_output ("../Result/PGun/PGun_16GeV_1fm_ElasSpec_leading_omega001_noangle.txt"); 
  
  const int nEvents = 100000; 
  double Emean = 0.; 
  double jetpTMin = 4., jetRadius = 0.4; 
  fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
  vector <fjcore::PseudoJet> fjInputs; 
  
  vector <double> EBin(51); 
  for (unsigned int j = 0; j < EBin.size(); j++)
  	EBin[j] = 4. + (double) j * 16. / 50.; 
  vector <double> cs(EBin.size()-1, 0.), err(EBin.size()-1, 0.); 
  vector <int> sqSum(EBin.size()-1, 0); 
  
  while (!reader->Finished())
    {
      reader->Next();
      vector <int> ct(EBin.size()-1, 0); 
      vector <shared_ptr <Parton> > partons; 
      partons = reader->GetPartons(); 
      double energy = 0.; 
      for (unsigned int i = 0; i < partons.size(); i++)
	if (partons[i]->e() > energy)	energy = partons[i]->e(); 
      Emean += energy; 
      for (unsigned int j = 0; j < EBin.size() - 1; j++)
      	if (energy>=EBin[j] && energy<EBin[j+1])
	{
		ct[j]++; 
		break; 
	}
      for (unsigned int j = 0; j < EBin.size() - 1; j++)
       {
  	   cs[j] += ct[j]; 
	   sqSum[j] += pow(ct[j], 2);
       }
     }
  Emean /= nEvents; 
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
