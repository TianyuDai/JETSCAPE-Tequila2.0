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
  
  auto reader=make_shared<JetScapeReaderAscii>("PGun_16GeV_3fm_elas.dat");
  std::ofstream jet_output ("../Result/PGun/PGun_16GeV_3fm_ElasSpec.txt"); 
  
  const int nEvents = 100000; 
  int jetNum = 0; 
  double Emean = 0.; 
  double jetpTMin = 4., jetRadius = 0.4; 
  fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
  vector <fjcore::PseudoJet> fjInputs; 
  
  vector <double> EBin(101); 
  for (unsigned int j = 0; j < EBin.size(); j++)
  	EBin[j] = 4. + (double) j * 16. / 100.; 
  vector <double> cs(EBin.size()-1, 0.), err(EBin.size()-1, 0.); 
  vector <int> sqSum(EBin.size()-1, 0); 
  
  while (!reader->Finished())
    {
      reader->Next();
      vector <int> ct(EBin.size()-1, 0); 
      cout << "The number of final partons "<< reader->GetPartons().size() << endl; 
      vector <fjcore::PseudoJet> inclusiveJets, sortedJets; 
      fjcore::ClusterSequence clustSeq(reader->GetPartonsForFastJet(), jetDef);
      inclusiveJets = clustSeq.inclusive_jets(jetpTMin); 
      sortedJets    = fjcore::sorted_by_pt(inclusiveJets); 
      jetNum += sortedJets.size(); 
      for (unsigned int i = 0; i < sortedJets.size(); i++)
      {
        Emean += sortedJets[i].perp(); 
      	for (unsigned int j = 0; j < EBin.size() - 1; j++)
      		if (sortedJets[i].perp()>=EBin[j] && sortedJets[i].perp()<EBin[j+1])
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
     }
  Emean /= jetNum; 
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
