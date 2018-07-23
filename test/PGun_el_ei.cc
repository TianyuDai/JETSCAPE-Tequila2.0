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
  
  auto reader=make_shared<JetScapeReaderAscii>("PGun_30GeV_5fm_collisional.dat");
  //std::ofstream jet_output ("../Result/PGun_5GeV_5fm_200MeV.txt"); 
  
  const int nEvents = 10000; 

  double npT = 0., sq = 0., ipT = 30.; 
  while (!reader->Finished())
    {
      reader->Next();
      //vector <int> ct(pTBin.size()-1, 0); 
      vector <shared_ptr <Parton> > partons = reader->GetPartons(); 
      cout << "The number of final partons "<< reader->GetPartons().size() << endl; 
      double fpT = 0.; 
      int k = 0; 
      for (unsigned int i = 0; i < partons.size(); i++)
      	if (abs(partons[i]->pid()) <= 3 && partons[i]->perp() > fpT)
      	{
      		fpT = partons[i]->perp(); 
      		k = i; 
      	}
      double dpT = ipT - fpT; 
      npT += dpT; 
      sq += pow(ipT, 2); 
     }
	
  cout<<"Finished!"<<endl; 
  
  npT /= nEvents; 
  double err = sqrt((sq / nEvents - pow(npT, 2)) / nEvents); 
  cout << npT << " " << err << endl; 
  reader->Close(); 
}
