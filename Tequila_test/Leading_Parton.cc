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
  double scale_list[1] = {1.}; 
  double alpha_list[1] = {0.3}; 
  double omegacut_list[1] = {1.}; 

  double jetpTMin = 20., jetRadius = 0.4, partonpTMin = 20.; 

  fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
  vector <fjcore::PseudoJet> fjInputs; 
  fjcore::Selector select_pt = fjcore::SelectorPtMin(partonpTMin);

  for (int ii = 0; ii < 1; ii++)
  {
  	  
	  double scale = scale_list[ii]; 
	  for (int jj = 0; jj < 1; jj++)
	  {
	  	double alpha_s = alpha_list[jj]; 
		for (int kk = 0; kk < 1; kk++)
		{
		  double omegacut = omegacut_list[kk]; 
		  auto reader=make_shared<JetScapeReaderAscii>("5GeV_gluon_muscale"+std::to_string(scale)+"alpha"+std::to_string(alpha_s)+"omega"+std::to_string(omegacut)+"_qgpbrick_elas_addConv.dat"); 
 		  // std::cout << "100GeV_gluon_muscale"+std::to_string(scale)+"alpha"+std::to_string(alpha_s)+"omega"+    std::to_string(omegacut)+"_qgpbrick_elas_put1.dat\n"; 
		  std::ofstream jet_output (("../../../../Result/Tequila2.0/elas/Toverp_test/compare/5GeV_gluon_muqperp"+std::to_string(scale)+"alpha"+std::to_string(alpha_s)+"muomega"+std::to_string(omegacut)+"_1fm_QGP_elas_addConv.txt").c_str()); 

		  while (!reader->Finished())
		  { 
			  reader->Next();
			  fjInputs.resize(0);          
		  	  vector <fjcore::PseudoJet> inclusiveJets, sortedJets; 
		  	  fjcore::ClusterSequence clustSeq(reader->GetPartonsForFastJet(), jetDef); 
   		  	  inclusiveJets = clustSeq.inclusive_jets(jetpTMin); 
			  vector <fjcore::PseudoJet> selected_jets = select_pt(inclusiveJets); 
		  	  sortedJets = fjcore::sorted_by_pt(selected_jets);
			  vector <shared_ptr <Parton> > partons; 
			  partons = reader->GetPartons(); 
			  double pT = 0.; 
			  double energy = 0.; 
			  double px;
			  double py; 
			  double pz; 
                          int id; 
                          int stat;  
			  for (unsigned int i = 0; i < partons.size(); i++)
			  {
			   	// if  (partons[i]->e() > energy && partons[i]->pid() == 21) {energy = partons[i]->e(); pT = partons[i]->pt(); px = partons[i]->px(); py = partons[i]->py(); pz = partons[i]->pz(); } 
			   	// if  (partons[i]->pt() > pT && partons[i]->pid() == 21)
				//	pT = partons[i]->pt();
				px = partons[i]->px(); 
                                py = partons[i]->py(); 
                                id = partons[i]->pid(); 
                                stat = partons[i]->pstat(); 
                                energy = partons[i]->e(); 
				if (abs(partons[i]->e() - 0.) > 1e-5)// jet_output << partons[i]->e() << "\n"; 
                                    jet_output << "energy: " << energy << ", px: " << px << ", py: " << py << ", id: " << id << ", stat: " << stat << endl; 
			  }
			  // if (energy != 0.) jet_output << energy << "\n"; 
		  }
		  reader->Close(); 
		}
	}
  }
}
