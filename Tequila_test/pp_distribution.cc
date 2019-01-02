#include <iostream>
#include <fstream>
std::ifstream fin("PbPb2760/sigmaGen"); 
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
	const int nEvents = 1000; 

	vector <double> pTBin{50., 63., 79., 100., 125., 158., 199., 251., 316.};
    // vector <double> pTHatBin{5., 10., 50., 100., 200., 400., 600., 800., 1000., 1380.};
	vector <double> pTHatBin{5., 10., 20., 40., 60., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.}; 
    vector <double> sigmaGen(pTHatBin.size()-1);  
    vector <double> sigmaErr(pTHatBin.size()-1);  
	// std::vector <double> jet_cs(pTBin.size()-1, 0.), err(pTBin.size()-1, 0.), sqSum(pTBin.size()-1, 0.); 
    //std::vector <int> jet_ct(pTBin.size()-1);
    for (unsigned int iSigma = 0; iSigma < pTHatBin.size() - 1; iSigma++)
    {
        double pTHat; 
        fin >> pTHat; 
        if (pTHat == pTHatBin[iSigma])
            fin >> sigmaGen[iSigma] >> sigmaErr[iSigma]; 
        else std::cout << "pTHatBin boundary does not match! "; 
    } 

    double jetpTMin = 10., jetRadius = 0.4, jetAbsRapMax = 2.8; 
    double deltaRap = jetAbsRapMax * 2.; 
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
    vector <fjcore::PseudoJet> fjInputs; 
    fjcore::Selector select_rapidity = fjcore::SelectorAbsRapMax(jetAbsRapMax); 

	std::ofstream jet_output ("../../Result/Tequila/AA/pp2760_hydro.txt");
   	std::vector <double> jet_cs(pTBin.size()-1, 0.), err(pTBin.size()-1, 0.), sqSum(pTBin.size()-1, 0.); 
    std::vector <int> jet_ct(pTBin.size()-1);
	for (unsigned int iBin = 0; iBin < pTHatBin.size() - 1; iBin++)
	{
		auto reader=make_shared<JetScapeReaderAscii>(("./pp2760/"+std::to_string(pTHatBin[iBin])+".dat").c_str());  
		jet_ct = std::vector<int>(pTBin.size()-1, 0); 
		while (!reader->Finished())
		{
			reader->Next();
			fjInputs.resize(0);          
			vector <fjcore::PseudoJet> inclusiveJets, sortedJets; 
			fjcore::ClusterSequence clustSeq(reader->GetPartonsForFastJet(), jetDef); 
			inclusiveJets = clustSeq.inclusive_jets(jetpTMin); 
			// std::cout << inclusiveJets.size() << endl; 
			vector <fjcore::PseudoJet> selected_jets = select_rapidity(inclusiveJets); 
			sortedJets = fjcore::sorted_by_pt(selected_jets); 
			for (unsigned int iJet = 0; iJet < sortedJets.size(); iJet++)
			{
				for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
					if (sortedJets[iJet].perp() > pTBin[ipT] && sortedJets[iJet].perp() < pTBin[ipT+1])
				    {
						// std::cout << "a jet" << endl; 
				    	jet_ct[ipT]++; 
				        break; 
				    }
			}
		}
		double sigma = sigmaGen[iBin] / nEvents; 
		for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
		{
			jet_cs[ipT] += jet_ct[ipT] * sigma; 
			sqSum[ipT] += jet_ct[ipT] * pow(sigma, 2);
		}
		reader->Close(); 
	}
	for (unsigned int ipT = 0; ipT < pTBin.size()-1; ipT++)
    {
        		err[ipT] = jet_cs[ipT] / sqrt(pow(jet_cs[ipT], 2) / sqSum[ipT]); 
				jet_output << (pTBin[ipT] + pTBin[ipT+1]) / 2 << " " << jet_cs[ipT] / (pTBin[ipT+1] - pTBin[ipT]) << " " << err[ipT] / (pTBin[ipT+1] - pTBin[ipT]) << endl; 
	}

	cout<<"Finished!"<<endl; 
}
