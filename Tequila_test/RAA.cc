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
	const int nEvents = 10000; 
	// double Emean = 0.; 

	vector <double> pTBin{50., 63., 79., 100., 125., 158., 199., 251., 316.};
    // vector <double> pTHatBin{5., 10., 50., 100., 200., 400., 600., 800., 1000., 1380.};
	vector <double> pTHatBin{20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.}; 
	// vector <double> pTHatBin{5., 1380.}; 
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

    double jetpTMin = 50., jetRadius = 0.4, jetAbsRapMax = 2.8, partonpTMin = 20.; 
    double deltaRap = jetAbsRapMax * 2.; 
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
    vector <fjcore::PseudoJet> fjInputs; 
    fjcore::Selector select_rapidity = fjcore::SelectorAbsRapMax(jetAbsRapMax); 
	fjcore::Selector select_pt = fjcore::SelectorPtMin(partonpTMin); 
	fjcore::Selector select_both = select_rapidity && select_pt; 

	double scale_list[1] = {1.}; 
  	double alpha_list[1] = {0.3}; 
  	double omegacut_list[2] = {0.5, 2.}; 
	for (int i = 0; i < 1; i++)
  	{
  	  
	  double scale = scale_list[i]; 
	  for (int j = 0; j < 1; j++)
	  {
	  	double alpha_s = alpha_list[j]; 
	  	for (int k = 0; k < 2; k++)
	  	{
	  		double omegacut = omegacut_list[k]; 
			std::ofstream jet_output (("../../Result/Tequila/AA/PbPb2760_omega"+std::to_string(omegacut)+"qperp"+std::to_string(scale)+".txt").c_str());
   			std::vector <double> jet_cs(pTBin.size()-1, 0.), err(pTBin.size()-1, 0.), sqSum(pTBin.size()-1, 0.); 
    		std::vector <int> jet_ct(pTBin.size()-1);
			for (unsigned int iBin = 0; iBin < pTHatBin.size() - 1; iBin++)
			{
				auto reader=make_shared<JetScapeReaderAscii>(("./PbPb2760/qperp"+std::to_string((int)scale)+"omega"+std::to_string((int)omegacut)+"/"+std::to_string(pTHatBin[iBin])+".dat").c_str());  
				jet_ct = std::vector<int>(pTBin.size()-1, 0); 
				int k = 0; 
				while (!reader->Finished())
				{
					k++; 
					reader->Next();
				    fjInputs.resize(0);          
					vector <fjcore::PseudoJet> inclusiveJets, sortedJets; 
				    fjcore::ClusterSequence clustSeq(reader->GetPartonsForFastJet(), jetDef); 
					// std::cout << reader->GetPartonsForFastJet().size() << "\n"; 
					inclusiveJets = clustSeq.inclusive_jets(jetpTMin); 
				    vector <fjcore::PseudoJet> selected_jets = select_both(inclusiveJets); 
				    sortedJets = fjcore::sorted_by_pt(selected_jets); 
					// std::cout << sortedJets.size() << "\n"; 
				    for (unsigned int iJet = 0; iJet < sortedJets.size(); iJet++)
				    {
						// std::cout << sortedJets[iJet].perp() << "\n"; 
				        for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
				            if (sortedJets[iJet].perp() > pTBin[ipT] && sortedJets[iJet].perp() < pTBin[ipT+1])
				            {
								if (iBin==0 && ipT != 0) std::cout << "event " << k; 
				                jet_ct[ipT]++; 
								// std::cout << jet_ct[iJet] << "\n"; 
				                break; 
				            }
				    }
				}
				std::cout << "pTHatBin " << iBin << " jet_abnormal " <<  jet_ct[3] << " jet normal " << jet_ct[0]<< " " << jet_ct[1] << " " << jet_ct[2]<< " " << jet_ct[4] << " " << jet_ct[5] << endl; 
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
		}
	  }
	}


	cout<<"Finished!"<<endl; 
}
