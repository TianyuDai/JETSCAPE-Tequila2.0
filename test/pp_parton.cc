#include <iostream>
#include <fstream>
#include <time.h>

#include "JetScape.h"
#include "JetScapeWriterStream.h"
#include "PythiaReader.h"

#include "PythiaGun.h"
#include "fjcore.hh"
#include "HardProcess.h"

#include <chrono>
#include <thread>

using namespace std;
using namespace Jetscape;
ofstream fout("test.out"); 

int main(int argc, char** argv)
{
  clock_t t; t = clock();
  time_t start, end; time(&start);

  JetScapeLogger::Instance()->SetInfo(true);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false); 
  JetScapeLogger::Instance()->SetVerboseLevel(0);

  const int nEvents = 10000;  
  vector <double> pTBin{80., 110., 160., 210., 260., 310., 400., 500., 600., 800.}; 
  vector <double> jet_cs(pTBin.size()-1, 0.), err(pTBin.size()-1, 0.); 
  vector <int> sqSum(pTBin.size()-1, 0); 

  double jetRadius = 0.4, jetrapMax = 2.8, jetpTMin = 60.; 
  fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
  vector <fjcore::PseudoJet> fjInputs; 
  fjcore::Selector select_rapidity = fjcore::SelectorRapMax(jetrapMax); 

  auto jetscape = make_shared<JetScape>("./pp_parton.xml", nEvents); 
  jetscape->SetId("primary");

  auto pythiaGun= make_shared<PythiaGun> ();
  jetscape->Add(pythiaGun);
 
  auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");
  jetscape->Add(writer);

  jetscape->Init();

  jetscape->Exec();
  
  jetscape->Finish(); 
  
  auto reader=make_shared<JetScapeReaderAscii>("test_out.dat");
  for (int iEvent = 0; iEvent < nEvents; iEvent++)
  {
	//auto reader=make_shared<JetScapeReaderAscii>("test_out.dat"); 
	int p = reader->Next(iEvent); 
	//auto plist = reader->GetPartonsForFastJet(); 
	//fout<<plist.size()<<endl; 
	//reader->Reset(); 
	//fout<<reader->GetSigmaGenFromFile(nEvents-1)<<endl; 
	vector <int> jet_ct(pTBin.size()-1, 0); 
    	vector <fjcore::PseudoJet> inclusiveJets, sortedJets; 
    	fjcore::ClusterSequence clustSeq(reader->GetPartonsForFastJet(), jetDef);
    	inclusiveJets = clustSeq.inclusive_jets(jetpTMin); 
    	vector <fjcore::PseudoJet> selected_jets = select_rapidity(inclusiveJets); 
    	sortedJets    = fjcore::sorted_by_pt(selected_jets); 
  	if (sortedJets.size() > 1)
  	{
		for (unsigned int j=0; j<pTBin.size()-1; j++)
			if (sortedJets[0].perp()>pTBin[j] && sortedJets[0].perp()<pTBin[j+1])
			{
				jet_ct[j]++; 
				break; 
			}
  	}
  	for (unsigned int j = 0; j < pTBin.size() - 1; j++)
  	{
  		jet_cs[j] += jet_ct[j]; 
  		sqSum[j] += pow(jet_ct[j], 2);
  	} 
  	reader->Reset(p); 
  }
  reader->Close(); 
  
  double sigmapb_weight = pythiaGun->GetSigmaGen() * 1.0e9 / nEvents; 
  for (unsigned int j=0; j<pTBin.size()-1; j++)
  {
	jet_cs[j] = jet_cs[j] * sigmapb_weight; 
	err[j] = sqrt((sqSum[j]*pow(sigmapb_weight, 2)*nEvents - pow(jet_cs[j], 2)) / nEvents); 
	fout << (pTBin[j] + pTBin[j+1]) / 2 << " " << jet_cs[j] / (pTBin[j+1] - pTBin[j]) << " " << err[j] / (pTBin[j+1] - pTBin[j])  << endl; 
  }
  
  INFO_NICE<<"Finished!";
  cout<<endl;

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));
   
  return 0;
}

