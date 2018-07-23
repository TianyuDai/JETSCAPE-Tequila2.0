#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"
#include "fjcore.hh"

#include <GTL/dfs.h>

using namespace std;
using namespace fjcore;
using namespace Jetscape;

int main(int argc, char** argv)
{
	JetScapeLogger::Instance()->SetDebug(true);
	JetScapeLogger::Instance()->SetRemark(true);
	JetScapeLogger::Instance()->SetVerboseLevel(0);
  
	auto reader=make_shared<JetScapeReaderAscii>("test_out.dat");
	std::ofstream dist_output ("MartiniHadronsJet.txt"); //Format is SN, PID, E, Px, Py, Pz, Eta, Phi
	vector < shared_ptr<Hadron> > hadrons;
	while (!reader->Finished())
	{
		reader->Next(); 
		hadrons = reader->GetHadrons();
		cout<<"Number of hadrons is: " << hadrons.size() << endl;
		for(unsigned int i=0; i<hadrons.size(); i++) 
		{
			dist_output<<i<<" "<<hadrons[i].get()->pid()<<" "<<hadrons[i].get()->pstat()<<" "<< hadrons[i].get()->e() << " "<< hadrons[i].get()->px()<< " "<< hadrons[i].get()->py() << " "<< hadrons[i].get()->pz()<<  endl;
		}
	}
	reader->Close();  
}
