#include <iostream>
#include <fstream>
#include <time.h>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

// User modules derived from jetscape framework clasess
#include "TrentoInitial.h"
#include "AdSCFT.h"
#include "Matter.h"
#include "Martini.h"
#include "Tequila.h"
#include "Brick.h"
#include "GubserHydro.h"
#include "HydroFromFile.h"
#include "PGun.h"
#include "PythiaGun.h"

#include <chrono>
#include <thread>

using namespace Jetscape;
// Forward declaration
void Show();

// -------------------------------------
void RunEvents(int N)
{

	//modify the init.xml file
  	  JetScapeXML::Instance()->OpenXMLFile("TrentoPythia.xml");

	  // JetScapeXML::Instance()->CloseXMLFile();
	  
	  double deltaT = 0.01; 
	  tinyxml2::XMLElement *dtxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("deltaT" );
  	  dtxml->SetText(deltaT);

	// vector <double> pTHatBin{20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.}; 
	vector <double> pTHatBin{160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.}; 
	// vector <double> pTHatBin{5., 1380.};

	for (unsigned int iBin = 0; iBin < pTHatBin.size() - 1; iBin++)
	{
		JetScapeXML::Instance()->OpenXMLFile("TrentoPythia.xml");
	  	tinyxml2::XMLElement *pthatminxml=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard" )->FirstChildElement("PythiaGun" )->FirstChildElement("pTHatMin"); 
	  	pthatminxml->SetText(pTHatBin[iBin]);
		tinyxml2::XMLElement *pthatmaxxml=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard" )->FirstChildElement("PythiaGun" )->FirstChildElement("pTHatMax"); 
	  	pthatmaxxml->SetText(pTHatBin[iBin+1]);
		auto jetscape = make_shared<JetScape>("TrentoPythia.xml", N); 
		jetscape->SetReuseHydro (true);
		jetscape->SetNReuseHydro (1000000);

		// Initial conditions and hydro
		auto pythiaGun = make_shared<PythiaGun> ();
		auto hydro = make_shared<HydroFromFile> ();
		jetscape->Add(pythiaGun);
		jetscape->Add(hydro);

		auto jlossmanager = make_shared<JetEnergyLossManager> ();
		auto jloss = make_shared<JetEnergyLoss> ();
		auto matter = make_shared<Matter> (); 

		jloss->Add(matter); 
		jlossmanager->Add(jloss);  
		jetscape->Add(jlossmanager);

		// Hadronization
		// This helper module currently needs to be added for hadronization.
		auto printer = make_shared<PartonPrinter> ();
		jetscape->Add(printer);

		// Output
		auto writer= make_shared<JetScapeWriterAscii> (("./pp2760/"+std::to_string(pTHatBin[iBin])+".dat").c_str());

		jetscape->Add(writer);

		jetscape->Init();
		jetscape->Exec();
		jetscape->Finish();
	}
}

int main(int argc, char** argv)
{
	clock_t t; t = clock();
	time_t start, end; time(&start);

	cout<<endl;

	// DEBUG=true by default and REMARK=false
	// can be also set also via XML file (at least partially)
	JetScapeLogger::Instance()->SetInfo(true);
	JetScapeLogger::Instance()->SetDebug(false);
	JetScapeLogger::Instance()->SetRemark(false);
	//SetVerboseLevel (9 adds a lot of additional debug output)
	//If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevle(9) or 10
	JetScapeLogger::Instance()->SetVerboseLevel(0);

  	Show();

	RunEvents(10000); 

	INFO_NICE<<"Finished!";
	cout<<endl;

	t = clock() - t;
	time(&end);
	printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
	printf ("Real time: %f seconds.\n",difftime(end,start));
/*	printf ("Init time: %f seconds.\n",((float)t1)/CLOCKS_PER_SEC);
	printf ("Exec time: %f seconds.\n",((float)t2)/CLOCKS_PER_SEC);
	printf ("Finish time: %f seconds.\n",((float)t3)/CLOCKS_PER_SEC);*/
	return 0;
}

// -------------------------------------

void Show()
{
	INFO_NICE<<"------------------------------------";
	INFO_NICE<<"| Brick Test JetScape Framework ... |";
	INFO_NICE<<"------------------------------------";
	INFO_NICE;
}
