#include <iostream>
#include <fstream>
#include <time.h>
std::ofstream fout("./PbPb2760/sigmaGen"); 

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
#include "LLBT.h"
#include "Brick.h"
#include "GubserHydro.h"
#include "HydroFromFile.h"
#include "PGun.h"
#include "PythiaGun.h"

#include <chrono>
#include <thread>

using namespace Jetscape;
bool flag = true; 
// Forward declaration
void Show();

// -------------------------------------
void RunEvents(double scale, double alpha_s, double omegacut, int N)
{

	//modify the init.xml file
  	  JetScapeXML::Instance()->OpenXMLFile("TrentoPythia.xml");
  	  tinyxml2::XMLElement *scalexml=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("LLBT" )->FirstChildElement("mu_scale"); 
  	  scalexml->SetText(scale);
  	  tinyxml2::XMLElement *alphaxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("LLBT" )->FirstChildElement("alpha_s");
  	  alphaxml->SetText(alpha_s);
  	  tinyxml2::XMLElement *omegaxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("LLBT" )->FirstChildElement("omega_over_T_cutoff");
  	  omegaxml->SetText(omegacut); 

	  // JetScapeXML::Instance()->CloseXMLFile();
	  
	  double deltaT = 0.01; 
	  tinyxml2::XMLElement *dtxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("deltaT" );
  	  dtxml->SetText(deltaT);

	// vector <double> pTHatBin{400., 600., 800., 1000., 1380.}; 
	// vector <double> pTHatBin{5., 10., 20., 40., 60., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.}; 
	vector <double> pTHatBin{5., 1380.};

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
		// auto trento = make_shared<TrentoInitial> ();
		auto pythiaGun = make_shared<PythiaGun> ();
		// auto pGun = make_shared<PGun> ();
		auto hydro = make_shared<HydroFromFile> ();
		// auto hydro = make_shared<Brick> ();
		// jetscape->Add(trento);
		jetscape->Add(pythiaGun);
		// jetscape->Add(pGun);
		jetscape->Add(hydro);

		// Energy loss
		auto jlossmanager = make_shared<JetEnergyLossManager> ();
		auto jloss = make_shared<JetEnergyLoss> ();
		auto matter = make_shared<Matter> (); 
		auto llbt = make_shared<LLBT> ();
		// auto martini = make_shared<Martini> (); 

		jloss->Add(matter); 
		jloss->Add(llbt);
		// jloss->Add(martini); 
		jlossmanager->Add(jloss);  
		jetscape->Add(jlossmanager);


		// Hadronization
		// This helper module currently needs to be added for hadronization.
		auto printer = make_shared<PartonPrinter> ();
		jetscape->Add(printer);

		// Output
		auto writer= make_shared<JetScapeWriterAscii> (("./test-PbPb2760/qperp"+std::to_string((int)scale)+"omega"+std::to_string((int)omegacut)+"/"+std::to_string(pTHatBin[iBin])+".dat").c_str());

		jetscape->Add(writer);

		jetscape->Init();
		jetscape->Exec();
		jetscape->Finish();
		if (flag)
			fout << pTHatBin[iBin] << " " << pythiaGun->GetSigmaGen() << " " << pythiaGun->GetSigmaErr() << endl; 
	}
	flag = false; 
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

	double scale_list[1] = {1.}; 
  	double alpha_list[1] = {0.3}; 
  	double omegacut_list[3] = {0.5, 1., 2.}; 
	for (int i = 0; i < 1; i++)
  	{
  	  
	  double scale = scale_list[i]; 
	  for (int j = 0; j < 1; j++)
	  {
	  	double alpha_s = alpha_list[j]; 
	  	for (int k = 0; k < 1; k++)
	  	{
	  		double omegacut = omegacut_list[k]; 
	  		RunEvents(scale, alpha_s, omegacut, 2); 
	  	}
	  }
  	}

/*		cout << "LLBT elastic time " << llbt->energy_loss_time_elas << " inelastic time " << llbt->energy_loss_time_inel << endl; 
		cout << "elastic omega sampling time " << llbt->omega_sample_time << " elastic qperp sampling time " <<  llbt->qperp_sample_time; 
		cout << " inelastic omega sampling time " << llbt->inel_sample_time << endl; 
		cout << "tequila time " << llbt->tequila_time << endl;
		cout << "martini time " << martini->martini_time << endl;
		// cout << "pythia time " << pythiaGun->pythia_time << endl;  
		cout << "matter time " << matter->matter_time << endl; 
		cout << "high vir time " << matter->matter_high_time << endl; 
		cout << "low vir time " << matter->matter_low_time << endl; 
		cout << "matter pre time " << matter->matter_pre_time << endl; 
		cout << "matter before time " << matter->matter_before_time << endl; 
		// cout << "hydro time " << hydro->hydro_time << endl; 
		cout << "energy loss time " << jloss->energy_loss_time << endl;
		cout << "manager time " << jlossmanager->manager_time << endl; 
		cout << "number of high vir particles in Matter: " << matter->k << endl; 
		cout << "radiated high vir num: " << llbt->vir_num << endl; 
		// For the future, cleanup is mostly already done in write and clear
		t3 = clock(); 
		t3 = clock() - t3; */

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
