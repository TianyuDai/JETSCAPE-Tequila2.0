#include "Langevin.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "tinyxml2.h"
#include "Srandom.h"
#include "FluidDynamics.h"

#define MAGENTA "\033[35m"
#define hbarc 0.197327053

using namespace Jetscape; 

const double QS = 1.0; 

Langevin::Langevin()
{
  SetId("Langevin");
  VERBOSE(8);
}

Langevin::~Langevin()
{
  VERBOSE(8);
}

double const tiny = 1e-10;

void Langevin::Init()
{
	INFO<<"Intialize Langevin ...";
	
  	// Redundant (get this from Base) quick fix here for now
  	tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  	if ( !eloss )     throw std::runtime_error("Eloss not properly initialized in XML file ...");

  	tinyxml2::XMLElement *langevin=eloss->FirstChildElement("Langevin");
  	// check that all is there
  	if ( !langevin )     throw std::runtime_error("Langevin not properly initialized in XML file ...");
  	if ( !langevin->FirstChildElement( "name" ) )     throw std::runtime_error("Langevin not properly initialized in XML file ...");
  	if ( !langevin->FirstChildElement( "alpha_s" ) )     throw std::runtime_error("Langevin not properly initialized in XML file ...");
  	if ( !langevin->FirstChildElement( "mass" ) )     throw std::runtime_error("Langevin not properly initialized in XML file ...");

  	string s = langevin->FirstChildElement( "name" )->GetText();
  	JSDEBUG << s << " to be initilizied ...";

  	alpha_s = 0.3;
  	langevin->FirstChildElement("alpha_s")->QueryDoubleText(&alpha_s);
    
  	M = 0.;
  	langevin->FirstChildElement("mass")->QueryDoubleText(&M);

	g = sqrt(4.*M_PI*alpha_s);
}

void Langevin::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
	VERBOSESHOWER(5)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<< pIn.size();
	int Id, newId;
  	double pAbs, px, py, pz;   // momentum for initial parton (pIn)
	double xx, yy, zz;         // position of initial parton (pIn)
  	FourVector pVec, pVecNew;  // 4 vectors for momenta before & after process
  	FourVector xVec;           // 4 vector for position (for next time step!)
  	double eta;                // pseudo-rapidity
  	
	// flow info
  	double vx, vy, vz;         // 3 components of flow velocity
	double T;                  // Temperature of fluid cell

	for (int i=0;i<pIn.size();i++) {
    	// Only accept low t particles
    		// if (pIn[i].t() > QS*QS + rounding_error) continue;
    		TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.
    
    		Id = pIn[i].pid();

    		px = pIn[i].px();
    		py = pIn[i].py();
    		pz = pIn[i].pz();
   		// In Langevin, particles are all massless and on-shell
    		pAbs = sqrt(px*px+py*py+pz*pz+M*M);
    		pVec = FourVector ( px, py, pz, pAbs );

    		xx = pIn[i].x_in().x();
    		yy = pIn[i].x_in().y();
    		zz = pIn[i].x_in().z();

    		eta = pIn[i].eta();

    		std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
    		GetHydroCellSignal(Time, xx, yy, zz, check_fluid_info_ptr);
    		VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "
	      		<<check_fluid_info_ptr->temperature;

    		vx = check_fluid_info_ptr->vx;
    		vy = check_fluid_info_ptr->vy;
    		vz = check_fluid_info_ptr->vz;
    		T = check_fluid_info_ptr->temperature;
    		
    		M = alpha_s * T / 3.; 

		xVec = FourVector( xx+px/pAbs*deltaT, yy+py/pAbs*deltaT, zz+pz/pAbs*deltaT, Time+deltaT );
		VERBOSE(8)<< MAGENTA
	      		<< "Time = " << Time << " Id = " << Id << " T = " << T
	      		<< " pAbs = " << pAbs << " " << px << " " << py << " " << pz 
	      		<< " | position = " << xx << " " << yy << " " << zz;
	      	fourvec FS; 
	      	Ito_update(deltaT/hbarc, T, {vx, vy, vz}, fourvec{pVec.t(), pVec.x(), pVec.y(), pVec.z()}, FS); 
	      	pVecNew = FourVector(FS.x(), FS.y(), FS.z(), FS.t()); 
	      	pOut.push_back(Parton(0, Id, 0, pVecNew, xVec)); 
	      	pOut[pOut.size()-1].set_form_time(0.);
	      	return; 
	 }
}

double Langevin::qpara(double E, double T){
	double mD = sqrt(std::pow(g*T, 2)*(nc/3. + nf/6.)); 
	double mu = 2.*g*T; 
	double Minf = sqrt(pow(mD, 2)/2.); 
	return std::pow(g*Minf, 2)*CR*T/(2.*M_PI)*log(mu/Minf);
	// return std::pow(T, 3) * 5.; 
}

double Langevin::qperp(double E, double T){
	double mD = sqrt(std::pow(g*T, 2)*(nc/3. + nf/6.)); 
	double mu = 2.*g*T; 
	return std::pow(g*mD, 2) * CR * T / (2.*M_PI) * log(mu/mD); 
	// return std::pow(T, 3) * 5.; 
}

void Langevin::Ito_update(double dt, double T, std::vector<double> v, 
						const fourvec & pIn, fourvec & pOut){
	// Boost pIn to medium frame
	auto pIn_cell = pIn.boost_to(v[0], v[1], v[2]); 
	// imaging rotating to a frame where pIn lies on z-axis
	double E0 = pIn_cell.t();
	double p0 = std::sqrt(E0*E0 - M*M + 1e-9);
	double kt = qperp(E0, T) / 2.;
	double kl = qpara(E0, T); 
	double drag = kl/(2.*p0*T)+1./(2.*p0*p0)*(kt*2.-kl*2.); //eta_D
		   
    	double white_noise_holder[3];
    	for (size_t i=0; i<3; ++i) 
		white_noise_holder[i] = Srandom::white_noise(Srandom::gen);

	double Ct = std::sqrt(kt*dt);
	double Cl = std::sqrt(kl*dt);

    	pOut.a[1] = Ct * white_noise_holder[0];
    	pOut.a[2] = Ct * white_noise_holder[1];
    	pOut.a[3] = p0 * (1. - drag * dt) + Cl * white_noise_holder[2];
    	pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
					+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// rotate back to the original frame
	pOut = pOut.rotate_back(pIn_cell);
	// boost back to lab frame
	pOut = pOut.boost_back(v[0], v[1], v[2]);
}


