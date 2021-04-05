#include "Tequila.h"
#include "IntTabulator.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "tinyxml2.h"
#include "Srandom.h"
#include "FluidDynamics.h"

#include "JetScapeWriter.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp2d.h"
#include "gsl/gsl_spline2d.h"

#define nB(k) 1./(exp(k)-1.)
#define nF(k) 1./(exp(k)+1.)

#define MAGENTA "\033[35m"
#define hbarc 0.197327053

using namespace std; 
ofstream fout_rate("dGamma_domega_qqg"); 
// ofstream fout_ratew("omegadGamma_domega"); 

using namespace Jetscape; 
Tequila::Tequila()
{
	SetId("Tequila");
	VERBOSE(8);
}

Tequila::~Tequila()
{
	VERBOSE(8);
}

void Tequila::Finish()
{
	free(za); 
	delete[] rate_ggg_p; 
	delete[] rate_gqq_p;
	delete[] rate_qqg_p;
	for(int ip=0;ip<nb_points_in_p; ip++)
	{
		for(int iomega=0;iomega<nb_points_in_omega; iomega++) {
			 delete[] differential_rate_ggg_p_omega_qperp[ip][iomega];
			 delete[] differential_rate_gqq_p_omega_qperp[ip][iomega];
			 delete[] differential_rate_qqg_p_omega_qperp[ip][iomega];
		}
		delete[] differential_rate_ggg_p_omega_qperp[ip];
		delete[] differential_rate_gqq_p_omega_qperp[ip];
		delete[] differential_rate_qqg_p_omega_qperp[ip];
	}
	delete[] differential_rate_ggg_p_omega_qperp;
	delete[] differential_rate_gqq_p_omega_qperp;
	delete[] differential_rate_qqg_p_omega_qperp;
}

void Tequila::Init()
{
	INFO<<"Intialize Langevin and Linear Boltzmann ...";
	
  	// Redundant (get this from Base) quick fix here for now
  	tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  	if ( !eloss )     throw std::runtime_error("Eloss not properly initialized in XML file ...");

  	tinyxml2::XMLElement *tequila=eloss->FirstChildElement("Tequila");
  	// check that all is there
  	if ( !tequila )     throw std::runtime_error("Tequila not properly initialized in XML file ...");
  	if ( !tequila->FirstChildElement( "name" ) )     throw std::runtime_error("Tequila - name not properly initialized in XML file ...");
  	if ( !tequila->FirstChildElement( "alpha_s" ) )     throw std::runtime_error("Tequila - alpha_s not properly initialized in XML file ...");
  	if ( !tequila->FirstChildElement( "mu_scale" ) )     throw std::runtime_error("Tequila - mu_scale not properly initialized in XML file ...");
  	if ( !tequila->FirstChildElement( "pcut" ) )     throw std::runtime_error("Tequila - pcut not properly initialized in XML file ...");
  	if ( !tequila->FirstChildElement( "Q0" ) )     throw std::runtime_error("Tequila - Q0 not properly initialized in XML file ...");
  	if ( !tequila->FirstChildElement( "hydro_Tc" ) )     throw std::runtime_error("Tequila - hydro_Tc not properly initialized in XML file ...");
  	if ( !tequila->FirstChildElement( "tables_path" ) )     throw std::runtime_error("Tequila - tables_path not properly initialized in XML file ...");
  	if ( !tequila->FirstChildElement( "omega_over_T_cutoff" ) )     throw std::runtime_error("Tequila - omega not properly initialized in XML file ...");


  	string s = tequila->FirstChildElement( "name" )->GetText();
  	JSDEBUG << s << " to be initilizied ...";

  	alpha_s = 0.3;
  	tequila->FirstChildElement("alpha_s")->QueryDoubleText(&alpha_s);

	mu_scale = 1.;
  	tequila->FirstChildElement("mu_scale")->QueryDoubleText(&mu_scale);
    
        pcut = 2.0;
  	tequila->FirstChildElement("pcut")->QueryDoubleText(&pcut);
  	
  	M = 0.;
  	tequila->FirstChildElement("mass")->QueryDoubleText(&M);

	g = sqrt(4.*M_PI*alpha_s); 
	// set muperp / T as muperp
	muperp = mu_scale; 
	
	
  	hydro_Tc = 0.16;
  	tequila->FirstChildElement("hydro_Tc")->QueryDoubleText(&hydro_Tc);
  	
  	Q0 = 1.0;
  	tequila->FirstChildElement("Q0")->QueryDoubleText(&Q0);



 	omega_over_T_cutoff = 1.;
  	tequila->FirstChildElement("omega_over_T_cutoff")->QueryDoubleText(&omega_over_T_cutoff);
	// omega_over_T_cutoff *= sqrt(g); 

  	// Path to additional data
  	path_to_tables=tequila->FirstChildElement( "tables_path" )->GetText();

	INFO << "Load elastic rate... "; 
	// Load elastic rate
	LoadElasticTables(); 
	// INFO << "Finished loading elastic rate! "; 
	
	ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
	
	allocate_memory_for_radiative_rate_table();

  	// Load differential rate
  	INFO << "Loading differential collinear rate...\n";
  	const double alpha_EM=1./137.; //Not currently a parameter
  	const int Nf=3;
  	//const std::string location_of_pretabulated_collinear_rates="./Tequila/";
  	// Either compute or read from file the collinear differential rates and load them into member array "differential_rate_p_omega_qperp[][][]"
  	load_differential_rate(alpha_s, alpha_EM, Nf, path_to_tables);
  	// Compute total rate from differential rate and save in member array "rate_p[]"
  	INFO << "Computing integrated collinear rate...\n";
  	evaluate_integrated_rate(omega_over_T_cutoff, differential_rate_qqg_p_omega_qperp,rate_qqg_p); 
  	evaluate_integrated_rate(omega_over_T_cutoff, differential_rate_ggg_p_omega_qperp,rate_ggg_p);
  	evaluate_integrated_rate(omega_over_T_cutoff, differential_rate_gqq_p_omega_qperp,rate_gqq_p);
  	INFO << "Done computing integrated collinear rate.\n";

}

void Tequila::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
        for (double w = 0.; w < 100./0.3; w += 0.01)
        {
            fout_rate << w << " " << differential_rate(100./0.3, w, 0.1, differential_rate_qqg_p_omega_qperp) << "\n"; 
        }
	VERBOSESHOWER(5)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<< pIn.size();
	int Id, newId = 1;
  	double pAbs, px, py, pz;   // momentum for initial parton (pIn)
  	double pRest, pxRest;      // momentum in the rest frame of fluid cell (pIn)
  	double pyRest, pzRest;

  	double k, kRest = 0.;           // momentum for radiated parton (pOut)
  	double pNewRest;
	double xx, yy, zz;         // position of initial parton (pIn)
  	FourVector pVec, pVecNew, pVecNewest;  // 4 vectors for momenta before & after process
  	FourVector xVec;           // 4 vector for position (for next time step!)
  	FourVector pVecRest; 
  	FourVector kVec, kVecRest, kVecNew, kVecNewest;           // 4 vector for momentum of radiated particle
  	double eta;                // pseudo-rapidity
	double velocity_jet[4];    // jet velocity for MATTER
  	
	// flow info
  	double vx, vy, vz;         // 3 components of flow velocity
	double T;                  // Temperature of fluid cell
	double beta, gamma; 
	double omega = 0., qperp = 0.; 

	for (int i=0;i<pIn.size();i++)
	{
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

		// Only accept low t particles
		if (pIn[i].t() > Q0*Q0 + rounding_error || T < hydro_Tc) continue; 
		TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.

		
		beta = sqrt( vx*vx + vy*vy + vz*vz );
		gamma = 1./sqrt(1.-beta*beta); 
		
		// Set momentum in fluid cell's frame
		fourvec pIn_4vec = fourvec{pAbs, px, py, pz}; 
		fourvec pIn_cell = pIn_4vec.boost_to(vx, vy, vz); 
		pVecRest = FourVector(pIn_cell.x(), pIn_cell.y(), pIn_cell.z(), pIn_cell.t()); 
		pRest = pVecRest.t(); 
		pxRest = pVecRest.x(); 
		pyRest = pVecRest.y(); 
		pzRest = pVecRest.z(); 
		if (pRest < pcut) continue; 

		VERBOSE(8) << MAGENTA
      		<< "Time = " << Time << " Id = " << Id << " T = " << T
      		<< " pAbs = " << pAbs << " " << px << " " << py << " " << pz 
      		<< " | position = " << xx << " " << yy << " " << zz;
      	
      	// Pass the address of deltaT to DetermineProcess, if deltaT is changed, the changed value could be passed back. Calculate xVec after passing back the real deltaT. 
      	// time step should be in the rest frame. 
	      	double deltaTRest = deltaT / gamma; 
		Lambda = std::min({pRest, 2.*T*wMax, 2.*sqrt(3.*pRest*T), 2.*pcut}); 
      		process_type process = DetermineProcess(pRest, T, deltaTRest, Id);
                // if (process != none) std::cout << "process is " << ProcessStrings[process] << "\n";  
      		xVec = FourVector( xx+px/pAbs*deltaT, yy+py/pAbs*deltaT, zz+pz/pAbs*deltaT, Time+deltaT );
		velocity_jet[0]=1.0;
    		velocity_jet[1]=pIn[i].jet_v().x();
    		velocity_jet[2]=pIn[i].jet_v().y();
    		velocity_jet[3]=pIn[i].jet_v().z();

		IntTabulator inttabulator; 
      		if (process == none)
      			pVecNew = pVecRest; 
      		else if (process == gg || process == gq || process == qg || process == qq || process == qqp || process == qqb)
      		{
      			omega = Energy_Transfer(pRest, T, process)*T; 
			qperp = 0.; 
                        // std::cout << omega << "\n"; 
			// qperp is actually \tilde{q_\perp} here, we should convert it to the real q_\perp. 
			// double q2 = omega*omega+qperp*qperp; 
			// double qz = qperp*qperp/(2.*pRest)+omega; 
			// qperp = sqrt(q2 - qz*qz); 
			pVecNew = Momentum_Update(omega, qperp, T, pVecRest); 
		}
      		else if (process == GqQg || process == GgQbq)
      		{
      			omega = Energy_Transfer(pRest, T, process)*T; 
			qperp = 0.; 
			pVecNew = Momentum_Update(omega, qperp, T, pVecRest); 
			// choose the Id of new qqbar pair. Note that we only deal with nf = 3
			double r = ZeroOneDistribution(*GetMt19937Generator());
			if (r < 1./6.) Id = 1;
			else if (r < 2./6.) Id = 2;
			else if (r < 3./6.) Id = 3;
			else if (r < 4./6.) Id = -1; 
			else if (r < 5./6.) Id = -2; 
			else Id = -3; 
		}
      		else if (process == QgGq || process == QbqGg)
      		{
      			omega = Energy_Transfer(pRest, T, process)*T; 
			qperp = 0.; 
			pVecNew = Momentum_Update(omega, qperp, T, pVecRest);
                        Id = 21;  
		}
		else if (process == gg_split || process == gq_split || process == qg_split || process == qq_split || process == qqp_split || process == qqb_split)
      		{
			double k0Rest = thermalDistributionSampling(process)*T; 
			// kRest = thermalDistributionSampling(process)*T; 
			double x = xSampling(pRest, k0Rest, T, process); 
      			// omega = Energy_Transfer_Split(pRest, T, process); 
			// omega = 0.; 
			omega = pRest * x; 
			// omega = (pRest - k0Rest) * x; 
			// double qplus = x * pRest; 
			// qperp = 2.*sqrt(kRest*qplus/pRest*(pRest-qplus));
			qperp = 0.; 
			// qperp = TransverseMomentum_Transfer_Split(pRest, omega, T, process); 
			pVecNew = Momentum_Update(omega, qperp, T, pVecRest); 
			// FourVector k0VecRest; 
			// k0VecRest.Set(-1.*(pxRest/pRest)*k0Rest, -1.*(pyRest/pRest)*k0Rest, -1.*(pzRest/pRest)*k0Rest, k0Rest); 
			// kVec = Momentum_Update(-1.*omega, qperp, T, k0VecRest); 	// Need to be revisited. 
			// kVec.Set( (pxRest/pRest)*kRest, (pyRest/pRest)*kRest, (pzRest/pRest)*kRest, kRest ); 
			if (process == gg_split || process == qg_split) newId = 21; 
			else if (process == gq_split)
			{
				double r = ZeroOneDistribution(*GetMt19937Generator());
				if (r < 1./6.) newId = 1;
				else if (r < 2./6.) newId = 2;
				else if (r < 3./6.) newId = 3;
				else if (r < 4./6.) newId = -1;
				else if (r < 5./6.) newId = -2;
				else newId = -3;
			}
			else if (process == qq_split)
				newId = Id; 
			else if (process == qqp_split)
			{
				double r = ZeroOneDistribution(*GetMt19937Generator());
				if (r < 1./4.) newId = (abs(Id)+1) % 3 + 1;
				else if (r < 2./4.) newId = (abs(Id)+2) % 3 + 1;
				else if (r < 3./4.) newId = -1*((abs(Id)+1) % 3 + 1);
				else newId = -1*((abs(Id)+2) % 3 + 1);
			}
			else if (process == qqb_split)
				newId = -1*Id; 
		}
      		/*else if (process == gqqg || process == gq_inel_conv)
      		{ 
      			double r = ZeroOneDistribution(*GetMt19937Generator());
			if (r < 1./6.) Id = 1;
			else if (r < 2./6.) Id = 2;
			else if (r < 3./6.)Id = 3;
			else if (r < 4./6.)Id = -1;
			else if (r < 5./6.)Id = -2;
			else Id = -3;
			pVecNew = pVecRest; 
		}
		else if (process == qggq || process == qg_inel_conv)
		{
      			Id = 21; 
      			pVecNew = pVecRest;  
      		}
      		else if (process == ggg)
      		{ 
      			if (pRest/T < AMYpCut) return;
			// WARN << process; 
    		// sample radiated parton's momentum
    			sample_dgamma_dwdq(pRest, T,differential_rate_ggg_p_omega_qperp, omega, qperp);
    			kRest = omega; 
        		if(kRest > pRest) return;

    		// final state parton's momentum
    			pNewRest = pRest - kRest;

    		// if pNew is smaller than pcut, final state parton is
    		// absorbed into medium
			pVecNew.Set( (pxRest/pRest)*pNewRest, (pyRest/pRest)*pNewRest, (pzRest/pRest)*pNewRest, pNewRest );

    			kVec.Set( (pxRest/pRest)*kRest, (pyRest/pRest)*kRest, (pzRest/pRest)*kRest, kRest );
			newId = 21; 
      		}
		else if (process == gqq)
		{
			if (pRest/T < AMYpCut) return;

			// sample radiated parton's momentum
			sample_dgamma_dwdq(pRest, T,differential_rate_gqq_p_omega_qperp, omega, qperp);		
			kRest = omega; 
        		if(kRest > pRest) return;

			// final state parton's momentum
			pNewRest = pRest - kRest;
			
			// choose the Id of new qqbar pair. Note that we only deal with nf = 3
			double r = ZeroOneDistribution(*GetMt19937Generator());
			if (r < 1./6.) Id = 1;
			else if (r < 2./6.) Id = 2;
			else if (r < 3./6.) Id = 3;
			else if (r < 4./6.) Id = -1; 
			else if (r < 5./6.) Id = -2; 
			else Id = -3; 
			
			// if pNew is smaller than pcut, final state parton is
			// absorbed into medium
			// *momentum of quark is usually larger than that of anti-quark
			pVecNew.Set( (pxRest/pRest)*pNewRest, (pyRest/pRest)*pNewRest, (pzRest/pRest)*pNewRest, pNewRest );

			kVec.Set( (pxRest/pRest)*kRest, (pyRest/pRest)*kRest, (pzRest/pRest)*kRest, kRest );
			newId = -1 * Id; 
		}
		else if (process == qqg)
		{
			if (pRest/T < AMYpCut) return;

			// sample radiated parton's momentum
			sample_dgamma_dwdq(pRest, T,differential_rate_qqg_p_omega_qperp, omega, qperp);		
			kRest = omega; 
			//kRest = getNewMomentumRad(pRest, T, process);
        		if(kRest > pRest) return;

			// final state parton's momentum
			pNewRest = pRest - kRest;

			// if pNew is smaller than pcut, final state parton is
			// absorbed into medium
			pVecNew.Set( (pxRest/pRest)*pNewRest, (pyRest/pRest)*pNewRest, (pzRest/pRest)*pNewRest, pNewRest );

			kVec.Set( (pxRest/pRest)*kRest, (pyRest/pRest)*kRest, (pzRest/pRest)*kRest, kRest );
			newId = 21; 
		}*/
		else pVecNew = pVecRest; 
      	// pVecNew = pVecRest; 
      		pVecNewest = Langevin_Update(deltaTRest / hbarc, T, pVecNew, Id); 

                // The conversion process for soft elastic interactions
                
                double minf2 = g*g*T*T*CF/4; 
                double dt = deltaTRest / hbarc;
                double soft_q2g = g*g*CF*minf2/8/M_PI/pRest*log(1+muperp*muperp/minf2); 
                double soft_g2q = dF/dA*soft_q2g; 
                if (Id == 21)
                {
                    double rand_conv = ZeroOneDistribution(*GetMt19937Generator()); 
                    if (rand_conv < soft_g2q*dt)
                    {
			double r = ZeroOneDistribution(*GetMt19937Generator());
			if (r < 1./6.) Id = 1;
			else if (r < 2./6.) Id = 2;
			else if (r < 3./6.) Id = 3;
			else if (r < 4./6.) Id = -1; 
			else if (r < 5./6.) Id = -2; 
			else Id = -3; 
                    }
                }
                else
                {
                    double rand_conv = ZeroOneDistribution(*GetMt19937Generator()); 
                    if (rand_conv < soft_q2g*dt) Id = 21; 
                }
                
		fourvec pOut_4vec, kOut_4vec; 
      		if (kVec.t() != 0)
		{
			kVecNewest = kVec; 
      			kVecNewest = Langevin_Update(deltaTRest / hbarc, T, kVec, newId); 
			kOut_4vec = fourvec{kVecNewest.t(), kVecNewest.x(), kVecNewest.y(), kVecNewest.z()}; 
      			kOut_4vec = kOut_4vec.boost_back(vx, vy, vz); 
		}
        // pVecNewest = pVecNew; 
      	pOut_4vec = fourvec{pVecNewest.t(), pVecNewest.x(), pVecNewest.y(), pVecNewest.z()}; 
      	pOut_4vec = pOut_4vec.boost_back(vx, vy, vz); 
        // std::cout << "before pushback " << pOut_4vec.t() << "\n"; 
      	if (pVecNewest.t() > pcut)
      	{
      		pOut.push_back(Parton(0, Id, 0, FourVector(pOut_4vec.x(), pOut_4vec.y(), pOut_4vec.z(), pOut_4vec.t()), xVec)); 
      		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); 
      	}
      	if (kVecNewest.t() > pcut)
      	{
      		pOut.push_back(Parton(0, newId, 0, FourVector(kOut_4vec.x(), kOut_4vec.y(), kOut_4vec.z(), kOut_4vec.t()), xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); 
      	}
      	return; 
	}
}

process_type Tequila::DetermineProcess(double pRest, double T, double deltaTRest, int Id)
{
	double dt = deltaTRest / hbarc; 
	double rateTotal; 
	double rate[nTotalProcess] = {0.}; 
	// Elastic rate
	
	for (int i = gg; i <= qqb; i++)
	{
	    rate[i] = elasticTable[i].rate * T; 
	    rate[i] += casimir[i] * pow(g*g*T, 2) / (2*wMax*T); 
            // Here we set the cutoff to be Lambda/2, which is different from Derek's note. 
	    rate[i] -= casimir[i] * pow(g*g*T, 2) / Lambda; 	// Add the cutoff Lambda
            rate[i] -= elasticTable[i].ratew / pRest * T * T; 
            rate[i] -= 1. / 96 / M_PI / pRest * pow(g*g*T, 2) * Toverp_c_ln[i]*log(wMax); 
            rate[i] += 1. / 96 / M_PI / pRest * pow(g*g*T, 2) * Toverp_c_ln[i]*log(Lambda/T/2);
            // std::cout << elasticTable[i].rate*T + casimir[i] * pow(g*g*T, 2) / (2*wMax*T) << " " << casimir[i] * pow(g*g, 2)*T/2 << " " << elasticTable[i].ratew*T+1./96/M_PI*pow(g*g, 2)*T*Toverp_c_ln[i]*log(wMax) << " " << 1. / 96 / M_PI * pow(g*g, 2)*T * Toverp_c_ln[i] << "\n"; 
            // std::cout << "large angle elas " << ProcessStrings[i] << " " << rate[i] << "\n";  
	}
        for (int i = GqQg; i <= QbqGg; i++)
        {
            rate[i] = elasticTable[i].rate * T * T / pRest; 
            rate[i] += 1. / 96 / M_PI * pow(g*g, 2) * T * Toverp_c_ln[i]*log(wMax) * T /pRest; 
            rate[i] -= 1. / 96 / M_PI * pow(g*g, 2) * T * Toverp_c_ln[i]*log(Lambda/T/2) * T / pRest;
            // std::cout <<  elasticTable[i].rate * T + 1. / 96 / M_PI * pow(g*g, 2) * T * Toverp_c_ln[i]*log(wMax) << " " << 1. / 96 / M_PI * pow(g*g, 2) * T * Toverp_c_ln[i] << "\n";     
            // std::cout << "conversion elas " << ProcessStrings[i] << " " << rate[i] << "\n"; 
        }
   
	// Conversion rate
	/*for (int i = gqqg; i <= qg_inel_conv; i++)
	{
		process_type iProcess = static_cast<process_type>(i); 
		rate[i] = rate_conv(iProcess, T, pRest);
	}
        */

	// Inelastic rate
	rate[qqg] = rate_inel(pRest, T, rate_qqg_p); 
	rate[gqq] = rate_inel(pRest, T, rate_gqq_p); 
	rate[ggg] = rate_inel(pRest, T, rate_ggg_p); 

	for (int i = gg_split; i <= qqb_split; i++)
	{
		
		rate[i] = splitting_c_lambda[i - gg_split]*pRest/Lambda + splitting_c_p[i - gg_split] + splitting_c_ln[i - gg_split]*log(2*pRest/Lambda); 
		if (i == gg_split || i == qg_split)
			rate[i] *= pow(g, 4)*pow(T, 2)/96/M_PI/pRest; 
		else
			rate[i] *= pow(g, 4)*pow(T, 2)/1536/M_PI/pRest; 
		if (i == qqp_split) rate[i] *= 4.; 
		if (i == gq_split) rate[i] *= 6; 
                // std::cout << "rate elas split " << ProcessStrings[i] << " " << rate[i] << "\n"; 
	}
	if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
	{
		double totalQuarkProb = 0.; 
		/*if (pRest/T > AMYpCut)
			totalQuarkProb += rate[qqg]*dt;*/
		totalQuarkProb += (rate[qg] + rate[qq] + rate[qqp] + rate[qqb] + rate[QgGq] + rate[QbqGg]/* + rate[qggq] + rate[qg_inel_conv]*/ + rate[qg_split] + rate[qq_split] + rate[qqp_split] + rate[qqb_split]) * dt; 
		// warn if total probability exceeds 0.1
		if (totalQuarkProb > 0.5)
      		WARN << " : Total Probability for quark processes exceeds 0.5 (" << totalQuarkProb << "). " << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01)."; 	
      		
  		double accumProb = 0.; 
  		double nextProb = 0.; 
  		double Prob = 0.; 
  		double randProb = ZeroOneDistribution(*GetMt19937Generator()); 
  		
  		if (randProb < totalQuarkProb)
  		{
  			Prob = rate[qg] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qg; 
  			
  			accumProb += Prob; 
  			Prob = rate[qq] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qq;

			accumProb += Prob; 
  			Prob = rate[qqp] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqp; 

			accumProb += Prob; 
  			Prob = rate[qqb] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqb; 

			accumProb += Prob; 
  			Prob = rate[QgGq] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return QgGq; 

			accumProb += Prob; 
  			Prob = rate[QbqGg] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return QbqGg; 
/*  			
  			accumProb += Prob; 
  			Prob = rate[qggq] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qggq; 

			accumProb += Prob; 
  			Prob = rate[qg_inel_conv] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qg_inel_conv; 
  			
  			accumProb += Prob; 
  			Prob = rate[qqg] * dt; 
  			if (pRest/T > AMYpCut && accumProb <= randProb && randProb < (accumProb + Prob)) return qqg; 
*/
			accumProb += Prob; 
			Prob = rate[qg_split] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qg_split; 
  			
  			accumProb += Prob; 
  			Prob = rate[qq_split] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qq_split; 

			accumProb += Prob; 
  			Prob = rate[qqp_split] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqp_split; 

			accumProb += Prob; 
  			Prob = rate[qqb_split] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqb_split; 
  		}
  		else
  			return none; 
  	}
  	else if (Id == 21)
  	{
  		double totalGluonProb = 0.; 
  		/*if (pRest/T > AMYpCut) 
			totalGluonProb += (rate[gqq] + rate[ggg])*dt;*/
		// WARN << totalGluonProb; 
  		totalGluonProb += (rate[gg] + rate[gq] + rate[GqQg] + rate[GgQbq]/* + rate[gqqg] + rate[gq_inel_conv]*/ + rate[gg_split] + rate[gq_split]) * dt; 
		if (totalGluonProb > 0.5)
  			WARN << " : Total Probability for gluon processes exceeds 0.5 (" << totalGluonProb << "). " << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01)."; 
  		
  		double accumProb = 0.; 
  		double nextProb = 0.; 
  		double Prob = 0.; 
  		double randProb = ZeroOneDistribution(*GetMt19937Generator()); 
  		
  		if (randProb < totalGluonProb)
  		{
  			Prob = rate[gg] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return gg; 
  			
  			accumProb += Prob; 
  			Prob = rate[gq] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return gq; 
  			
  			accumProb += Prob; 
  			Prob = rate[GqQg] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return GqQg; 
  			
  			accumProb += Prob; 
  			Prob = rate[GgQbq] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return GgQbq; 
/*  			
  			accumProb += Prob; 
  			Prob = rate[gqqg] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return gqqg; 

			accumProb += Prob; 
  			Prob = rate[gq_inel_conv] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return gq_inel_conv; 

  			accumProb += Prob; 
  			Prob = rate[ggg] * dt; 
  			if (pRest > AMYpCut && accumProb <= randProb && randProb < (accumProb + Prob)) return ggg; 

			accumProb += Prob; 
  			Prob = rate[gqq] * dt; 
  			if (pRest > AMYpCut && accumProb <= randProb && randProb < (accumProb + Prob)) return gqq; 
*/
			accumProb += Prob; 
			Prob = rate[gg_split] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return gg_split; 
  			
  			accumProb += Prob; 
  			Prob = rate[gq_split] * dt; 
  			if (accumProb <= randProb && randProb < (accumProb + Prob)) return gq_split; 
		}
  		else return none; 
  	}
  	return none; 
}

FourVector Tequila::Momentum_Update(double omega, double qperp, double T, FourVector pVec)
{
	double p = pVec.t(); 
	double pp = p - omega; 
	double sintheta = fabs(qperp/pp); 
	double phi = 2.*M_PI*ZeroOneDistribution(*GetMt19937Generator()); 
	fourvec pOut, pIn; 
	FourVector pVecNew; 
	pIn = fourvec{pVec.t(), pVec.x(), pVec.y(), pVec.z()}; 
	// Assume pVec is along z axis. 
	pOut = fourvec{fabs(pp), fabs(pp)*sintheta*cos(phi), fabs(pp)*sintheta*sin(phi), pp*sqrt(1-sintheta*sintheta)}; 
	pOut = pOut.rotate_back(pIn); 
	pVecNew = FourVector(pOut.x(), pOut.y(), pOut.z(), pOut.t()); 
	return pVecNew; 
}

double Tequila::qpara(double E, double T, int id)
{
	double CR; 
	if (id == 21) CR = CA; 
	else if (std::abs(id) == 1 || std::abs(id) == 2 || std::abs(id) == 3) CR = CF; 
	else {WARN << "Strange particle! ID = " << id; CR = CF; }
	double mD = sqrt(std::pow(g*T, 2)*(nc/3. + nf/6.)); 
	// double mD = sqrt(std::pow(g*T, 2)*nc/3.); 
	double Minf = sqrt(pow(mD, 2)/2.); 
	// return 0.; 
	return std::pow(g*Minf, 2)*CR*T/(2.*M_PI)*log(1.+pow(muperp*T/Minf, 2))/2.;
	// 		+ std::pow(g, 4)*CR*CA*std::pow(T, 3)*omega_over_T_cutoff*(2-ln2)/(4.*std::pow(M_PI, 3));
	// return std::pow(g*Minf, 2)*CR*T/(2.*M_PI)*log(muperp*T/Minf);
}

double Tequila::qperp(double E, double T, int id)
{
	double CR; 
	if (id == 21) CR = CA; 
	else if (std::abs(id) == 1 || std::abs(id) == 2 || std::abs(id) == 3) CR = CF; 
	else {WARN << "Strange particle! ID = " << id; CR = CF; }
	double mD = sqrt(std::pow(g*T, 2)*(nc/3. + nf/6.)); 
	// double mD = sqrt(std::pow(g*T, 2)*nc/3.); 
	return std::pow(g*mD, 2) * CR * T / (2.*M_PI) * log(1.+pow(muperp*T/mD, 2))/2.;
	// return 0.; 
}

FourVector Tequila::Langevin_Update(double dt, double T, FourVector pIn, int id)
{
	// imaging rotating to a frame where pIn lies on z-axis
	double E0 = pIn.t();
	double p0 = std::sqrt(E0*E0 - M*M + 1e-9);
	double kt = qperp(E0, T, id) / 2.;
	double kl = qpara(E0, T, id); 
	double drag = kl/(2.*p0*T)+1./(2.*p0*p0)*(kt*2.-kl*2.); //eta_D
		   
    double white_noise_holder[3];
    for (size_t i=0; i<3; ++i) 
		white_noise_holder[i] = Srandom::white_noise(Srandom::gen);

	double Ct = std::sqrt(kt*dt);
	double Cl = std::sqrt(kl*dt);
	
	fourvec pOut; 
    pOut.a[1] = Ct * white_noise_holder[0];
    pOut.a[2] = Ct * white_noise_holder[1];
    pOut.a[3] = p0 * (1. - drag * dt) + Cl * white_noise_holder[2];
    pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
				+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// rotate back to the original frame
	pOut = pOut.rotate_back(fourvec{pIn.t(), pIn.x(), pIn.y(), pIn.z()});
	return FourVector(pOut.x(), pOut.y(), pOut.z(), pOut.t()); 
}

double Tequila::TransverseMomentum_Transfer(double pRest, double omega, double T, process_type process)
{
	pRest /= T; 
	omega /= T; 
	// if (fabs(pRest-omega)*T < pcut || (pRest-omega) < muperp) return 0; 	
	// double qperpMax = 2.*sqrt(kMax*kMax + kMax*omega); 
	double limitMax = std::min(qperpMax, fabs(pRest-omega));
	// limitMax = std::min(limitMax, fabs(omega)); 
	if (limitMax*T < pcut || limitMax < muperp) return 0.;  

	// Rejection method
	const int max_try = 10000; 
	int n_try = 0; 
	double qperp = 0., qperp_test, rate_qperp_test, max_rate = 0., max_rate2 = 0., r; 

	int i_omega_index = (atan(omega)-atan(omegaMin))/(atan(wMax)-atan(omegaMin))*Nw; 
	if (omega < 0) i_omega_index++; 
	iqperp0 = (int)((log(muperp)-log(muperp0))/(log(qperpMax)-log(muperp0))*Nq); 
	for (size_t j = iqperp0; j < Nq; j++)
	{
		double i_rate = elasticTable[process].rate_qperp[i_omega_index][j];
		if (max_rate < i_rate) max_rate = i_rate; 
	}
	max_rate *= 1.5; 
	while (n_try < max_try)
	{
		qperp_test = muperp+ZeroOneDistribution(*GetMt19937Generator())*(limitMax-muperp);
		rate_qperp_test = Interpolator_dGamma_domega_qperp(omega, qperp_test, process); 
		r = max_rate * ZeroOneDistribution(*GetMt19937Generator());
		if (rate_qperp_test > max_rate) WARN << "The sampled rate is larger than the maximum rate we assumed in elastic 2d part!! " << rate_qperp_test << " " << max_rate; 
		if (r < rate_qperp_test)
		{
			qperp = qperp_test; 
			break; 
		}
		n_try++; 
		if (n_try == max_try) WARN << "cannot find a proper qperp in elastic part!!!"; 
	}
	// The qperp in the interpolator is qperp/T
	return qperp; 
}

double Tequila::Energy_Transfer(double pRest, double T, process_type process)
{
	// Rejection method
	const int max_try = 10000; 
	int n_try = 0; 
	double omega = 0., omega_test=0., rate_omega_test, max_rate = 0., max_omega, r; 

	for (size_t i = 0; i < Nw; i++)
	{
		double i_rate = exp(elasticTable[process].y[i])/(elasticTable[process].x[i]*elasticTable[process].x[i]+1.); 
		// if (max_rate < i_rate && (elasticTable[process].x[i] > muperp || elasticTable[process].x[i] < -muperp)) max_rate = i_rate; 
		if (max_rate < i_rate) max_rate = i_rate; 
	}
	max_rate *= 1.5;
	while (n_try < max_try)
	{
		// while (omega_test < muperp && omega_test > -muperp)
		omega_test = omegaMin+ZeroOneDistribution(*GetMt19937Generator())*(Lambda/(2.*T)-omegaMin);
		// omega_test = omegaMin+ZeroOneDistribution(*GetMt19937Generator())*(kMax-2*muperp-omegaMin);
		// if (omega_test > -muperp) omega_test += 2*muperp;
		if (process == gg || process == gq || process == qg || process == qq || process == qqp || process == qqb) 
          	    rate_omega_test = Interpolator_dGamma_domega(omega_test, process)*(1.-omega_test*T/pRest); 
                else
                    rate_omega_test = Interpolator_dGamma_domega(omega_test, process); 
                // std::cout << 1.-omega_test*T/pRest << "\n"; 
                // std::cout << omega_test << " " << rate_omega_test << " " << max_rate << "\n"; 
		if (rate_omega_test > max_rate) WARN << "The sampled rate is larger than the maximum rate we assumed in elastic 1d part!!"; 
		r = max_rate * ZeroOneDistribution(*GetMt19937Generator());
		if (r < rate_omega_test)
		{
			omega = omega_test; 
			break; 
		}
		n_try++; 
                // std::cout << n_try << " " << omega << "\n"; 
		if (n_try == max_try) WARN << "cannot find a proper omega in elastic part!!!"; 
	}
	// The omega returned by the interpolator is omega/T
	// std::cout << omega << "\n"; 
	return omega; 
}

double Tequila::TransverseMomentum_Transfer_Split(double pRest, double omega, double T, process_type process)
{
	// if (fabs(pRest-omega)*T < pcut || (pRest-omega) < muperp) return 0; 	
	// double qperpMax = 2.*sqrt(kMax*kMax + kMax*omega); 
	// double limitMax = std::min(qperpMax, fabs(pRest/T-omega));
	double limitMax = fabs(pRest/T-omega);
	if (limitMax*T < pcut || limitMax < muperp) return 0.;  
/*	// Metropolis method
	// Randomly select initial values of qperp and its rate. 
	double qperp = ZeroOneDistribution(*GetMt19937Generator())*limitMax;
	double rate_qperp = splittingRateQperp(pRest, qperp, omega, T, process); 
	double qperpNew, rate_qperpNew, ratio; 
	
	for (int step = 0; step < Nsteps; step++)
	{
		qperpNew = muperp+ZeroOneDistribution(*GetMt19937Generator())*(limitMax-muperp);
		rate_qperpNew = splittingRateQperp(pRest, qperpNew, omega, T, process); 
		ratio = rate_qperpNew / rate_qperp; 
		if (ZeroOneDistribution(*GetMt19937Generator()) < ratio)
		{
			qperp = qperpNew; 
			rate_qperp = rate_qperpNew; 
		}
	}*/

	// Rejection method
	const int max_try = 1000000; 
	int n_try = 0; 
	double qperp = 0., qperp_test, rate_qperp_test, max_rate = 0., max_rate2 = 0., r; 

	int i_omega_index = (atan(omega)-atan(omegaMin))/(atan(wMax)-atan(omegaMin))*Nw; 
	if (omega < 0) i_omega_index++; 
	iqperp0 = (int)((log(muperp)-log(muperp0))/(log(qperpMax)-log(muperp0))*Nq); 
	for (size_t j = iqperp0; j < Nq; j++)
	{
		double i_rate = elasticTable[process].rate_qperp[i_omega_index][j];
		if (max_rate < i_rate) max_rate = i_rate; 
	}
	max_rate *= 1.5; 
	while (n_try < max_try)
	{
		qperp_test = muperp+ZeroOneDistribution(*GetMt19937Generator())*(limitMax-muperp);
		rate_qperp_test = Interpolator_dGamma_domega_qperp(omega, qperp_test, process); 
		r = max_rate * ZeroOneDistribution(*GetMt19937Generator());
		if (rate_qperp_test > max_rate) WARN << "The sampled rate is larger than the maximum rate we assumed in elastic 2d part!! " << rate_qperp_test << " " << max_rate; 
		if (r < rate_qperp_test)
		{
			qperp = qperp_test; 
			break; 
		}
		n_try++; 
		if (n_try == max_try) WARN << "cannot find a proper qperp in elastic part!!!"; 
	}
	// The qperp returned by the interpolator is qperp/T
	return qperp; 
}

double Tequila::thermalDistributionSampling(process_type process)
{
	double k, k_new, rate_k, rate_k_new, ratio; 
	// int n_try = 0, max_try = 10000; 
	int Nsteps = 500; 
	double k_max = 16., k_min = 1e-3; 
	if (process == gg_split || process == qg_split)
	{
		k = k_min + ZeroOneDistribution(*GetMt19937Generator())*(k_max-k_min); 
		rate_k = nB(k); 
		for (int step = 0; step < Nsteps; step++)
		{
			k_new = k_min + ZeroOneDistribution(*GetMt19937Generator())*(k_max-k_min); 
			rate_k_new = nB(k_new); 
			ratio = rate_k_new / rate_k; 
			double rn = ZeroOneDistribution(*GetMt19937Generator()); 
			if (rn < ratio)
			{
				k = k_new; 
				rate_k = rate_k_new; 
			}
		}
/*		
		max_rate = nB(k_min); 
		while (n_try < max_try)
		{
			k_test = k_min + ZeroOneDistribution(*GetMt19937Generator())*(k_max-k_min); 
			r = max_rate * ZeroOneDistribution(*GetMt19937Generator());		
			if (r < nB(k_test))
			{
				k = k_test; 
				break; 
			}
			n_try++; 
			if (n_try == max_try) WARN << "cannot find a proper omega in elastic part!!!"; 
		}*/
	}
	else
	{
		k = k_min + ZeroOneDistribution(*GetMt19937Generator())*(k_max-k_min); 
		rate_k = nF(k); 
		for (int step = 0; step < Nsteps; step++)
		{
			k_new = k_min + ZeroOneDistribution(*GetMt19937Generator())*(k_max-k_min); 
			rate_k_new = nF(k_new); 
			ratio = rate_k_new / rate_k; 
			double rn = ZeroOneDistribution(*GetMt19937Generator()); 
			if (rn < ratio)
			{
				k = k_new; 
				rate_k = rate_k_new; 
			}
		}
/*
		max_rate = nF(k_min); 
		while (n_try < max_try)
		{
			k_test = k_min + ZeroOneDistribution(*GetMt19937Generator())*(k_max-k_min); 
			r = max_rate * ZeroOneDistribution(*GetMt19937Generator());		
			if (r < nF(k_test))
			{
				k = k_test; 
				break; 
			}
			n_try++; 
			if (n_try == max_try) WARN << "cannot find a proper omega in elastic part!!!"; 
		}*/
	}
	// WARN << "thermal distribution: " << k; 
	return k; 
}

double Tequila::xSampling(double pRest, double kRest, double T, process_type process)
{
	// double x = (Lambda/(2.*T)-kRest)/(pRest-kRest) + ZeroOneDistribution(*GetMt19937Generator()) * (1.-2.*(Lambda/(2.*T)-kRest)/(pRest-kRest)); 
	double x = Lambda/(2.*pRest) + ZeroOneDistribution(*GetMt19937Generator()) * (1.-Lambda/pRest); 
	// double x = Lambda/(2.*pRest) + ZeroOneDistribution(*GetMt19937Generator()) * (0.5-Lambda/2/pRest);
	double rate_x = splittingRateOmega(pRest, x, T, process); 
	double x_new, rate_x_new, ratio; 
	
	for (int step = 0; step < Nsteps; step++)
	{
		// x_new = (Lambda/(2.*T)-kRest)/(pRest-kRest) + ZeroOneDistribution(*GetMt19937Generator()) * (1.-2.*(Lambda/(2.*T)-kRest)/(pRest-kRest)); 
		x_new = Lambda/(2.*pRest) + ZeroOneDistribution(*GetMt19937Generator()) * (1.-Lambda/pRest); 
		rate_x_new = splittingRateOmega(pRest, x_new, T, process); 
		ratio = rate_x_new / rate_x; 
		double rn = ZeroOneDistribution(*GetMt19937Generator()); 
		if (rn < ratio)
		{
			x = x_new; 
			rate_x = rate_x_new; 
		}
	}
	return x; 
}

double Tequila::Energy_Transfer_Split(double pRest, double T, process_type process)
{
	// Metropolis method
	double omega = Lambda/(2.*T)+ZeroOneDistribution(*GetMt19937Generator())*(pRest/T-Lambda/(2.*T)-Lambda/(2.*T));
	double rate_omega = splittingRateOmega(pRest, omega, T, process); 
	double omegaNew, rate_omegaNew, ratio; 
	
	for (int step = 0; step < Nsteps; step++)
	{
		omegaNew = Lambda/(2.*T)+ZeroOneDistribution(*GetMt19937Generator())*(pRest/T-Lambda/(2.*T)-Lambda/(2.*T));
		rate_omegaNew = splittingRateOmega(pRest, omegaNew, T, process); 
		ratio = rate_omegaNew / rate_omega; 
		double rn = ZeroOneDistribution(*GetMt19937Generator()); 
		if (rn < ratio)
		{
			omega = omegaNew; 
			rate_omega = rate_omegaNew; 
		}
	}
	/*const int max_try = 10000; 
	int n_try = 0; 
	double omega = 0., omega_test=0., rate_omega_test, max_rate = 0., max_omega, r; 

	max_rate = splittingRateOmega(pRest, Lambda/(2.*T), T, process); 
	while (n_try < max_try)
	{
		// while (omega_test < muperp && omega_test > -muperp)
		omega_test = Lambda/(2.*T)+ZeroOneDistribution(*GetMt19937Generator())*(pRest/T-Lambda/(2.*T)-Lambda/(2.*T));
		// omega_test = omegaMin+ZeroOneDistribution(*GetMt19937Generator())*(kMax-2*muperp-omegaMin);
		// if (omega_test > -muperp) omega_test += 2*muperp; 
		rate_omega_test = splittingRateOmega(pRest, omega_test, T, process); 
		// INFO << n_try << " sampled omega " << omega_test << " " << rate_omega_test<< " " << max_rate; 
		if (rate_omega_test > max_rate) WARN << "The sampled rate is larger than the maximum rate we assumed in elastic 1d part!! " << "process "<< process << " omega test " << omega_test << " omega min " << Lambda/(2.*T); 
		r = max_rate * ZeroOneDistribution(*GetMt19937Generator());
		// INFO << BOLDYELLOW << "new omega " << omega_test << " rate " << rate_omega_test << " r " << r << " ratio " << rate_omega_test/max_rate; 		
		if (r < rate_omega_test)
		{
			omega = omega_test; 
			// INFO << BOLDYELLOW << "the number of trials is " << n_try; 
			break; 
		}
		n_try++; 
		if (n_try == max_try) WARN << "cannot find a proper omega in elastic part!!!"; 
	}*/
	// The omega returned by the interpolator is omega/T
	return omega; 
}

bool Tequila::is_empty(std::ifstream& pFile)
{
	return pFile.peek() == std::ifstream::traits_type::eof();
}

bool Tequila::is_exist (const std::string& name)
{
	struct stat buffer;   
	return (stat (name.c_str(), &buffer) == 0); 
}
/*
double Tequila::rate_conv(process_type process, double T, double pRest)
{
	double m_inf2 = pow(g*T, 2) * CF / 4.; 
	double m_D2 = pow(g*T, 2) * (nc/3. + nf/6.); 
	if (process == gqqg || process == qggq)
	{
		double rate_base = pow(g, 2) * m_inf2 * CF / (8. * M_PI * pRest) * log(1.+pow(muperp*T, 2)/m_inf2)/2.; 
		if (process == gqqg) return dF/(double)dA*rate_base; 
		else if (process == qggq) return rate_base; 
		else {WARN << "Invalid conversion process " << std::to_string(process); return 0.; }
	}
	else
	{
		double rate_base = pow(g, 4) * pow(T, 2) * CF * omega_over_T_cutoff * (2-log(2)) * m_D2 / (32. * pow(M_PI, 3) * m_inf2); 
		// INFO << BOLDYELLOW << "p " << pRest << " rate " << rate_base; 
		if (process == gq_inel_conv) return rate_base/2.; 
		else if (process == qg_inel_conv) return rate_base*CF; 
		else {WARN << "Invalid conversion process " << std::to_string(process); return 0.; }
	}
	// else {WARN << "Invalid conversion process " << std::to_string(process); return 0.; }
}
*/
double Tequila::splitting_gg1(double x)
{
    return 16.*CA*CA*(3.+(1.-x)/pow(x, 2)+x/pow(1.-x, 2)-x*(1.-x));
}

double Tequila::splitting_gg2(double x)
{
    return 8.*CA*(1.+pow(x, 4)+pow(1.-x, 4))/(pow(x, 2)*pow(1.-x, 2))*(CA/2*pow(x, 2)+CA/2*(1.+pow(1.-x, 2)));
}

double Tequila::splittingF(double x, process_type process)
{
	double F = 0.; 
	switch(process)
	{
		case gg_split: F = 16.*CA*CA*(3.+(1.-x)/pow(x, 2)+x/pow(1.-x, 2)-x*(1.-x)); 
		// case gg_split: F = 16.*CA*CA*(3+(1.-x)/pow(x, 2)); 
		//	 break; 
		// case gg_split: F = 8.*CA*(1.+pow(x, 4)+pow(1.-x, 4))/(pow(x, 2)*pow(1.-x, 2))*(CA/2*pow(x, 2)+CA/2*(1.+pow(1.-x, 2))); 
		         break; 
		case gq_split: F = 4.*(1.+pow(1.-x, 2))/(pow(x, 2)*(1.-x))*(CF*pow(x, 2)+CA*(1.-x))*6./2; 
		 	 break; 
		// case gq_split: F = 4.*CA*(1.+pow(x, 2))/(x*pow(1.-x, 2))*((CF-CA/2)*pow(1.-x, 2)+CA/2*(1.+pow(x, 2)))*6./2; 
		// 	 break; 
		case qg_split: F = 8.*CF*(1.+pow(1.-x, 2))/(pow(x, 2)*(1.-x))*(CF*pow(x, 2)+CA*(1.-x));
		// 	WARN << "x is " << x;  
			 break; 
		case qq_split: F = (4.*CF*((1.+pow(1.-x, 2))/pow(x, 2)+(1.+pow(x, 2))/pow(1.-x, 2))+16.*CF*(CF-CA/2)/x/(1.-x))/2/2; 
			 break; 
		case qqp_split: F = 4.*CF*(1.+pow(1.-x, 2))/pow(x, 2)*4./2; 
			 break; 
		case qqb_split: F = 4.*CF*((1.+pow(1.-x, 2))/pow(x, 2)+pow(x, 2)+pow(1.-x, 2))-16.*CF*(CF-CA/2)*pow(1.-x, 2)/x/2; 
			 break; 
		default: WARN << "The process determination is wrong! The process is " << process; 
			 break; 
	}
	return F; 
}

double Tequila::splittingRateOmega(double pRest, double x, double T, process_type process)
{
	return pow(g, 4)*pow(T, 2)/(96.*8.*M_PI*pow(pRest, 2)) * splittingF(x, process); 
}

/*double Tequila::splittingRateQperp(double pRest, double qperp, double omega, double T, process_type process)
{

	double q = sqrt(pow(qperp, 2)+pow(omega, 2)); 
	double limit = (q-omega)/2; 
	double rate = 0; 
	if (process == gg_split || process == qg_split)
		rate = qperp/q * (limit-log(1-exp(limit))) * splittingF(pRest, omega, T, process); 
	else if (process == gq_split || process == qq_split || process == qqp_split || process == qqb_split)
		rate = qperp/q * (log(exp(limit)+1)-limit) * splittingF(pRest, omega, T, process); 
	return rate; 
}*/

void Tequila::LoadElasticTables()
{
	IntTabulator inttabulator; 
	gsl_interp2d *interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, Nw+1, Nq+1); 

	// Iterate over all the elastic processes	
	for (int iProcess = gg; iProcess <= qqb; iProcess++)
	{
		INFO << BOLDYELLOW << "process is " << inttabulator.GetProcessString(iProcess); 
		process_type process = static_cast<process_type>(iProcess); 
		tables iTables; 
		if (!is_exist((path_to_tables+"elastic_rate_table"+inttabulator.GetProcessString(process)+".dat").c_str())) inttabulator.Tabulator_dGamma_domega_qperp(path_to_tables, process); 
		std::ifstream table2d_in((path_to_tables+"elastic_rate_table"+inttabulator.GetProcessString(process)+".dat").c_str()); 
		if (is_empty(table2d_in)) inttabulator.Tabulator_dGamma_domega_qperp(path_to_tables, process); 
		double z; 
    	for (size_t iomega = 0; iomega <= Nw; iomega++)
    	{
    		iTables.xa[iomega] = tan(((double)iomega)*(atan(wMax)-atan(omegaMin))/Nw+atan(omegaMin)); 
			// double qperpMax = 2.*sqrt(kMax*kMax + kMax*iTables.xa[iomega]); 
		for (size_t iqperp = 0; iqperp <= Nq; iqperp++)
    		{
    			iTables.ya[iqperp] = exp(((double)iqperp)*(log(qperpMax)-log(muperp0))/Nq+log(muperp0)); 
			table2d_in >> z; 
			iTables.rate_qperp[iomega][iqperp] = z; 
			gsl_interp2d_set(interp, &za[iProcess*(Nw+1)*(Nq+1)], iomega, iqperp, log(z)); 
    		}
    	}
	table2d_in.close(); 

	for (size_t iomega = 0; iomega <= Nw; iomega++)
    	{
		iTables.x[iomega] = iTables.xa[iomega]; 
		iTables.y[iomega] = 0.;
		// double qperpMax = 2.*sqrt(kMax*kMax + kMax*iTables.x[iomega]); 
		iqperp0 = (int)((log(muperp)-log(muperp0))/(log(qperpMax)-log(muperp0))*Nq); 
		for (size_t iqperp = iqperp0; iqperp < Nq; iqperp++)
    		{
			double dy = exp(((double)(iqperp+1))*(log(qperpMax)-log(muperp0))/Nq+log(muperp0)) - exp(((double)iqperp)*(log(qperpMax)-log(muperp0))/Nq+log(muperp0)); 
		        iTables.y[iomega] += dy * (iTables.rate_qperp[iomega][iqperp] + iTables.rate_qperp[iomega][iqperp+1]) / 2; 
		}
	}

	for (size_t iomega = 0; iomega < Nw; iomega++)
    	{
    		double dx = iTables.x[iomega+1] - iTables.x[iomega]; 
			iTables.rate += (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2;
			iTables.ratew += (iTables.y[iomega+1] + iTables.y[iomega]) * dx * (iTables.x[iomega]+iTables.x[iomega+1]) / 2 / 2;

			iTables.y[iomega] = log(iTables.y[iomega]*(iTables.x[iomega]*iTables.x[iomega]+1)); 
		}
		iTables.y[Nw] = log(iTables.y[Nw]*(iTables.x[Nw]*iTables.x[Nw]+1)); 
		iTables.rate *= pow(g, 4);
                std::cout << "total rate is " << iTables.rate << "\n";  
		iTables.ratew *= pow(g, 4); 
		INFO << BOLDGREEN << "muperp is " << muperp << " rate is " << iTables.rate; 
		elasticTable.push_back(iTables); 
	}
	for (int iProcess = GqQg; iProcess <= QbqGg; iProcess++)
	{
		INFO << BOLDYELLOW << "process is " << inttabulator.GetProcessString(iProcess); 
		process_type process = static_cast<process_type>(iProcess); 
		tables iTables; 
                std::cout << "i process is " << iProcess << "\n"; 
		if (!is_exist((path_to_tables+"elastic_rate_table"+inttabulator.GetProcessString(process)+"for_wdGdw.dat").c_str())) inttabulator.Tabulator_omega_dGamma_domega_qperp(path_to_tables, process); 
		std::ifstream table2d_in_2((path_to_tables+"elastic_rate_table"+inttabulator.GetProcessString(process)+"for_wdGdw.dat").c_str()); 
		if (is_empty(table2d_in_2)) inttabulator.Tabulator_omega_dGamma_domega_qperp(path_to_tables, process); 
		double z; 
    	for (size_t iomega = 0; iomega <= Nw; iomega++)
    	{
    		iTables.xa[iomega] = ((double)iomega)*(wMax-omegaMin)/Nw+omegaMin; 
			// double qperpMax = 2.*sqrt(kMax*kMax + kMax*iTables.xa[iomega]); 
		for (size_t iqperp = 0; iqperp <= Nq; iqperp++)
    		{
    			iTables.ya[iqperp] = ((double)iqperp)*(qperpMax-muperp0)/Nq+muperp0; 
			table2d_in_2 >> z; 
			iTables.rate_qperp[iomega][iqperp] = z; 
			gsl_interp2d_set(interp, &za[iProcess*(Nw+1)*(Nq+1)], iomega, iqperp, log(z)); 
    		}
    	}
	table2d_in_2.close(); 

	for (size_t iomega = 0; iomega <= Nw; iomega++)
    	{
		iTables.x[iomega] = iTables.xa[iomega]; 
		iTables.y[iomega] = 0.;
		// double qperpMax = 2.*sqrt(kMax*kMax + kMax*iTables.x[iomega]); 
		iqperp0 = (int)((muperp-muperp0)/(qperpMax-muperp0)*Nq); 
		for (size_t iqperp = iqperp0; iqperp < Nq; iqperp++)
    		{
			double dy = (qperpMax-muperp0)/Nq; 
			iTables.y[iomega] += dy * (iTables.rate_qperp[iomega][iqperp] + iTables.rate_qperp[iomega][iqperp+1]) / 2; 
		}
	}

	for (size_t iomega = 0; iomega < Nw; iomega++)
    	{
    		double dx = iTables.x[iomega+1] - iTables.x[iomega]; 
			iTables.rate += (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2;
			iTables.ratew += (iTables.y[iomega+1] + iTables.y[iomega]) * dx * (iTables.x[iomega]+iTables.x[iomega+1]) / 2 / 2;
		}
		iTables.rate *= pow(g, 4);
                std::cout << "total rate is " << iTables.rate << "\n";  
		iTables.ratew *= pow(g, 4); 
		INFO << BOLDGREEN << "muperp is " << muperp << " rate is " << iTables.rate; 
		elasticTable.push_back(iTables); 
	}
	gsl_interp2d_free(interp); 
}

double Tequila::Interpolator_dGamma_domega(double omega, process_type process)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp *interp = gsl_interp_alloc (gsl_interp_linear, Nw+1);

	gsl_interp_init (interp, elasticTable[process].x, elasticTable[process].y, Nw+1); 
	double result; 
	result = gsl_interp_eval (interp, elasticTable[process].x, elasticTable[process].y, omega, acc);  
	gsl_interp_free (interp); 
	gsl_interp_accel_free (acc); 
	return exp(result)/(omega*omega+1.); 
}

double Tequila::Interpolator_dGamma_domega_qperp(double omega, double qperp, process_type process)
{	
	gsl_interp2d *interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, Nw+1, Nq+1);
  	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  	gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  	// gsl_interp2d_init(interp, elasticTable[process].xa, elasticTable[process].ya, elasticTable[process].za, Nw+1, Nq+1);
	gsl_interp2d_init(interp, elasticTable[process].xa, elasticTable[process].ya, &za[process*(Nw+1)*(Nq+1)], Nw+1, Nq+1);
	// if (omega == elasticTable[iProcess].xa[iqperp] && qperp == elasticTable[iProcess].ya[iqperp]) return z; 
    double result; 
    // result = gsl_interp2d_eval(interp, elasticTable[process].xa, elasticTable[process].ya, elasticTable[process].za, omega, qperp, xacc, yacc); 
	result = gsl_interp2d_eval(interp, elasticTable[process].xa, elasticTable[process].ya, &za[process*(Nw+1)*(Nq+1)], omega, qperp, xacc, yacc); 
    gsl_interp2d_free(interp);
  	gsl_interp_accel_free(xacc);
  	gsl_interp_accel_free(yacc);
    return exp(result); 
}

void Tequila::allocate_memory_for_radiative_rate_table() {
	// Allocate memory for differential and integrated rate
	rate_ggg_p = new double [nb_points_in_p];
	rate_gqq_p = new double [nb_points_in_p];
	rate_qqg_p = new double [nb_points_in_p];
	//maximum_differential_rate = new double [nb_points_in_p];
	differential_rate_ggg_p_omega_qperp = new double ** [nb_points_in_p];
	differential_rate_gqq_p_omega_qperp = new double ** [nb_points_in_p];
	differential_rate_qqg_p_omega_qperp = new double ** [nb_points_in_p];
	for(int ip=0;ip<nb_points_in_p; ip++) {
		differential_rate_ggg_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		differential_rate_gqq_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		differential_rate_qqg_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		for(int iomega=0;iomega<nb_points_in_omega; iomega++) {
			 differential_rate_ggg_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
			 differential_rate_gqq_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
			 differential_rate_qqg_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
		}
	}	
	/*rate_ggg_p = new double [nb_points_in_p];
	rate_gqq_p = new double [nb_points_in_p];
	rate_qqg_p = new double [nb_points_in_p];
	//maximum_differential_rate = new double [nb_points_in_p];
	differential_rate_ggg_p_omega_qperp = new double ** [nb_points_in_p];
	differential_rate_gqq_p_omega_qperp = new double ** [nb_points_in_p];
	differential_rate_qqg_p_omega_qperp = new double ** [nb_points_in_p];
	for(int ip=0;ip<nb_points_in_p; ip++) {
		differential_rate_ggg_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		differential_rate_gqq_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		differential_rate_qqg_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		for(int iomega=0;iomega<nb_points_in_omega; iomega++) {
			 differential_rate_ggg_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
			 differential_rate_gqq_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
			 differential_rate_qqg_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
		}
	}*/

}

//Load rate differential in incoming parton energy p and radiated parton energy omega (not differential in transverse momentum q)
//The table contains dGamma/(dx domega)= gs^4*use_table(), where use_table() is the function from Moore's code that tabulates the collinear rate
void Tequila::load_differential_rate(const double alpha_s, const double alpha_EM, const int Nf, const std::string location_of_collinear_rates) {

	// Check if the tabulated rate is already available for these parameters
	// The format is "dGamma_dp_domega_NfX_alphasYYY" where "X" is the number of flavours and "YYY" is the value of alpha_s
	std::stringstream filename;
	filename << location_of_collinear_rates << "dGamma_dp_domega_Nf" << Nf << "_alphaEM" << alpha_EM << "_alphaS" << alpha_s;
	//Open file
	std::ifstream rate_file;
	rate_file.open(filename.str().c_str(),std::fstream::in);

/*	if (!rate_file.is_open()) {
		std::cout << "Can't open rate table file, aborting...\n";
		exit(1);
	}*/

	//Struct [from Moore's code] to store rates
	dGammas Moore_rate_arrays;
	Gamma_info Moore_rate_info;
	// Save file with proper name
	sprintf(Moore_rate_info.in_fname, "%s", filename.str().c_str());

	// If not, run Guy's program given parameters
	if (!rate_file.is_open()) {

		//
		std::cout << "Pre-tabulated collinear rate not available. Tabulating it now...\n";
	
		//Pre-initialized in data: Nc, Nf, alpha_s, alpha
		Moore_rate_info.Nf=Nf;
		Moore_rate_info.Nc=3;
		Moore_rate_info.alpha_s=alpha_s;
		Moore_rate_info.alpha=alpha_EM;
		build_table(&Moore_rate_info , &Moore_rate_arrays);
		
		std::cout << "... and writing it into file \"" << filename.str().c_str() << "\" for future use.\n";
		write_table(&Moore_rate_info , &Moore_rate_arrays);

	}
	else {
		std::cout << "Collinear rate available for value of alpha_s, alpha_EM and Nf. Reading from file.\n";

		rate_file.close();
		// Read rate file in temporary array differential in p and omega using Moore's functions
		read_table(&Moore_rate_info,&Moore_rate_arrays);
		std::cout << "Collinear rates read.\n";
	}
	const double gs4=alpha_s*alpha_s*(16*M_PI*M_PI);


	// Populate member array "differential_rate_p_omega_qperp[]" with p and omega rate from Moore's code and Gaussian distribution for qperp
	// Some useful parameters for the (temporary) qperp distrib
	// The (temporary) Gaussian envelope used in the q_perp direction (since the q_perp distribution is not currently known)
	const double envelope_sigma_sqr= 4.*M_PI*alpha_s/3.*(3 + Nf*0.5); //Debye mass
	for(int ip=0;ip<nb_points_in_p; ip++) {
		const double p_over_T_val=get_p_over_T_from_index(ip);
		for(int iomega=0;iomega<nb_points_in_omega; iomega++) {
			const double omega_over_T_val=get_omega_over_T_from_index(ip,iomega);
			//Object "Moore_rate_arrays" does not actually contains the rate, but rather the rate stripped of various factors 
			//Moore's function "use_table()" returns the actual rate, after multiplications by the proper factors
			//Using "use_table()" to do this is kind-of sub-optimal, but speed shouldn't be too much of an issue here
			//Need to multiply by g_s^4 since the factor is stripped from the rate in Moore's program
			const double tmp_rate_ggg=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma_ggg , 2 );
			const double tmp_rate_gqq=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma_gqq , 1 );
			const double tmp_rate_qqg=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma, 0 );

			for(int iq=0; iq<nb_points_in_qperp;iq++) {	
				// q_perp envelope function: Gaussian of width sigma \propto m_D?
				// Integrate[ (2 Pi r) Exp[-r^2/sigma^2]/Sqrt[Pi sigma^2]^2, {r, 0, Infinity}]=1
				//jacobian between E dG/dq^3 and dG/(dqperp domega)???
				const double qperp_over_T=get_qperp_over_T_from_index(iq);
				//const double envelope=exp(-qperp_over_T*qperp_over_T/envelope_sigma_sqr)/(M_PI*envelope_sigma_sqr);
				// Hack to make the integral over q_perp unity
				//const double jacobian=2*M_PI*qperp_over_T;
				const double envelope=1./(qperp_over_T_max()); 
				const double rate_ggg=tmp_rate_ggg*envelope;
				const double rate_gqq=tmp_rate_gqq*envelope;
				const double rate_qqg=tmp_rate_qqg*envelope;
				//if (ip ==0 && iomega == 0) std::cout << "qperp_over_T=" << qperp_over_T << " & rate=" << rate << "\n";
				//assign rate
				differential_rate_ggg_p_omega_qperp[ip][iomega][iq]=rate_ggg;
				differential_rate_gqq_p_omega_qperp[ip][iomega][iq]=rate_gqq;
				differential_rate_qqg_p_omega_qperp[ip][iomega][iq]=rate_qqg;
				//differential_rate_p_omega_qperp[ip][iomega][iq]=envelope;
			}
		}
	}

}

//Find position in array corresponding to value of p_over_T and omega_over_T
//Matches Moore's code (it must!)
double Tequila::get_index_from_omega_over_T(const double p_over_T, const double omega_over_T) {

	double b;
	
	if ( omega_over_T < 2 )
	{
		if ( omega_over_T < -1 )
		{
			if ( omega_over_T < -2 )
				b = 60 + 5*omega_over_T;
			else
				b = 70+10*omega_over_T;
		}
		else
		{
			if ( omega_over_T < 1 )
				b = 80 + 20*omega_over_T;
			else
				b = 90 + 10*omega_over_T;
		}
	}
	else if ( omega_over_T < p_over_T-2 )
	{ /* This is that tricky middle ground. */
		b = 190 - 10*log ( 1.000670700260932956l / 
				( 0.0003353501304664781l + (omega_over_T-2) / (p_over_T-4) ) - 1 );
	}
	else
	{
		if ( omega_over_T < p_over_T+1 )
		{
			if ( omega_over_T < p_over_T-1 )
				b = 290 + 10*(omega_over_T-p_over_T);
			else
				b = 300 + 20*(omega_over_T-p_over_T);
		}
		else
		{
			if ( omega_over_T < p_over_T+2 )
				b = 310 + 10*(omega_over_T-p_over_T);
			else
				b = 320 + 5*(omega_over_T-p_over_T);
		}
	}

	return b;

}


//Find the value of p_over_T and omega_over_T for a given position in array
//Matches Moore's code (it must!)
double Tequila::get_omega_over_T_from_index(const int index_p_over_T, const int index_omega_over_T) {

	//Need p_over_T to get omega_over_T
	double p_over_T=get_p_over_T_from_index(index_p_over_T);
	double omega_over_T;

	if ( index_omega_over_T < 50 )        /* spaced by 0.2  from -12 to -2 */
		omega_over_T = -12 + index_omega_over_T * 0.2;
	else if ( index_omega_over_T < 60 )   /* spaced by 0.1  from -2  to -1 */
		omega_over_T = -2 + (index_omega_over_T-50) * 0.1;
	else if ( index_omega_over_T < 100 )  /* spaced by 0.05 from -1  to +1 */
		omega_over_T = -1 + (index_omega_over_T-60) * 0.05;
	else if ( index_omega_over_T < 110 )  /* spaced by 0.1  from +1  to +2 */
		omega_over_T = 1 + (index_omega_over_T-100) * 0.1;
	else if ( index_omega_over_T < 270 )  /* spaced complicated, +2 to p_over_T-2 */
	{
		omega_over_T = 0.1 * (index_omega_over_T-190);
		omega_over_T = 2 + (p_over_T-4) * ( -0.0003353501304664781l
				+ 1.000670700260932956l / (1+exp(-omega_over_T)) );
	}
	else if ( index_omega_over_T < 280 )  /* spaced by 0.1  from p_over_T-2 to p_over_T-1 */
		omega_over_T = p_over_T - 2 + 0.1 * (index_omega_over_T-270);
	else if ( index_omega_over_T < 320 )  /* spaced by 0.05 from p_over_T-1 to p_over_T+1 */
		omega_over_T = p_over_T + 0.05 * (index_omega_over_T - 300);
	else if ( index_omega_over_T < 330 )  /* spaced by 0.1  from p_over_T+1 to p_over_T+2 */
		omega_over_T = p_over_T + 0.1 * (index_omega_over_T - 310);
	else                   /* spaced by 0.2  from p_over_T+2 to p_over_T+12 */
		omega_over_T = p_over_T + 0.2 * (index_omega_over_T - 320);

	return omega_over_T;

}

//Evaluated integrated rate, with cut-off "omega_over_T_cut" on omega, from differential rate stored in member object "differential_rate_p_omega_qperp[][]"
//Also save the maximum value of the rate
void Tequila::evaluate_integrated_rate(double omega_over_T_cut, double *** differential_rate, double * integrated_rate) {

	//current omega/T integration not very good if omega_over_T_cut<0.1
	if (omega_over_T_cut<0.1) std::cout << "Warning: omega/T integration is not very good for omega/T cut-off smaller than 0.1\n";
	//loop over all values of "p"
	for(int ip=0;ip<nb_points_in_p; ip++) {

		const double p_over_T=get_p_over_T_from_index(ip);

		// integrate omega/T from -infinity to omega_over_T_cut, and omega_over_T_cut to p/2T
		// the discretization is not uniform in omega/T
		// let's just use the GSL interpolation routine to integrate, since performance is not an issue here
		double integral_omega=0.0;

		// Maximum point in omega/T necessary
		int pOver2T_cut_pos_position_int=ceil(get_index_from_omega_over_T(p_over_T,p_over_T/2.0))+2;

		//Arrays to store the omega/T rates
		double * omega_rate_array = new double [pOver2T_cut_pos_position_int];
		double * omega_position_array = new double [pOver2T_cut_pos_position_int];
		// WARN << "before the loop "; 
		// Fill an array with the qperp-integrated values for each omega points
		for(int iomega=0;iomega<pOver2T_cut_pos_position_int; iomega++) {

			omega_position_array[iomega]=get_omega_over_T_from_index(ip,iomega);

			// integrate in q_perp
			double integral_qperp=0.0;

			// track the position of various timesteps to compute the weight of each step properly with a minimum number of call to function "get_qperp_over_T_from_index()"
			// trapeze integration not very good for radial coordinate. might require improvements.
			double prev_qperp=0.0;
			double curr_qperp=get_qperp_over_T_from_index(0);
			double next_qperp;
			for(int iq=0; iq<nb_points_in_qperp;iq++) {	
				next_qperp=get_qperp_over_T_from_index(iq+1);
				// get weight
				const double weight_qperp=non_uniform_trapeze_weights(iq,nb_points_in_qperp,prev_qperp,curr_qperp,next_qperp);
				// disable jacobian for now
				// jacobian for radial integration
				//const double jacobian=2*M_PI*curr_qperp;
				const double jacobian=1.0;

				// add contribution from this timestep
				//integral_qperp+=weight_qperp*jacobian*differential_rate_p_omega_qperp[ip][iomega][iq];
				integral_qperp+=weight_qperp*jacobian*differential_rate[ip][iomega][iq];
				//update positions for q_perp
				prev_qperp=curr_qperp;
				curr_qperp=next_qperp;
			}

			//std::cout << "ip=" << ip << ", iomega=" << iomega << ", rate=" << integral_qperp << "\n";

			omega_rate_array[iomega]=integral_qperp;

		}
		// WARN << "after the loop. "; 
		// initialise GSL interpolator
		gsl_interp * interp = gsl_interp_alloc( gsl_interp_akima, pOver2T_cut_pos_position_int);
	//	gsl_interp * interp = gsl_interp_alloc( gsl_interp_linear, pOver2T_cut_pos_position_int);
	//	gsl_interp * interp = gsl_interp_alloc(gsl_interp_cspline, pOver2T_cut_pos_position_int);
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_interp_init(interp,  omega_position_array, omega_rate_array, pOver2T_cut_pos_position_int); 

		// integral in omega from -infinity to -omega_over_T_cut
		// if (omega_over_T_min(p_over_T) < -1.*omega_over_T_cut) 
		// WARN << "minimum " << std::min(omega_over_T_min(p_over_T), -1.*omega_over_T_cut); 
			integral_omega+=gsl_interp_eval_integ(interp, omega_position_array, omega_rate_array, omega_over_T_min(p_over_T), std::max(omega_over_T_min(p_over_T), -1.*omega_over_T_cut), acc);
		// integral in omega from omega_over_T_cut to p_over_T/2
		integral_omega+=gsl_interp_eval_integ(interp, omega_position_array, omega_rate_array , std::min(omega_over_T_cut, p_over_T/2.), p_over_T/2., acc);
		// free memory
		gsl_interp_free(interp);
		gsl_interp_accel_free(acc);
		delete [] omega_position_array;
		delete [] omega_rate_array;


		//rate_p[ip]=integral_omega;
		integrated_rate[ip]=integral_omega;
	}
}

// weights for the trapezoid rule on a non-uniform grid
// 1/2 Sum[f[getQfromIndex[position]]Which[0==position,(getQfromIndex[position+1]-getQfromIndex[position]),size_of_array-1==position,(getQfromIndex[position]-getQfromIndex[position-1]),True,(getQfromIndex[position+1]-getQfromIndex[position-1])],{position,0,Qnum-1}]
// 1/2 Sum[f[getQfromIndex[position]]Which[0==position,(next_x-curr_x),size_of_array-1==i,(curr_x-prev_x),True,(next_x-prev_x)],{i,0,size_of_array-1}]
double Tequila::non_uniform_trapeze_weights(int position, int size_of_array, double prev_x, double curr_x, double next_x) {

	double weight=0.5;
	
	if (0 == position) {
		weight*=(next_x-curr_x);	
	}
	else if (size_of_array -1 == position) {
		weight*=(curr_x-prev_x);	
	}
	else {
		weight*=(next_x-prev_x);
	}
	return weight;

}

//Returns T^2 dGamma/(dx domega d^2 qperp)
double Tequila::differential_rate(const double p_over_T, const double omega_over_T, const double qperp_over_T, double *** differential_rate_p_omega_qperp) {


	//tri-linear interpolation
	//somehow I thought this function would be shorter...

	if ((p_over_T<p_over_T_min())||(p_over_T>p_over_T_max())) return 0.0;
	if ((omega_over_T<omega_over_T_min(p_over_T))||(omega_over_T>omega_over_T_max(p_over_T))||(qperp_over_T>qperp_over_T_max())) return 0.0;

	//first, get position in grid of where rate is located
	const double tmp_p=get_index_from_p_over_T(p_over_T);
	const double tmp_omega=get_index_from_omega_over_T(p_over_T,omega_over_T);
	const double tmp_qperp=get_index_from_qperp_over_T(qperp_over_T);

	const int pos_array_p_over_T=floor(tmp_p);
	const int pos_array_omega_over_T=floor(tmp_omega);
	const int pos_array_qperp_over=floor(tmp_qperp);

	//actual positions of the grid points around the desired value
	const double p_over_T_low_val=get_p_over_T_from_index(pos_array_p_over_T);
	const double omega_over_T_low_val=get_omega_over_T_from_index(pos_array_p_over_T,pos_array_omega_over_T);
	const double qperp_over_T_low_val=get_qperp_over_T_from_index(pos_array_qperp_over);
	const double p_over_T_high_val=get_p_over_T_from_index(pos_array_p_over_T+1);
	const double omega_over_T_high_val=get_omega_over_T_from_index(pos_array_p_over_T,pos_array_omega_over_T+1);
	const double qperp_over_T_high_val=get_qperp_over_T_from_index(pos_array_qperp_over+1);

	//value of the rate at the above gridpoints
	const double v000=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T][pos_array_qperp_over];
	const double v001=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T][pos_array_qperp_over+1];
	const double v010=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T+1][pos_array_qperp_over];
	const double v011=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T+1][pos_array_qperp_over+1];
	const double v100=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T][pos_array_qperp_over];
	const double v101=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T][pos_array_qperp_over+1];
	const double v110=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T+1][pos_array_qperp_over];
	const double v111=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T+1][pos_array_qperp_over+1];

	//fraction of each corner to use
	const double frac_p_over_T=(p_over_T-p_over_T_low_val)/(p_over_T_high_val-p_over_T_low_val);
	const double frac_omega_over_T=(omega_over_T-omega_over_T_low_val)/(omega_over_T_high_val-omega_over_T_low_val);
	const double frac_qperp_over_T=(qperp_over_T-qperp_over_T_low_val)/(qperp_over_T_high_val-qperp_over_T_low_val);

	//get the value
	const double v00=v000*(1-frac_qperp_over_T)+v001*frac_qperp_over_T;
	const double v01=v010*(1-frac_qperp_over_T)+v011*frac_qperp_over_T;
	const double v10=v100*(1-frac_qperp_over_T)+v101*frac_qperp_over_T;
	const double v11=v110*(1-frac_qperp_over_T)+v111*frac_qperp_over_T;

	const double v0=v00*(1-frac_omega_over_T)+v01*frac_omega_over_T;
	const double v1=v10*(1-frac_omega_over_T)+v11*frac_omega_over_T;

	double res=v0*(1-frac_p_over_T)+v1*frac_p_over_T;

	return res;

}

////Rate
//double ERateColl::rate(const struct ERateParam &rate_params, const Pythia8::Vec4 &p0, const int &id0)
double Tequila::rate_inel(double energy, double temp, double * rate_p)
{

//	const double temp=rate_params.T();
//
//	if (id0!=fId0) { 
//		// the incoming particle is not a gluon and this rate doesn't apply
//		return 0. ;
//	} 
//	if (temp/p0.e() >  ERateParam::MaxTOverP ) {
//		return 0. ;  // Don't evolve if the energy is too small
//	}

	//energy/temperature ratio of the incoming jet
	//const double energy_over_T=p0.e()/temp;
	const double energy_over_T=energy/temp;

	//if outside of tabulated range, set to 0 (that is, assume untabulated rate is tiny and can be ignored)
	if ((energy_over_T<p_over_T_min())||(energy_over_T>p_over_T_max())) return 0.0;
	//get the real-valued index
	double a = get_index_from_p_over_T(energy_over_T);
	//get the actual index of the array
	int n_p = int(a);
	//get the remainder
	a -= n_p;
	//get rate for energy/T
	double result = (1-a) * rate_p[n_p] + a * rate_p[n_p+1];

	//a factor of temperature is missing from the integrated rate
	result*=temp;


	return result;

}

//Given a value of parton momentum p, sample the rate in omega and qperp
void Tequila::sample_dgamma_dwdq(double p, double T, double *** differential_rate_p_omega_qperp, double &w, double &q) {
	//   double lam = pp.lambda(p0) ;
	int ntry = 0  ;
	const int ntry_max = 1000000;   
	
	//helper variables
	//const double qperp_over_T_val_min=0.0;
	//const double qperp_over_T_val_max=qperp_over_T_max();
	const double p_over_T=p/T;
	const double omega_over_T_neg_min=omega_over_T_min(p_over_T);
	const double omega_over_T_neg_max=-omega_over_T_cutoff;
	const double omega_over_T_pos_min=omega_over_T_cutoff;
	const double omega_over_T_pos_max=p_over_T/2.0;
	const double omega_over_T_neg_span=(omega_over_T_neg_max-omega_over_T_neg_min);
	const double omega_over_T_pos_span=(omega_over_T_pos_max-omega_over_T_pos_min);

	const double max_rate=maximum_rate_p(p_over_T, differential_rate_p_omega_qperp);

	while (ntry < ntry_max) {
		//      fRateTable.sample_in_x1range(-GSL_DBL_MAX,lam/(2.*T), rng, params, w, q) ;
		//double r = max_rate*rng(params) ;
		double r = max_rate*ZeroOneDistribution(*GetMt19937Generator()) ;
		//      double v = (1. -  T*w/p0.e())/rratio(w) ;
		//const double qx_over_T_tmp=qperp_over_T_val_max*rng(params);
		//const double qy_over_T_tmp=qperp_over_T_val_max*rng(params);
		//const double q_over_T_test=sqrt(qx_over_T_tmp*qx_over_T_tmp+qy_over_T_tmp*qy_over_T_tmp);

		// sample omega_over_T
		// a bit more complicated due to the omega_over_T_cutoff around 0
		// With[{tmp = (rnd*(xnegspan + xposspan))},
		// If[tmp > xnegspan, xposmin + (tmp - xnegspan), xnegmin + tmp]
		//  ]
		//const double rnd2=rng(params);
		const double rnd2=ZeroOneDistribution(*GetMt19937Generator());
		const double tmp=rnd2*(omega_over_T_neg_span+omega_over_T_pos_span);
		double omega_over_T_test;
		if (tmp < omega_over_T_neg_span) omega_over_T_test=omega_over_T_neg_min+tmp;
		else omega_over_T_test=omega_over_T_pos_min+(tmp-omega_over_T_neg_span);

		// sample q_over_T, remembering it is a radial variable
		// also, q_over_T must be smaller than abs(omega_over_T_test)
		//const double q_over_T_sqr=omega_over_T_test*omega_over_T_test*rng(params);
		//const double q_over_T_test=sqrt(q_over_T_sqr);
		const double q_over_T_test=0.0;  //let's use qperp/T=0 for now

		double v = differential_rate(p_over_T, omega_over_T_test, q_over_T_test, differential_rate_p_omega_qperp);

		//safety check: did we pick the maximum rate correctly??
		if (v > max_rate) {
			std::cout << "Function \"maximum_rate_p()\" apparently does not return the maximum of the rate... This is bad. Fix it.\n";
			std::cout << "current guess for maximum=" << max_rate << " & sampled rate at omega/T="<< omega_over_T_test << "&q_perp/T=" << q_over_T_test << " is " << v << "\n";
			exit(1);
			//assert(false);
		}

		if (r < v) {
			w=omega_over_T_test*T;
			q=q_over_T_test*T;
			return ;
		}
		else {
			ntry++ ;
		}
	}
	w=omega_over_T_cutoff*T;
	q=0.0;
	std::cout << "*** ERateColl::sample_dgamma_dwdq *** Failed to find a sample "
		"after "<< ntry  << " iterations! Returning with w,q = " << w<<","<< q << std::endl;
}

