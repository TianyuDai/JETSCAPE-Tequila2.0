#ifndef TEQUILA_H
#define TEQUILA_H

#include "JetEnergyLossModule.h"
#include "JetScapeConstants.h"
#include "gsl/gsl_integration.h"
#include "IntTabulator.h"
#include "lorentz.h"
#include <vector>
#include "random.h"
#include <math.h>
#include <gsl//gsl_spline.h>
#include "Tequila/evolve9_text_noExp.h"

using namespace Jetscape;

class Tequila : public JetEnergyLossModule<Tequila> //, public std::enable_shared_from_this<Tequila>
{
  private: 
	double alpha_s;
	double g;
   	double pcut;
   	double M; 
   	double muperp; 
   	double mu_scale; 
   	// AMY rates are calculated in p/T > AMYpCut
  	static constexpr double AMYpCut = 0.1;

  	double Q0;
  	double alpha_EM;
  	double hydro_Tc;
  	double omega_over_T_cutoff;
	double Lambda; 

  	const double nc = 3.; 
	const double ln2 = log(2); 

	const int Nsteps = 500; 
	const static size_t nElasProcess = 6; 
	const static size_t nTotalProcess = 20; 
	
	int iqperp0 = 0; 
	double max_omega_rate; 
	double max_qperp_rate; 
	double casimir[nElasProcess] = {CA*CA/(24.*M_PI), dF*CF*CA/dA/(8.*M_PI), CF*CA/(24.*M_PI), dF*CF*CF/dA/(72.*M_PI), dF*CF*CF/dA/(18.*M_PI), dF*CF*CF/dA/(72.*M_PI)};		// The coefficient for total rate integration of large omega

	// gg, gq, qg, qq, qqp, qqb
	// const static size_t nSplittingProcess = 9; 
	double splitting_c_lambda[nElasProcess] = {4.*CA*CA, 16.*CA, 4.*CF*CA, 16.*CF, 16.*CF, 16.*CF}; 
	double splitting_c_p[nElasProcess] = {5./6*CA*CA, 2.*CF-4.*CA, CF*CF/2-CF*CA, -4.*CF, -4.*CF, 4.*CF*(6.*CF-3.*CA-1./3)}; 
	double splitting_c_ln[nElasProcess] = {-2.*CA*CA, 4.*CF-8.*CA, CF*CF-2.*CF*CA, 8.*CF*(2.*CF-CA-1.), -8.*CF, -8.*CF*(2.*CF-CA+1.)}; 
	struct tables
	{
		double rate = 0.; 
		double x[Nw+1], y[Nw+1]; 
		double xa[Nw+1], ya[Nq+1]; 
		double rate_qperp[Nw+1][Nq+1];  
	}; 
	vector <tables> elasticTable; 
	double *za = (double*) malloc(nElasProcess * (Nw+1) * (Nq+1) * sizeof(double)); 
	// double *za = new double(nElasProcess * (Nw+1) * (Nq+1)); 

	bool is_empty(std::ifstream& pFile); 
	bool is_exist (const std::string& name); 
	
	// size of the array containing the differential rate in p, omega and q_\perp
  	// those number in p and omega are, and must be, the same used in Guy Moore's code that is used to evaluate the collinear rates
	// the number of point in q_\perp is arbitrary at this point
	static const int nb_points_in_p=NP, nb_points_in_omega=NK;
	// Addionnal grid in q_perp, not defined in Guy Moore's code
	static const int nb_points_in_qperp=4;
	// the minimum and maximum parton energy p and parton energy loss omega are set in Guy Moore's code as well and must be the same here
	double p_over_T_min() { return 4.01; };
	double p_over_T_max() { return (p_over_T_min()*pow(1.04119,nb_points_in_p-1)); };
	double omega_over_T_min(const double p_over_T) { return -12.0; };
	double omega_over_T_max(const double p_over_T) { return p_over_T + 0.2 * (nb_points_in_omega -1 - 320); };
	// for future use
	//double qperp_over_T_min(const double p_over_T, const double omega_over_T) { return -5.0; };
	//double qperp_over_T_min() { return -10.0; };
	double qperp_over_T_max() { return 10.0; };
	//given the above tabulation of p/T, this function returns the index for a given p/T
	//the index is a real number, such that "int(index)" gives the position in the array and "index-int(index)" gives the remainder
	double get_index_from_p_over_T(const double p_over_T) { return 24.7743737154026 * log( p_over_T * .2493765586034912718l ); };
	double get_index_from_omega_over_T(const double p_over_T, const double omega_over_T); //Big function, defined elsewhere instead
	double get_index_from_qperp_over_T(const double qperp_over_T) { 
		const double qmin=0.0;
		const double qmax=qperp_over_T_max();
		return (((nb_points_in_qperp-1)*(qperp_over_T-qmin))/(qmax-qmin));
	}
	double get_p_over_T_from_index(const int index_p_over_T) { return (p_over_T_min()*pow(1.04119,index_p_over_T)); };
	double get_omega_over_T_from_index(const int index_p_over_T, const int index_omega_over_T); //Big function, defined elsewhere instead
	// Assume uniform grid in q_perp for now
	double get_qperp_over_T_from_index(const int index_qperp_over_T) { 
		const double qmin=0.0;
		const double qmax=qperp_over_T_max();
		return qmin+(qmax-qmin)/(nb_points_in_qperp-1)*index_qperp_over_T;
	};
	// arrays containing the differentail and integrated rates
	//double differential_rate_p_omega_qperp[nb_points_in_p][nb_points_in_omega][nb_points_in_qperp];
	//double rate_p[nb_points_in_p];
	double *** differential_rate_qqg_p_omega_qperp, *** differential_rate_ggg_p_omega_qperp, *** differential_rate_gqq_p_omega_qperp;
	double *rate_qqg_p, * rate_ggg_p, *rate_gqq_p;
	//maximum of the differential rate, used to sample the rate

	void allocate_memory_for_radiative_rate_table();

	double non_uniform_trapeze_weights(int position, int size_of_array, double prev_x, double curr_x, double next_x);

	// load differential rate from file into member array "differential_rate_p_omega_q"
	void load_differential_rate(const double alpha_s, const double alpha_EM, const int Nf, const std::string location_of_collinear_rates);
	// fill member array "rate_p" by integrating the content of "differential_rate_p_omega_q"
	//void evaluate_integrated_rate(double omegaMin);
	void evaluate_integrated_rate(double omega_over_T_cutoff, double *** differential_rate_gqq_p_omega_qperp, double * rate_gqq_p);

	//  double differential_rate(const double p_over_T, const double omega_over_T, const double qperp_over_T);
	double differential_rate(const double p_over_T, const double omega_over_T, const double qperp_over_T, double *** differential_rate_p_omega_qperp);
	double rate_inel(double energy, double temp, double * rate_p);
	void sample_dgamma_dwdq(double p, double T, double *** differential_rate_p_omega_qperp, double &w, double &q);
	//  void sample_dgamma_dwdq(const struct ERateParam &pp, const Pythia8::Vec4 &p0, double (*rng)(void*params), void *params,  double &w, double &q);
	//
	//educated guess on the highest value of the rate
	//double *maximum_differential_rate;
	double maximum_rate_p(double p_over_T, double *** differential_rate_p_omega_qperp) {
		//the collinear rate should presumably be maximum at qperp/T=0 and either omega/T=+/-omega_over_T_cut 
		const double val1=2.*differential_rate(p_over_T,omega_over_T_cutoff,0,differential_rate_p_omega_qperp);
		const double val2=2.*differential_rate(p_over_T,-1*omega_over_T_cutoff,0,differential_rate_p_omega_qperp);
		return val1 > val2 ? val1 : val2;
	};

	

  public:

	Tequila();
	virtual ~Tequila();

	double rate_conv(process_type process, double T, double pRest); 
  	void LoadElasticTables(); 
	double Interpolator_dGamma_domega(double omega, process_type process); 
	double Extrapolator_dGamma_domega(double omega, process_type process); 
	double Interpolator_dGamma_domega_qperp(double omega, double qperp, process_type process); 
	double splittingF(double x, process_type process);
        double splitting_gg1(double x); 
        double splitting_gg2(double x);  
	double splittingRateOmega(double pRest, double x, double T, process_type process); 
	// double splittingRateQperp(double pRest, double qperp, double omega, double T, process_type process); 

  	void Init();
	void Finish(); 
  	void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  	process_type DetermineProcess(double pRest, double T, double deltaTRest, int Id); 
  	FourVector Momentum_Update(double omega, double qperp, double T, FourVector pVec); 
  	double TransverseMomentum_Transfer(double pRest, double omega, double T, process_type process); 
  	double Energy_Transfer(double pRest, double T, process_type process); 
	double TransverseMomentum_Transfer_Split(double pRest, double omega, double T, process_type process); 
  	double Energy_Transfer_Split(double pRest, double T, process_type process); 
	double xSampling(double pRest, double kRest, double T, process_type process); 
	double thermalDistributionSampling(process_type process); 
  	
  	double qpara(double E, double T, int id); 
  	double qperp(double E, double T, int id); 
  	FourVector Langevin_Update(double dt, double T, FourVector pIn, int Id); 
  	
  protected:
  	uniform_real_distribution<double> ZeroOneDistribution;
	string PathToTables;
	string path_to_tables;
}; 

double dGamma_domega(double omega, void *params); 
double dGamma_domega_qperp(double qperp, void *params); 
double dGamma_domega_qperp_k(double k, void *params); 
double dGamma_domega_qperp_k_phi(double phi, void *params); 

#endif

