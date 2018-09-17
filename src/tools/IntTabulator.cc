#include "IntTabulator.h"
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

IntTabulator::IntTabulator()
{
	Space_phi = gsl_integration_workspace_alloc(NWorkSpace); 
	Space_k = gsl_integration_workspace_alloc(NWorkSpace); 
	Space_qperp = gsl_integration_workspace_alloc(NWorkSpace); 
	Space_omega = gsl_integration_workspace_alloc(NWorkSpace); 
}

IntTabulator::~IntTabulator()
{
	gsl_integration_workspace_free(Space_phi); 
	gsl_integration_workspace_free(Space_k); 
	gsl_integration_workspace_free(Space_qperp); 
	gsl_integration_workspace_free(Space_omega); 
}

std::string IntTabulator::GetProcessString(int enumVal)
{
  	return ProcessStrings[enumVal];
}
	
void IntTabulator::Gamma(double muperp, std::string path, process_type process)
{
	struct f_params p; 
	// cls_muperp = muperp; 
	p.f_muperp = muperp; 
	p.f_process = process; 
	gsl_function F; 
	F.function = dGamma_domega; 
	F.params = &p; 
	double result, err; 
	gsl_integration_qagi(&F, ErrAbs, ErrRel, NWorkSpace, Space_omega, &result, &err); 
	std::ofstream table0d_out((path+"rate0d_table"+"_muperp"+std::to_string(muperp)+GetProcessString(process)+".dat").c_str()); 
	table0d_out << result << "\n"; 
	return; 
}

double IntTabulator::dGamma_domega_forTab(double omega, double muperp, process_type process)
{
	struct f_params p; 
	p.f_omega = omega; 
	p.f_process = process; 
	gsl_function F; 
	F.function = dGamma_domega_qperp; 
	F.params = &p; 
	double result, err; 
	gsl_integration_qagiu(&F, muperp, ErrAbs, ErrRel, NWorkSpace, Space_qperp, &result, &err); 
	return result; 
}

void IntTabulator::Tabulator_dGamma_domega(double muperp, std::string path, process_type process)
{
	std::ofstream table1d_out((path+"rate1d_table"+"_muperp"+std::to_string(muperp)+GetProcessString(process)+".dat").c_str()); 
	double w; 
	for (size_t i = 0; i <= Nw; i++)
	{
		w = (double)i*(kMax-omegaMin)/Nw+omegaMin; 
		// Set the interpolation function to be log(dGamma_domega*(w^2+1)) to make the function linear-like
		table1d_out << w << " " << log(dGamma_domega_forTab(w, muperp, process)*(w*w+1)) << "\n"; 
	}
}

double dGamma_domega(double omega, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	p->f_omega = omega; 
	gsl_function F; 
	F.function = dGamma_domega_qperp; 
	F.params = p; 
	double result, err; 
	// double qperpMax = (2.*sqrt(cls.kMax*cls.kMax + cls.kMax*omega)); 
	// if (p->f_muperp > cls.qperpMax) return 0.; 
	gsl_integration_qagiu(&F, p->f_muperp, cls.ErrAbs, cls.ErrRel, cls.NWorkSpace, cls.Space_qperp, &result, &err); 
	return result; 
}

double IntTabulator::dGamma_domega_qperp_forTab(double omega, double qperp, process_type process)
{
	struct f_params p; 
	p.f_qperp = qperp; 
	p.f_omega = omega; 
	p.f_process = process; 
	gsl_function F; 
	F.function = dGamma_domega_qperp_k; 
	F.params = &p; 
	double result, err; 
	double q = sqrt(qperp*qperp + omega*omega); 
	double lowLimit = (q - omega) / 2.; 
	gsl_integration_qagiu(&F, lowLimit, ErrAbs, ErrRel, NWorkSpace, Space_k, &result, &err); 
	return result; 
}

void IntTabulator::Tabulator_dGamma_domega_qperp(double muperp, std::string path, process_type process)
{
	// cls_muperp = muperp; 
	std::ofstream table2d_out((path+"rate2d_table"+"_muperp"+std::to_string(muperp)+GetProcessString(process)+".dat").c_str()); 
	double wp, w, q, qp; 
	for (size_t i = 0; i <= Nw; i++)
	{
		wp = ((double)i)*(atan(kMax)-atan(omegaMin))/Nw+atan(omegaMin); 
		w = tan(wp); 
		for (size_t j = 0; j <= Nq; j++)
		{
			qp = ((double)j)*(log(qperpMax)-log(muperp))/Nq+log(muperp); 
			q = exp(qp); 
			table2d_out << log(dGamma_domega_qperp_forTab(w, q, process)) << " "; 
		}
		table2d_out << "\n"; 
	}
}

double dGamma_domega_qperp(double qperp, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	p->f_qperp = qperp; 
	gsl_function F; 
	F.function = dGamma_domega_qperp_k; 
	F.params = p; 
	double result, err, lowLimit; 
	double omega = p->f_omega; 
	double q = sqrt(qperp*qperp + omega*omega); 
	lowLimit = (q - omega) / 2.; 
	gsl_integration_qagiu(&F, lowLimit, cls.ErrAbs, cls.ErrRel, cls.NWorkSpace, cls.Space_k, &result, &err); 
	return result; 
}

double dGamma_domega_qperp_k(double k, void *params)
{	
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	p->f_k = k; 
	gsl_function F; 
	switch(p->f_process)
	{
		case gg: F.function = dGamma_domega_qperp_k_phi_gggg; 
			 break; 
		case gq: F.function = dGamma_domega_qperp_k_phi_gqgq; 
			 break; 
		case qg: F.function = dGamma_domega_qperp_k_phi_qgqg; 
			 break; 
		case qq: F.function = dGamma_domega_qperp_k_phi_qqqq; 
			 break; 
		default: std::cout << "The process determination is wrong! The process is " << p->f_process; 
			 break; 
	}
	F.params = p; 
	double result, err; 
	gsl_integration_qag(&F, 0, 2.*M_PI, cls.ErrAbs, cls.ErrRel, cls.NWorkSpace, cls.fKey, cls.Space_phi, &result, &err); 
	return result/(2.*M_PI); 
}

double dGamma_domega_qperp_k_phi_gggg(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp, omega; 
	qperp = p->f_qperp; 
	k = p->f_k; 
	omega = p->f_omega;  
	kp = k + omega; 
	q = sqrt(qperp*qperp + omega*omega); 
	t = -1.*pow(qperp, 2); 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./4./pow(2.*M_PI, 3)*qperp/q; 
	M2 = (double)(cls.CA*cls.CA)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nB(k)*(1+nB(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_gqgq(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp, omega; 
	qperp = p->f_qperp; 
	k = p->f_k; 
	omega = p->f_omega; 
	kp = k + omega; 
	q = sqrt(qperp*qperp + omega*omega); 
	t = -1.*pow(qperp, 2); 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./4./pow(2.*M_PI, 3)*qperp/q; 
	M2 = (double)(cls.dF*cls.CF*cls.CA/cls.dA*6)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_qgqg(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp, omega; 
	qperp = p->f_qperp; 
	k = p->f_k; 
	omega = p->f_omega;  
	kp = k + omega; 
	q = sqrt(qperp*qperp + omega*omega); 
	t = -1.*pow(qperp, 2); 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./4./pow(2.*M_PI, 3)*qperp/q; 
	M2 = 2.*(double)(cls.CF*cls.CA)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nB(k)*(1+nB(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_qqqq(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp, omega; 
	qperp = p->f_qperp; 
	k = p->f_k; 
	omega = p->f_omega;  
	kp = k + omega; 
	q = sqrt(qperp*qperp + omega*omega); 
	t = -1.*pow(qperp, 2); 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./4./pow(2.*M_PI, 3)*qperp/q; 
	M2 = (double)(cls.dF*cls.CF*cls.CF/cls.dA*4)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

