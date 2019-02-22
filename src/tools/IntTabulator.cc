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

void IntTabulator::Tabulator_dGamma_domega_qperp(std::string path, process_type process)
{
	std::ofstream table2d_out((path+"elastic_rate_table"+GetProcessString(process)+".dat").c_str()); 
	double wp, w, q, qp; 
	for (size_t i = 0; i <= Nw; i++)
	{
		wp = ((double)i)*(atan(wMax)-atan(omegaMin))/Nw+atan(omegaMin); 
		w = tan(wp); 
		// double qperpMax = 2.*sqrt(kMax*kMax + kMax*w); 
		for (size_t j = 0; j <= Nq; j++)
		{
			qp = ((double)j)*(log(qperpMax)-log(muperp0))/Nq+log(muperp0); 
			q = exp(qp); 
			table2d_out << dGamma_domega_qperp_forTab(w, q, process) << " "; 
		}
		table2d_out << "\n"; 
	}
}

double dGamma_domega_qperp_k(double k, void *params)
{	
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	p->f_k = k; 
	gsl_function F; 
	switch(p->f_process)
	{
		case gg: F.function = dGamma_domega_qperp_k_phi_gg; 
			 break; 
		case gq: F.function = dGamma_domega_qperp_k_phi_gq; 
			 break; 
		case qg: F.function = dGamma_domega_qperp_k_phi_qg; 
			 break; 
		case qq: F.function = dGamma_domega_qperp_k_phi_qq; 
			 break; 
		case qqp: F.function = dGamma_domega_qperp_k_phi_qqp; 
			 break; 
		case qqb: F.function = dGamma_domega_qperp_k_phi_qqb; 
			 break; 
		default: std::cout << "The process determination is wrong! The process is " << p->f_process; 
			 break; 
	}
	F.params = p; 
	double result, err; 
	gsl_integration_qag(&F, 0, 2.*M_PI, cls.ErrAbs, cls.ErrRel, cls.NWorkSpace, cls.fKey, cls.Space_phi, &result, &err); 
	return result/(2.*M_PI); 
}

double dGamma_domega_qperp_k_phi_gg(double phi, void *params)
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
	M2 = (double)(CA*CA)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nB(k)*(1+nB(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_gq(double phi, void *params)
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
	M2 = (double)(dF*CF*CA/dA*6)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_qg(double phi, void *params)
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
	M2 = (double)(CF*CA)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nB(k)*(1+nB(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_qqp(double phi, void *params)
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
	M2 = (double)(dF*CF*CF/dA*4)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_qq(double phi, void *params)
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
	M2 = (double)(dF*CF*CF/dA)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_qqb(double phi, void *params)
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
	M2 = (double)(dF*CF*CF/dA)*2.*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(2.*nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

