#ifndef evolve9
#define evolve9

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#define PI          3.14159265358979l
#define Gamma_E     0.57721566490153286061l
#define SIXTH       0.16666666666666667l
//#define NP          230
//#define NK          381
const int NP=230;
const int NK=381;

/* NP and NK are the number of p and k values needed to build a table
   of dGamma(p,k)/dkdx sufficiently refined to do interpolation.  Don't
   change them unless you know what you are doing. */
#define ABS(x)      ( (x) > 0 ? (x) : -(x) )
#define PauliBlock(x) ( (x)>0 ? 1/(1+exp(-(x))) : exp(x)/(exp(x)+1) )
#define BoseStim(x)   ( (x)>0 ? 1/(1-exp(-(x))) : exp(x)/(exp(x)-1) )
#define  READ_LETTER(x,y) while ( (( (x) = getc (y) ) < 'a' ||       \
                    (x) > 'z' ) && ( (x) < 'A' || (x) > 'Z' ) )    
#define READ_TO_EOLN(ff) {char junkc; while ( (junkc = getc(ff))!= '\n' );}

/* Structure to store information about couplings, momentum range and
   discretization, group theoretic factors, etc */
typedef struct
{
  double df;
  double da;
  double cf;
  double ca;
  int    Nc;
  int    Nf;
  int    Bethe_Heitler;
  int    BDMPS;
  int    include_gluons;
  int    photon_22;
  int    do_collisional;
  int    collisional_only;
  int    do_compton;
  int    do_comptong;
  char   in_fname[1000]; /* Stores name of file where pre_saved
			   info on dGamma/dxdk is stored. */

  double alpha_s;
  double alpha;
  double delta_x;
  double dx_max;
  double dp;
  double p_min;
  double p_max;
  long   n_p;
  long   n_pmin;
  double k_min;
  double k_max;
  long   n_k;
  long   n_kmin;

} Gamma_info;

/* Structure to store information about splitting functions, to interpolate
   from */
typedef struct
{
  double dGamma[NP][NK];
  double dGamma_gqq[NP][NK];
  double dGamma_ggg[NP][NK];
  double dGamma_em[NP][NK];
} dGammas;


  
void build_table ( Gamma_info * dat , dGammas *Gam );
void write_table ( Gamma_info *dat , dGammas *Gam );
void read_table ( Gamma_info *dat , dGammas *Gam );
double use_table ( double p , double k , double dGamma[NP][NK] , int which_kind );

#endif
