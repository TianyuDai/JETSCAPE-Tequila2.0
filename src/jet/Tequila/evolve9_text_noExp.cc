//This is Guy Moore's code

/* This version is modified to incorporate an approximate kinematic
   cutoff to discard scatterings which are not small-angle.  Namely,
   choosing p to be the largest and k to be the smallest energy of an
   external state, we require k_perp < k with k_perp the component
   of k transverse to p, or h^2 < p^2 k^2.  We do this by imposing a cutoff
   on the h integration, \int h dh h*F -> \int h dh h*F exp(-h^2/k^2p^2)
   (note that \int d^2 h exp(-h^2/k^2 p^2) = pi k^2 p^2 = area of circle
   with h < kp.  So this choice roughly treats the size of the region
   allowed in the integral correctly.)
   In impact parameter space, this means that
   f(b->0) should be replaced by
   int b db f(b) (k^2 exp(-k^2 b^2/4) / 2)

   The remaining challenge is the units; h is given units of g T^2
   whereas pk has units of T^2, so one must use an explicit value of g
   in evaluating the kinematically cut off matrix elements, therefore
   we have to input a value for alpha_s before building the table
   of gammas. */

/* Program to evaluate the parton energy loss rate on traversing a
   medium, and to evolve the Fokker-Planck equation for the probability 
   distribution with energy as a function of time in the medium (coarse 
   graining the set of allowed energies).

   This version is intended for evolution through a series of slabs of 
   medium, each of different temperature and boost factor.  Since the
   temperature changes, we have to put all p,k momenta in GeV rather than
   in "natural" units of T.

   As a temporary stop-gap way to include the destructive interference with
   vacuum radiation indicated by a long formation time, this version allows
   the inclusion of a formation-time dependent suppression factor, which
   is 
   max_tau^2 / ( max_tau^2 + tau^2 )
   The user inputs max_tau in units of fermis--it should be comparable to
   the system age.  At present this is a "hack" and a more careful treatment
   should be performed.

   It works by solving the (Migdal) integral equation for the
   differential splitting rate, in impact parameter space.  This version
   in addition includes Compton/pair type processes, in the approximation
   that the exchange energy vanishes and they are simply identity 
   changing interactions.  It also includes elastic scattering.  These
   options are turned on or off by user input, no recompilation required.

   Starting from Arnold, Moore, Yaffe, the splitting rate is

   dGamma(p,k) / dk dx

   dGamma/dk dx = 1/(2 p^2 d_s) (1 +- f(k))(1 +- f(p-k)) 
                  (d_s C_s g^4)/(8 PI) [FACTOR1] [FACTOR2]

   FACTOR1 = 1/(p^3 (p-k)^3 k^3) times | ( k^4 + p^4 + (p-k)^4 ) g->gg
                                       | k(p-k)(k^2 + (p-k)^2)   g->qq
				       | p(p-k)(p^2 + (p-k)^2)   q->gq

   The 2 d_s is a degeneracy; the 1/2 in the paper's expression is a
   symmetry factor handled by our always making k the gluon.
   FACTOR2 = solution of integral equation, which is (p^4 m_D^4)/(g^4 T)
   times solution to a dimensionless integral equation,

   \int d^2 h/(2 PI)^2 2 h * Re F, 

   2 h = i (dE/g^2 T) F + \int d^2 q/(2 PI)^2 ( 1/q^2 - 1/(q^2 + 1) ) X
           ( C_1 (F(h)-F(h+q)) + C_2 (F(h)-F(h+q k/p)) + C_3 (same,(p-k)/p))

   C_1,C_2,C_3 = C_(particle) - C_A/2
   dE/g^2 T = h^2 p T/(2 k (p-k)) + A
   A = -(m_p^2 / g^2 T^2) 1/2p + (m_k^2 / g^2 T^2) 1/2k + (same, (p-k))

   Turn it into the differential equation,

   0 = i ( A - B \nabla^2 ) f(b) + D K(b) f
   B = p T / (k (p-k) )
   D K(b) = (1/2 PI) ( C_1 K(b) + C_2 K(b k/p) + C_3 K(b (p-k)/p) )
   K(b) = K_0(b) + ln(b/2) + Gamma_E  [ K_0 the modified Bessel function ]

   The 1-D version, writing vector{f} = \vec{b} f, is (b!= 0)

   0 = i A f - i B [3/b] f' - i B f''  + D K(b) f
   
   or
   f'' = - [3/b] f' + (A/B) f - i (DK(b)/B) f
   
   In this form, we want 4/Pi times real part if Im part shows 1/b^2
   behavior, times 1/B to make up for rescaling of the equation.

   The result of the integral, [ dGamma / dk dx ], is then used to
   evolve the probability distribution of particle energy,

   dP(p)/dt = \int dk (  P(p+k) [dGamma(p+k,k)/dk dx] 
                       - P(p)   [dGamma(p,k)/dk dx] )

   [here x and t are used interchangeably....]
   where the two terms are classic gain and loss terms due to emission
   of a quantum of momentum k.

   For the case of g->gg, we have to restrict to $k < p/2$ to avoid
   double counting of final states; this can be done for g->qqbar as
   well, as "quarks" means "q plus qbar", but there needs to be a 
   compensating factor of 2 as q and qbar are distinguishable.
   When including gluons, there are extra gain terms (emitted gluon,
   emitted q and qbar, both emitted gluons in g->gg).

   We accelerate the solution of the Migdal equation by storing enough 
   values to use for interpolation, and invoking an interpolation 
   procedure.  The interpolation extracts the dominant behavior,

   dGamma(p,k) \sim (population functions) * [ 1/k, 1/p, or p/k(p-k) ]
   and fits the residual.

   To solve the equation dP(p)/dx = ..., we do the following;

   First, the allowed values of p are discretized with discretization 
   Delta.  Then, [dGamma(p)/dk dx] is evaluated at each discrete value of
   p, within some range [p_min,p_max], and at k equaling ODD multiples
   of Delta.  The evolution is replaced by,

   dP(p)/dt = sum_n 2 Delta ( P(p+n Delta) [dGamma(p+n Delta,k)/dk dx]
                            - P(p) [dGamma(p,k)/dk dx] )
   
   When gluons are also tracked, a similar expression exists for gluons,
   with a gluon gain term due to the radiated gluon from the quark, and 
   a quark gain term from the splitting of a gluon into two quarks.

   The use of odd multiples of Delta for values of k, is necessary to
   avoid the quadratic singularity in dGamma(p,k)/dk at k=0, which
   requires that the point k=0 be excluded, without causing a linear
   in Delta error; having evenly spaced values of k is required because
   the region close to k=0 is important to the momentum diffusion
   coefficient.

   The use of odd multiples means that the emitted gluon always has a
   momentum which is an odd multiple of Delta.  This is "fixed" by
   smearing the final momentum of the gluon, so that rather than
   equaling k, it is k-Delta with 25% weight, k with 50% weight,
   and k+Delta with 25% weight.  The error due to this and due to 
   discretization vanishes, for smooth momentum distributions,
   as Delta^2.

   The infinitesimal time dP(p,t)/dt result is converted into a finite
   difference time evolution algorithm as follows;
   1) dP(p,t)/dt is determined;
   2) a maximum "safe" dt is determined;
   3) P is evolved to leading order accuracy to time t+dt
   4) dP(p_proj,t+dt)/dt is evaluated
   5) the evolution is taken to be,
      P(p,t+dt) = P(p,t) + (dt/2) ( dP(p,t)/dt + dP(p_proj,t+dt)/dt )
   which is just the second-order Runge-Kutta algorithm.

   The criteria for acceptable dt are, that it satisfies both
   a) dt * dP(p,t) < (NUMBER) * P(p,t) AT EACH p
   b) the same holds as an integral relation but with a tighter NUMBER.
   The choice of the NUMBER is determined by the user; the error
   shrinks quadratically as it is made small.

   When a series of checkpoints or printout points are specified, the
   time steps are also chosen to go up exactly to the edge of a checkpoint
   or printout point, with the last dt being shrunk to fit.

   The program also includes the photons produced by 2 <-> 2 processes
   (Compton and pair annihilation) when asked to by the user.
   They are included in the approximation that
   the photon energy does not differ from the energy of the particle
   producing it.  This is done by adding a special contribution to
   dGamma(p,k)/dk dx   right at k=p (or k = p +- delta_p if k=p is
   not allowed by the "staggering" of the k values), of magnitude equal
   to the total rate, scaled correctly by the k spacing.

   2 <-> 2 Compton-type quark gluon scattering is included in the same way.
*/
/* #define  NOISY */
#include "evolve9_text_noExp.h"

/* ----------My favorite file opener.---------------- */

FILE *openfile ( char filename[] , std::string mode )
     /*  Gets a file's name and opens it in the mode given. */
{
  FILE *result;
  long  i;

  i = 0;
  while ( ( ( filename[i] = getchar() ) == '\n' ) || filename[i] == ' ' );
  while ( ( filename[++i] = getchar () ) != '\n' )
    if ( filename[i] == ' ' ) i--; /* read to end of line, throwing 
				      out empty spaces. */
  filename[i] = '\0';
  result = fopen ( filename , mode.c_str() );

  return ( result );
}

/* ============================================================= */
/* Equipment for evaluating dGamma/dkdx the splitting kernel.    */
/* ============================================================= */

double K ( double z , int just_K0 )
     /* Determines K(z) = K_0(z) + ln(z/2) + Gamma_E either by
	power series or by asymptotic series.  10 digit accurate */
     /* I = I_0(z) - 1.  (or I_0(z) to get just_K0).*/
     /* This is not very elegant.  It uses the power series except 
	for large z, where it uses an asymptotic series.  The power
	series suffers from big, cancelling terms, hence only 10 
	digit accuracy.  Is there a better way using recurrence relations? */
{
  double sofar , I , coeff , coeff2 , logplus;
  int    j;
  
  if ( z < 0 ) return ( K ( -z , just_K0 ) );
  if ( z == 0 ) return ( 0 ); /* Safety first. */
  if ( z < 11 ) /* Use power series */
    {
      logplus = log ( 0.5 * z ) + Gamma_E;
      z *= 0.25 * z; /* z*z/4 is what series is in. */
      I = ( just_K0 ? 1 : 0 );
      coeff2 = 1; /* keeps track of 1+1/2+1/3+... which is added to logplus. */
      j = 1;
      coeff = z;  /* will become z^j / (j!)^2 */
      sofar = 0;
      while ( coeff > 1.0e-12 )
	{
	  I += coeff;
	  sofar += coeff * coeff2;
	  j++;
	  coeff *= z / ( j * j ); /* z^j / (j!)^2 */
	  coeff2 += 1.0l/j;       /* 1 + 1/2 + 1/3 + .. 1/j */
	}
      sofar -= I * logplus;
      
      return ( sofar );
    }
  /* That failing, use asymptotic series. */
  logplus = ( just_K0 ? 0 : log ( 0.5 * z ) + Gamma_E );
  I = sqrt(PI/(2*z)) * exp(-z); /* leading order asymptotic of K0 */
  z = 1.0l/z;
  sofar = ( ( ( -0.0732421875 * z + 0.0703125 ) * z - 0.125 ) * z + 1 );
  /* the numbers are 1, -1/8, 1*9/8*8*2!, -1*9*25/8*8*8*3! */
  sofar = I * sofar + logplus;
  return ( sofar );
}

double K_BDMPS ( double z )
     /* Evaluates the BDMPS approximation version of the above,
	that is, 0.5 * ( 1 - z K_1(z) ), K_1 the modified Bessel
	function. */
{
  double sofar , I , coeff1 , psi , the_log;
  int    j;

  if ( z < 0 ) return ( K_BDMPS ( -z ) );
  if ( z == 0 ) return ( 0 );
  if ( z < 9 )
    {
      /* Note, I will really be (z/2) I_0. */

      the_log = log ( 0.5 * z );
      z *= 0.25 * z; /* Series is in z^2 / 4. */
      coeff1 = z;
      I = 0; sofar = 0;
      psi = - 2 * Gamma_E + 1; /* psi(1) + psi(2) */
      j = 0;
      while ( coeff1 > 1.0e-13 )
	{
	  I += coeff1;
	  sofar += coeff1 * psi;
	  j++;
	  coeff1 *= z / ( j * (j+1) );
	  psi += 1.0l/j + 1.0l/(j+1);
	}
      sofar = 0.5 * sofar - the_log * I;
      /* This is (1/2) (1 - z K_1(z)) */
    }
  else
    {
      sofar = sqrt ( 0.125 * z * PI ) * exp(-z);
      I = 1;
      coeff1 = 1;
      z = 0.125 / z;
      for ( j = 1 ; j < 12 ; j++ )
	{
	  coeff1 *= (4 - (2*j-1)*(2*j-1)) * ( z / j );
	  I += coeff1;
	}
      sofar = 0.5 - I * sofar;
    }
  return ( sofar );
}

/* ------------------------------------------------------------- */

/* Now, function to evaluate derivative, for solving differential
   equation (for determining splitting kernel) by method of Runge-Kutta.  
   
   f[0] = Re f
   f[1] = Im f
   f[2] = Re f'
   f[3] = Im f'
   Therefore
   fpr[0] = f[2]
   fpr[1] = f[3]
   fpr[2] = -3/b f[2] + A f[0] + DK f[1]
   fpr[3] = -3/b f[3] + A f[1] - DK f[0] */

#define get_fpr(f,fpr,A,DK,b) {double b_coeff;\
        b_coeff=3/(b); fpr[0]=f[2];fpr[1]=f[3];\
        fpr[2]=-b_coeff*f[2]+A*f[0]+DK*f[1];\
        fpr[3]=-b_coeff*f[3]+A*f[1]-DK*f[0];}

/* -----------------------------------------------------------  */

void RK_step ( double f[4] , double A , double D , double *b , double *k ,  
  	       double C[3] , double r[3] , int BDMPS )
     /* Performs a 4'th order Runge-Kutta step.  
  	Internally determines appropriate step length. 
  	Note that I am stepping from larger to smaller b.   */
     /* C are the Casimirs for the three terms in the collision piece, 
  	r are the rescalings of b for those three terms.   */
{ 
  double fpr1[4] , fpr2[4] , fpr3[4] , fpr4[4] , ftmp[4]; 
  double step; 
  int    i; 
  
  /* First, determine how large a step is permitted.  */
  step = 3.0l / *b; 
  step = 0.03l / sqrt ( step * step + sqrt ( A*A + D*D* *k * *k ) ); 
  /*  The quantity in the square root is an estimate of the square of scale 
      of variation of the function f.  (scale)/30 is good enough with 
      Runge-Kutta to give 10^{-8} errors.  Tighter doesn't help (checked)  */
  
  /* Recall that RK update is: 
     Get fpr1 = f'[b,f[b]] 
     fpr2 = f'[b+step/2,f[b]+step fpr1 /2] 
     fpr3 = f'[b+step/2,f[b]+step fpr2 /2] 
     fpr4 = f'[b+step,f[b]+step fpr3] 
     fprime = (1/6) [fpr1 + fpr4 + 2 fpr2 + 2 fpr3] 
     The f[b]+... is an estimate of the value at b+step/2 or b+step 
     using the fprime's just derived.   */
  
  get_fpr ( f , fpr1 , A , D* *k , *b ); 
  for ( i = 0 ; i < 4 ; i++ ) ftmp[i] = f[i] - 0.5 * fpr1[i] * step; 
  /* fpr1 = slope at b; ftmp = first midpt extrapolation.  */
  *b -= 0.5 * step; 
  if ( BDMPS )
    *k = C[0] * K_BDMPS (*b * r[0] ) + C[1] * K_BDMPS (*b * r[1] ) 
      + C[2] * K_BDMPS (*b * r[2] ); 
  else
    *k = C[0] * K (*b * r[0] , 0 ) + C[1] * K (*b * r[1] , 0 ) 
      + C[2] * K (*b * r[2] , 0 ); 
  get_fpr ( ftmp , fpr2 , A , D* *k , *b ); 
  for ( i = 0 ; i < 4 ; i++ ) ftmp[i] = f[i] - 0.5 * fpr2[i] * step; 
  /* fpr2 = first guess slope at b-step/2; ftmp=second mid extrap.  */
  get_fpr ( ftmp , fpr3 , A , D* *k , *b ); 
  /* fpr3 = second guess slope at b-step/2.  */
  *b -= 0.5 * step; 
  if ( BDMPS )
    *k = C[0] * K_BDMPS (*b * r[0] ) + C[1] * K_BDMPS (*b * r[1] ) 
      + C[2] * K_BDMPS (*b * r[2] ); 
  else
    *k = C[0] * K (*b * r[0] , 0 ) + C[1] * K (*b * r[1] , 0 ) 
      + C[2] * K (*b * r[2] , 0 ); 
  for ( i = 0 ; i < 4 ; i++ ) ftmp[i] = f[i] - fpr3[i] * step; 
  get_fpr ( ftmp , fpr4 , A , D* *k , *b ); 
  /* ftmp=endpt extrap and fpr4 is guess slope there.  */
  for ( i = 0 ; i < 4 ; i++ ) 
    f[i] -= step * ( fpr1[i] + 2*fpr2[i] + 2*fpr3[i] + fpr4[i] ) * SIXTH; 
} 

double find_f_atzero ( double A , double D , double C[3] , double r[3] ,
		       int BDMPS , double kmaxsq )
     /* solves the differential equation with exponentially shrinking 
  	large value boundary conditions to find the finite piece 
  	orthogonal to the divergent piece at the origin. 
  	Here C are the Casimirs and r the rescalings for the 3 terms 
  	in the collision piece.   */
     /* Here kmax is the scale for cutting of the transverse momentum
        which means we use the envelope function 
	(kmax*kmax/2)*exp(-kmax*kmax*b*b/4)         */
{ 
  double f[4] , oldf[4] , b , b_start , KD , b_min , result; 
  double bess_coeff[4] , j , jpr , y , ypr , k; 
  double delta_b;
  double b_int[2] , envelope , envelope_pr , old_envelope , old_envelope_pr;

  /* First task:  what b is big enough?  Based on rate of exponential  
     falloff in the tail--ignore 1/b term and b dependence of K.   */
  b = 2;  
  do 
    { 
      b_start = b; 
      if ( BDMPS )
	KD = D*( C[0] * K_BDMPS (b * r[0] ) + C[1] * K_BDMPS (b * r[1] ) +  
		 C[2] * K_BDMPS (b * r[2] ) ); 
      else
	KD = D*( C[0] * K (b * r[0] , 0 ) + C[1] * K (b * r[1] , 0 ) +  
		 C[2] * K (b * r[2] , 0 ) ); 
      b = 28 / sqrt ( A + sqrt ( A*A + KD * KD ) ); 
      /* denominator = real part of negative eigenvalue /sqrt(2)  */
      b_start /= b; 
    } 
  while ( b_start < 0.95 || b_start > 1.05 ); 
  f[0] = 1; 
  f[1] = 0; 
  f[2] = -sqrt ( 0.5 * ( sqrt ( A*A + KD*KD ) + A ) ); 
  f[3] =  sqrt ( 0.5 * ( sqrt ( A*A + KD*KD ) - A ) ); 
  /* Approximate shrinking solution--assumed 
     1) you can neglect f'/b term 
     2) you can take K constant.  Both valid far enough out in tail.  */
  /* Now find how small b must go before I accurately see b~=0 behavior  */
  b_min = 0.005 / sqrt ( 1 + sqrt ( A*A+D*D ) ); 
  if ( BDMPS )
    k = C[0] * K_BDMPS (b * r[0] ) + C[1] * K_BDMPS (b * r[1] ) 
      + C[2] * K_BDMPS (b * r[2] ); 
  else
    k = C[0] * K (b * r[0] , 0 ) + C[1] * K (b * r[1] , 0 ) 
      + C[2] * K (b * r[2] , 0 ); 
  /* Now run differential equation from b to b_min  */
  b_int[0] = 0; b_int[1] = 0;
  envelope = b * kmaxsq * 0.5 * exp ( - kmaxsq * b * b * 0.25 );
  envelope_pr = ( kmaxsq * 0.5 * exp ( -kmaxsq * b * b * 0.25 )
		  * ( 1 - kmaxsq * b * b * 0.5  ) );
  while ( b > b_min ) 
    { 
      /* Perform Runge-Kutta.  But also evaluate the integral of
	 (f*envelope) bdb.  I want fourth-order accuracy but only
	 have info at boundaries; however
	 int_a^b g(x) dx = (b-a)(g(a)+g(b))/2
	                  -(b-a)^2 (g'(b)-g'(a))/12 + O((b-a)^5) */

      delta_b = b; /* Need to find out how much b is changing */
      oldf[0] = f[0]; oldf[1] = f[1]; oldf[2] = f[2]; oldf[3] = f[3]; 
      old_envelope = envelope; old_envelope_pr = envelope_pr;
      /* I need this because I need 
	 (f*envelope)' = f' * envelope + f * envelope' */
      RK_step ( f , A , D , &b , &k , C , r , BDMPS ); 
      delta_b -= b;
      envelope = b * kmaxsq * 0.5 * exp ( - kmaxsq * b * b * 0.25 );
      envelope_pr = ( kmaxsq * 0.5 * exp ( -kmaxsq * b * b * 0.25 )
		      * ( 1 - kmaxsq * b * b * 0.5  ) );
      b_int[0] += ( ( envelope * f[0] + old_envelope * oldf[0] ) -
		    ( old_envelope_pr * oldf[0] + old_envelope * oldf[2]
		      - envelope_pr * f[0] - envelope * f[2] )
		    * SIXTH * delta_b ) * 0.5 * delta_b;
      b_int[1] += ( ( envelope * f[1] + old_envelope * oldf[1] ) -
		    ( old_envelope_pr * oldf[1] + old_envelope * oldf[3]
		      - envelope_pr * f[1] - envelope * f[3] )
		    * SIXTH * delta_b ) * 0.5 * delta_b;



      if ( f[0]*f[0] + f[1]*f[1] > 100 ) /* Apply harmless rescaling  */
  	{ 
	  f[0] *= 0.125; f[1] *= 0.125; 
	  f[2] *= 0.125; f[3] *= 0.125; 
	  b_int[0] *= 0.125; b_int[1] *= 0.125;
	} 
    } 
  /* Find projection of solution onto divergent and finite Bessel 
     type (D=0) solutions.  Then find the finite part, orthogonal to the 
     divergent one.  bess_coeff[0] and [1] are the Re, Im parts of 
     divergent, and bess_coeff[2] and [3] of finite, solutions, 
     normalized to to 1/b^2 and 1.   */
  
  /* fourth-order finite and divergent solutions with derivatives
     Note that "j" = 2 j_1/b and "y" = y_1/b plus some (A dependent) 
     multiple of j to make the logs simpler.  This doesn't matter 
     since I want "j" in the direction orthogonal to where "y" != 0.  */
  j = 1 + 0.125 * A*b*b * ( 1 + .04166666666*A*b*b ); 
  jpr = 0.25 * A * b * (1+.083333333333*A*b*b); 
  y = 1.0/(b*b) + 0.5 * A * ( log(b) - 0.5 + b*b*A*0.125* 
  			      (log(b)-1.25) ); 
  ypr = -2.0/(b*b*b) + 0.5 * A / b 
    + 0.125*A*A*b*(log(b)-0.75); 
  
  bess_coeff[0] = 0.5 * b * b * b * ( jpr * f[0] - j * f[2] ); 
  bess_coeff[1] = 0.5 * b * b * b * ( jpr * f[1] - j * f[3] ); 
  bess_coeff[2] = 0.5 * b * b * b * ( -ypr* f[0] + y * f[2] ); 
  bess_coeff[3] = 0.5 * b * b * b * ( -ypr* f[1] + y * f[3] ); 
  /* The 0.5 b^3 is a Wronskian for the solutions . . .   */
  result = (bess_coeff[0]*b_int[1]-bess_coeff[1]*b_int[0]) / 
    (bess_coeff[0]*bess_coeff[0]+bess_coeff[1]*bess_coeff[1]); 
  b = kmaxsq * b * b * 0.25; /* Argument of exponential */
  A /= 2 * kmaxsq;
  result +=  (bess_coeff[0]*bess_coeff[3]-bess_coeff[1]*bess_coeff[2]) / 
    (bess_coeff[0]*bess_coeff[0]+bess_coeff[1]*bess_coeff[1])
    * ( 1 + A - ( 1 + A + A * b ) * exp ( - b ) );
  /* Integral from old b down to 0:  the first two lines are the
     coefficient of the constant term; the second line is
     int_0^b  b' db' exp(-kmaxsq b'^2/4) (1+Ab^2/8)    A,b=original A,b
     = int_0^b dx exp(-x) (1+Ax)    A,b=new A,b defined above   */
  
  return ( result ); 
} 

//JF: p and k are in units of temperature in this function
double solve_int_eqn ( double p , double k ,  
  		       double df , double da , double cf , double ca , 
  		       int  n_quark_flavors , int p_is_gluon , 
		       int all_are_gluons , int emit_photon , 
		       int Bethe_Heitler ,
		       int BDMPS , double mDsq ) 
     /* Solves the integral equations to give the production  
  	rate at p,k.  Contains all the annoying overall factors,
        except a factor of Nf in g->qq.  */
{ 
  double A , B , D , sofar , C[3] , r[3]; 
  double md2_over_g2T2 , kappa , kmaxsq;
  
  if ( k*k < 1.0e-9 || (k-p)*(k-p) < 1.0e-9 ) 
    { 
      /* Do average on either side of p, to avoid p=0 complications 
  	 in the following procedure.  */
      sofar = (solve_int_eqn( p , k-.001 , df , da , cf , ca , 
			      n_quark_flavors , p_is_gluon , all_are_gluons
			      , emit_photon , Bethe_Heitler , 
			      BDMPS , mDsq )
	       +solve_int_eqn ( p , k+.001 , df , da , cf , ca , 
				n_quark_flavors,p_is_gluon , all_are_gluons
				, emit_photon , Bethe_Heitler , 
				BDMPS , mDsq )
		) * 0.5;
      return ( sofar );
    }
  md2_over_g2T2 = ( ca + n_quark_flavors * cf * df / da ) / 3.0l;
  kappa = 0.25 * cf / md2_over_g2T2;
  if ( p_is_gluon ) 
    { /* Then gluon goes either to 2 quarks or 2 gluons. */
      A = -1 / ( 4 * p );
      if ( all_are_gluons )
	A += 1 / ( 4.0 * k ) + 1 / ( 4.0 * (p-k) );
      else
	A += kappa * ( 1 / ( 2.0 * k ) + 1 / ( 2.0 * (p-k) ) );
    }
  else
    { /* p and (p-k) are quarks and k is a gluon. */
      A = kappa * ( 1 / ( 2.0 * ( p - k ) ) - 1 / ( 2 * p ) );
      if ( ! emit_photon )  /* Photon is approximately massless */
	A += 1 / ( 4 * k ); /* Gluon mass contribution. */
    }
  A *= md2_over_g2T2;
  B = md2_over_g2T2 * p / ( 2 * k * (p-k) );
  D = 1.0 / ( 2 * PI );
  
  /* Respectively, const and p^2 coefficients of delta E, and
     coefficient of collision term. */

  A /= B;
  D /= B;
  if ( BDMPS ) 
  A = 0; /* Because BDMPS forgot about thermal masses. */

  r[0] = 1;
  r[1] = ( ( k > 0 ) ? k / p : -k / p );
  r[2] = ( (p-k) > 0 ? (p-k)/p : (k-p)/p );
  if ( ABS(p) > ABS(k) )
    if ( ABS(p-k) > ABS(k) )
      kmaxsq = k * k / mDsq;
    else
      kmaxsq = (p-k) * (p-k) / mDsq;
  else 
    if ( ABS(p) > ABS(p-k) )
      kmaxsq = (p-k) * (p-k) / mDsq;
    else
      kmaxsq = p * p / mDsq;
  if ( p_is_gluon )
    {
      if ( all_are_gluons )
	{
	  C[0] = ca / 2;
	  C[1] = ca / 2;
	  C[2] = ca / 2;
	}
      else
	{
	  C[0] = cf - ca / 2;
	  C[1] = ca / 2;
	  C[2] = ca / 2;
	}
    }
  else if ( emit_photon )
    {
      C[0] = 0;
      C[1] = cf;
      C[2] = 0;
    }
  else
    {
      C[0] = ca/2;
      C[1] = cf - ca/2;
      C[2] = ca/2;
    }
  kmaxsq=1.0e30;
  /* Solving everything inside the p_parallel integration. */
  if ( Bethe_Heitler )
    {
      /* In this case we want to work perturbatively in D.  Do this
	 by determining things at D= (correct D)/100 and (correct D)/200
	 and use these results to extrapolate to zero D. */
      /* This is the "bad way."  I should fix this some day. */

      sofar = 400  * find_f_atzero ( A , D/200 , C , r , 
				     BDMPS , kmaxsq );
      sofar -= 100 * find_f_atzero ( A , D/100 , C , r , 
				     BDMPS , kmaxsq );
      sofar *= 4 / ( PI * B );
    }
  else
    {
      sofar = 4 * find_f_atzero ( A , D , C , r , 
				  BDMPS , kmaxsq ) / ( PI * B );
    }  
  sofar *= p*p / (16 * PI ) * md2_over_g2T2 * md2_over_g2T2;
  if ( p_is_gluon )
    {
      if ( all_are_gluons )
	{
	  sofar *= ca * ( p*p*p*p + k*k*k*k + (p-k)*(p-k)*(p-k)*(p-k) )
	    / ( p*p*p*k*k*k*(p-k)*(p-k)*(p-k) )
	    * BoseStim(k) * BoseStim((p-k));
	}
      else
	{
	  sofar *= cf * ( k*k + (p-k)*(p-k) ) / ( k*k*(p-k)*(p-k)*p*p*p )
	    * PauliBlock(k) * PauliBlock((p-k));
	}
    }
  else
    {
      sofar *= ( p*p + (p-k)*(p-k) ) / ( p*p*(p-k)*(p-k)*k*k*k )
	* PauliBlock((p-k));
      if ( emit_photon ) 
	{
	  if ( k < 0 )
	    sofar = 0; /* Because there are no photons to "pick up". */
	}
      else
	{
	  sofar *= cf * BoseStim(k);
	}
    }
  return ( sofar );
  /* Rate per unit k, except for g^4 T factor (and Nf for g->qqbar)! */
}

/* ================================================================= */
/* Conclusion of code to evaluate dGamma(p,k)/dk dx and beginning of */
/* code for applying this to the evolution problem.                  */
/* ================================================================= */

/* Now, equipment for building a table for interpolation, and for
   storing, reading, and using that table. */

/* NORMALIZATION:

   In
     build_table
     write_table
     read_table
     use_table
     prepare_table
     prep_dGamma
     and in the stored tables
   The normalization is, that the dGamma(p,k)/dk dx given, is
   in units of 1/g^4 except for the photon one, which is in units
   of 1/g^2 e^2.  That means, that to convert dGamma(p,k)/dk dx
   to the rate of emission per momentum range $k$ and distance or
   time $x$, one should multiply by g^4 or g^2 e^2.
   In other words, the dGamma, integrated \int dk/T, gives the
   likelihood to make an emission in a length 1/(g^4 T) of plasma.

   In
     find_dP
     evolve_one_step
   the normalization of dx or dt is, that it is 1/g^4 times an inverse
   energy (in GeV).

   It is in evolve_in_medium that the conversion to Fermis is performed,
   meaning the inclusion of the 1/g^4 and of the relation between a
   Fermi and a GeV^-1.  
*/

//JF: pretty sure p and k are in units of temperature in this function
void build_table ( Gamma_info * dat , dGammas *Gam )
     /* Determines the ranges of p,k and spacing dp to consider,
	and loads up an array of values of dGamma/dkdx needed for 
	updating. */
{
  double k , p , rate , mDsq;
  int    i_p , i_k , i_g;
//  char   c;
 
  //Assume Nc is already set
  //printf ( "Assuming number of colors is 3.  To change, you have to\n" );
  //printf ( "actually edit the code.\n" );
  //dat->Nc = 3;
  dat->df = dat->Nc;
  dat->da = dat->Nc * dat->Nc - 1;
  dat->cf = dat->da / ( 2 * dat->df );
  dat->ca = dat->Nc;
  //Assumes Nf is already set in "dat"
  //printf ( "Enter the number of flavors.\n" );
  //scanf ( "%d" , &(dat->Nf) ); READ_TO_EOLN ( stdin );
  //printf ( "Enter b for Bethe-Heitler, any other letter for full LPM.\n" );
  //READ_LETTER ( c , stdin );
  //dat->Bethe_Heitler = ( c == 'b' ? 1 : 0 );
  dat->Bethe_Heitler = 0; //No Bethe-Heitler approx by default
  //printf ( "Enter b for BDMPS approximations, other letters for " );
  //printf ( "the correct treatment.\n" );
  //READ_LETTER ( c , stdin );
  //dat->BDMPS = ( c == 'b' ? 1 : 0 );
  dat->BDMPS = 0; //No BDMPS by default
  dat->include_gluons = 1;
  /* Gluon inclusion will be the default, since we are including
     photons as well! */
  //printf ( "Enter alpha_s (needed for k_perp k relation!)\n" );
  //scanf ( "%lf" , &mDsq );
  mDsq=dat->alpha_s;
  mDsq *= 4 * PI * ( 2 * dat->Nc +dat->Nf ) / 6.0l;
  /* mD^2 in units of temperature. */

  for ( i_p = 0 ; i_p < NP ; i_p++ )
    for ( i_k = 0 ; i_k < NK ; i_k++ )
      {
	Gam->dGamma[i_p][i_k] = 0; 
	Gam->dGamma_gqq[i_p][i_k] = 0;
	Gam->dGamma_ggg[i_p][i_k] = 0; /* Initialize them all to zero. */
	Gam->dGamma_em[i_p][i_k] = 0;
      }
  for ( i_g = 0 ; i_g < 4 ; i_g++ )
    for ( i_p = 0 ; i_p < NP ; i_p++ )
      {
	if ( i_p == 0 )
	  p = 4.01;
	else
	  p = p * 1.04119; /* spaced so 6---1000 is 0--127 */
	printf ( "Starting p %e\n" , p );
	for ( i_k = 0 ; i_k < ( ( (i_g % 3) == 0 ) ? NK : NK - 160 ) ; i_k++ )
	  {
	    if ( i_k < 50 )        /* spaced by 0.2  from -12 to -2 */
	      k = -12 + i_k * 0.2;
	    else if ( i_k < 60 )   /* spaced by 0.1  from -2  to -1 */
	      k = -2 + (i_k-50) * 0.1;
	    else if ( i_k < 100 )  /* spaced by 0.05 from -1  to +1 */
	      k = -1 + (i_k-60) * 0.05;
	    else if ( i_k < 110 )  /* spaced by 0.1  from +1  to +2 */
	      k = 1 + (i_k-100) * 0.1;
	    else if ( i_k < 270 )  /* spaced complicated, +2 to p-2 */
	      {
		k = 0.1 * (i_k-190);
		k = 2 + (p-4) * ( -0.0003353501304664781l
				  + 1.000670700260932956l / (1+exp(-k)) );
	      }
	    else if ( i_k < 280 )  /* spaced by 0.1  from p-2 to p-1 */
	      k = p - 2 + 0.1 * (i_k-270);
	    else if ( i_k < 320 )  /* spaced by 0.05 from p-1 to p+1 */
	      k = p + 0.05 * (i_k - 300);
	    else if ( i_k < 330 )  /* spaced by 0.1  from p+1 to p+2 */
	      k = p + 0.1 * (i_k - 310);
	    else                   /* spaced by 0.2  from p+2 to p+12 */
	      k = p + 0.2 * (i_k - 320);
	    
	    if ( ABS(k) > 0.001 )
	      {
		rate = 
		  solve_int_eqn ( p , k , dat->df , dat->da , dat->cf , 
				  dat->ca , dat->Nf , i_g%3 , (i_g == 2) , 
				  (i_g == 3) , dat->Bethe_Heitler , 
				  dat->BDMPS , mDsq );
		switch ( i_g )
		  {
		  case 0:
		    rate *=  k;
		    if ( k < 20 ) rate *= 1 - exp(-k);
		    if ( k > p - 20 ) rate *= 1 + exp(k-p);
		    Gam->dGamma[i_p][i_k] = rate;
		    break;
		  case 1:
		    rate *= p;
		    if ( k < 20 ) rate *= 1 + exp(-k);
		    if ( k > p - 20 ) rate *= 1 + exp(k-p);
		    Gam->dGamma_gqq[i_p][i_k] = rate;
		    break;
		  case 2:
		    rate *= k * (p-k) / p;
		    if ( k < 20 ) rate *= 1 - exp(-k);
		    if ( k > p - 20 ) rate *= 1 - exp(k-p);
		    Gam->dGamma_ggg[i_p][i_k] = rate;
		    break;
		  case 3:
		    rate *= k;
		    if ( k < 0 ) rate = 0; /* No photon pickup allowed */
		    if ( k > p-20 ) rate *= 1 + exp(k-p);
		    Gam->dGamma_em[i_p][i_k] = rate;
		    break;
		  }
	      }
	  }
	switch ( i_g )
	  {
	  case (0):
	    Gam->dGamma[i_p][80] = 0.5 * ( Gam->dGamma[i_p][79] 
					   + Gam->dGamma[i_p][81] );
	    break;
	  case(1):
	    Gam->dGamma_gqq[i_p][80] = 0.5 * ( Gam->dGamma_gqq[i_p][79] 
					       + Gam->dGamma_gqq[i_p][81] );
	    break;
	  case(2):
	    Gam->dGamma_ggg[i_p][80] = 0.5 * ( Gam->dGamma_ggg[i_p][79] 
					       + Gam->dGamma_ggg[i_p][81] );
	  case(3):
	    Gam->dGamma_em[i_p][80] = 0.5 * ( Gam->dGamma_em[i_p][79] 
					      + Gam->dGamma_em[i_p][81] );
	  }
      }
}

void write_table ( Gamma_info *dat , dGammas * Gam )
     /* Writes the table, in binary, to a file. */
{
  FILE *wfile;
//  int  i;
  printf ("df=%g\n", (dat->df)); 
  printf ("da=%g\n", (dat->da)); 
  printf ("cf=%g\n", (dat->cf)); 
  printf ("ca=%g\n", (dat->ca)); 
  printf ("nc=%d\n", (dat->Nc)); 
  printf ("nf=%d\n", (dat->Nf)); 
  printf ("bh=%d\n", (dat->Bethe_Heitler));
  printf ("bdmps=%d\n", (dat->BDMPS) );
  printf ("includegluons=%d\n", (dat->include_gluons) );
  printf ("dGggg=%g\n", Gam->dGamma_ggg[0][0] );

  //printf ( "Enter the file name to write Gamma information into.\n" );
  //wfile = openfile ( dat->in_fname , "w" );
  wfile=fopen( dat->in_fname, "wb" );
  if (wfile==NULL) printf("can't open file!!\n");
  fwrite ( (&dat->df) , sizeof ( double ) , 1 , wfile );
  fwrite ( (&dat->da) , sizeof ( double ) , 1 , wfile );
  fwrite ( (&dat->cf) , sizeof ( double ) , 1 , wfile );
  fwrite ( (&dat->ca) , sizeof ( double ) , 1 , wfile );
  fwrite ( (&dat->Nc) , sizeof ( int ) , 1 , wfile );
  fwrite ( (&dat->Nf) , sizeof ( int ) , 1 , wfile );
  fwrite ( (&dat->Bethe_Heitler) , sizeof ( int ) , 1 , wfile );
  fwrite ( (&dat->BDMPS) , sizeof ( int ) , 1 , wfile );
  fwrite ( (&dat->include_gluons) , sizeof ( int ) , 1 , wfile );
  fwrite ( Gam->dGamma , sizeof ( double ) , NP * NK , wfile );
  fwrite ( Gam->dGamma_gqq , sizeof ( double ) , NP * NK , wfile );
  fwrite ( Gam->dGamma_ggg , sizeof ( double ) , NP * NK , wfile );
  fwrite ( Gam->dGamma_em , sizeof ( double ) , NP * NK , wfile );

  //fprintf(wfile,"DF=%g \n",dat->df);
  //fprintf(wfile,"CF=%g \n",dat->cf);
  //fprintf(wfile,"CA=%g \n",dat->ca);
  //fprintf(wfile,"Nc=%d \n",dat->Nc);
  //fprintf(wfile,"Nf=%d \n",dat->Nf);
  //fprintf(wfile,"Bether-Heither=%d \n",dat->Bethe_Heitler);
  //fprintf(wfile,"BDMPS=%d \n",dat->BDMPS);
  //fprintf(wfile,"Include gluons?=%d \n",dat->include_gluons);
  //fprintf(wfile,"p/T k/T dGamma_ggg\n");
//  double p=4.01;
//  double k;
//  for(int ip=0;ip<NP;ip++) {
//	  for(int i_k=0;i_k<NK;i_k++) {
//
//		  if ( i_k < 50 )        /* spaced by 0.2  from -12 to -2 */
//			  k = -12 + i_k * 0.2;
//		  else if ( i_k < 60 )   /* spaced by 0.1  from -2  to -1 */
//			  k = -2 + (i_k-50) * 0.1;
//		  else if ( i_k < 100 )  /* spaced by 0.05 from -1  to +1 */
//			  k = -1 + (i_k-60) * 0.05;
//		  else if ( i_k < 110 )  /* spaced by 0.1  from +1  to +2 */
//			  k = 1 + (i_k-100) * 0.1;
//		  else if ( i_k < 270 )  /* spaced complicated, +2 to p-2 */
//		  {
//			  k = 0.1 * (i_k-190);
//			  k = 2 + (p-4) * ( -0.0003353501304664781l
//					  + 1.000670700260932956l / (1+exp(-k)) );
//		  }
//		  else if ( i_k < 280 )  /* spaced by 0.1  from p-2 to p-1 */
//			  k = p - 2 + 0.1 * (i_k-270);
//		  else if ( i_k < 320 )  /* spaced by 0.05 from p-1 to p+1 */
//			  k = p + 0.05 * (i_k - 300);
//		  else if ( i_k < 330 )  /* spaced by 0.1  from p+1 to p+2 */
//			  k = p + 0.1 * (i_k - 310);
//		  else                   /* spaced by 0.2  from p+2 to p+12 */
//			  k = p + 0.2 * (i_k - 320);
//
//	  	  //getting the value of p and k is non-trivial!
//		  fprintf(wfile,"%.17g %.17g %g \n",p,k,Gam.dGamma_ggg[ip][i_k]);
//	}
//	p = p * 1.04119; /* spaced so 6---1000 is 0--127 */
//  }
  //fwrite ( (char *) Gam.dGamma_ggg , sizeof ( double ) , NP * NK , wfile );

  //fwrite ( (char *) (&dat->da) , sizeof ( double ) , 1 , wfile );
  //fwrite ( (char *) (&dat->cf) , sizeof ( double ) , 1 , wfile );
  //fwrite ( (char *) (&dat->ca) , sizeof ( double ) , 1 , wfile );
  //fwrite ( (char *) (&dat->Nc) , sizeof ( int ) , 1 , wfile );
  //fwrite ( (char *) (&dat->Nf) , sizeof ( int ) , 1 , wfile );
  //fwrite ( (char *) (&dat->Bethe_Heitler) , sizeof ( int ) , 1 , wfile );
  //fwrite ( (char *) (&dat->BDMPS) , sizeof ( int ) , 1 , wfile );
  //fwrite ( (char *) (&dat->include_gluons) , sizeof ( int ) , 1 , wfile );
  //fwrite ( (char *) Gam.dGamma , sizeof ( double ) , NP * NK , wfile );
  //fwrite ( (char *) Gam.dGamma_gqq , sizeof ( double ) , NP * NK , wfile );
  //fwrite ( (char *) Gam.dGamma_ggg , sizeof ( double ) , NP * NK , wfile );
  //fwrite ( (char *) Gam.dGamma_em , sizeof ( double ) , NP * NK , wfile );
  fclose ( wfile );
}

void read_table ( Gamma_info *dat , dGammas *Gam )
     /* Reads in the binary stored file of dGamma values. */
{
  FILE *rfile;
  
  rfile = fopen ( dat->in_fname , "rb" );
  if (rfile == NULL) printf("can't read file '%s', can't open it\n", dat->in_fname);
  fread (  (&dat->df) , sizeof ( double ) , 1 , rfile );
  fread (  (&dat->da) , sizeof ( double ) , 1 , rfile );
  fread (  (&dat->cf) , sizeof ( double ) , 1 , rfile );
  fread (  (&dat->ca) , sizeof ( double ) , 1 , rfile );
  fread (  (&dat->Nc) , sizeof ( int ) , 1 , rfile );
  fread (  (&dat->Nf) , sizeof ( int ) , 1 , rfile );
  fread (  (&dat->Bethe_Heitler) , sizeof ( int ) , 1 , rfile );
  fread (  (&dat->BDMPS) , sizeof ( int ) , 1 , rfile );
  fread (  (&dat->include_gluons) , sizeof ( int ) , 1 , rfile );
  fread (  Gam->dGamma , sizeof ( double ) , NP * NK , rfile );
  fread (  Gam->dGamma_gqq , sizeof ( double ) , NP * NK , rfile );
  fread (  Gam->dGamma_ggg , sizeof ( double ) , NP * NK , rfile );
  fread (  Gam->dGamma_em , sizeof ( double ) , NP * NK , rfile );
  fclose ( rfile );
}

double use_table ( double p , double k , double dGamma[NP][NK] , 
		   int which_kind )
     /* Uses the lookup table and simple interpolation to get the value
	of dGamma/dk dx at some value of p,k. */
     /* This works by inverting the relations between (p,k) and (n_p,n_k)
	used in building the table, to find out what continuous values
	of n_p, n_k should be considered; then linearly interpolates. */
{
  double a , b , result;     /* fraction of way from corner of box. */
  int    n_p , n_k; /* location of corner of box. */

  if ( (p<4.01) | (p>46000) | (k<-12) | (k>p+12) ) return ( 0 );
  /* Out of range. */
  if ( ( which_kind % 3 ) && ( k > p/2 ) )
    k = p - k;  /* Take advantage of symmetry in these cases */
  a = 24.7743737154026 * log ( p * .2493765586034912718l );
  n_p = (int) a;
  a -= n_p;
  if ( k < 2 )
    {
      if ( k < -1 )
	{
	  if ( k < -2 )
	    b = 60 + 5*k;
	  else
	    b = 70+10*k;
	}
      else
	{
	  if ( k < 1 )
	    b = 80 + 20*k;
	  else
	    b = 90 + 10*k;
	}
    }
  else if ( k < p-2 )
    { /* This is that tricky middle ground. */
      b = 190 - 10*log ( 1.000670700260932956l / 
			 ( 0.0003353501304664781l + (k-2) / (p-4) ) - 1 );
    }
  else
    {
      if ( k < p+1 )
	{
	  if ( k < p-1 )
	    b = 290 + 10*(k-p);
	  else
	    b = 300 + 20*(k-p);
	}
      else
	{
	  if ( k < p+2 )
	    b = 310 + 10*(k-p);
	  else
	    b = 320 + 5*(k-p);
	}
    }
  n_k = (int) b;
  b -= n_k;
  result = (1-a) * ( (1-b) * dGamma[n_p][n_k] + b * dGamma[n_p][n_k+1] )
    +        a * ( (1-b) * dGamma[n_p+1][n_k] + b * dGamma[n_p+1][n_k+1] );
  if ( ABS(k) > 0.001 ) /* Avoid division by 0, should never get asked for */
    {
      switch ( which_kind )
	{
	case 0:
	  result /= k;
	  if ( k < 20 )
	    result /= 1 - exp(-k);
	  if ( k > p - 20 )
	    result /= 1 + exp(k-p);
	  break;
	case 1:
	  result /= p;
	  if ( k < 20 )
	    result /= 1 + exp(-k);
	  if ( k > p - 20 )
	    result /= 1 + exp(k-p);
	  break;
	case 2:
	  result /= k * (p-k) / p;
	  if ( k < 20 )
	    result /= 1 - exp(-k);
	  if ( k > p - 20 )
	    result /= 1 - exp(k-p);
	  break;
	case 3:
	  result /= k;
	  if ( k < 0 ) result = 0;
	  if ( k > p-20 )
	    result /= 1 + exp(k-p);
	  break;
	}
    }
  return ( result );
}

//void prepare_table ( Gamma_info *dat , dGammas *Gam )
//     /* Either reads in or generates gamma table.  Simple driver. */
//{
//  char  c;
//  FILE  *infile; 
//
//  printf ( "If a data file with a table exists, enter y, otherwise n.\n" );
//  READ_LETTER ( c , stdin );
//  if ( c == 'y' )
//    {
//      printf ( "Enter the input file name.\n" );
//      infile = openfile ( dat->in_fname , "r" );
//      fclose ( infile );
//      read_table ( dat , Gam );
//    }
//  else
//    {
//      build_table ( dat , Gam );
//      write_table ( dat , *Gam );
//    }
//}
//
//
///* Here is the one which learns everything from the user or a file */
//
//void prep_equipment ( Gamma_info *dat , dGammas *Gam , 
//		      double ** dGam1 , double ** dGam2 , 
//		      double ** dGam3 , double ** dGam4 ,
//		      double ** P , double ** Pg , double ** Pem )
//     /* Learns the basic information needed in all that follows . . . . */
//     /* When that info is in a file, its storage order needs to be,
//	
//	Name of file with dGamma info
//	dp spacing of data
//	Maximum p
//	Minimum p
//	dx
//	alpha_s
//	alpha_EM
//	1 or 0:  1 include 22 photon, 0 do not
//	quark_num gluon_num photon_num (iterated over the number of p's)
//     */
//
//{
//  double junk;
//  char   c , in_fname[100];
//  FILE   *infile , *table_file;
//  int    i;
//
//  printf ( "If all information exists in a file, enter y, otherwise n.\n" );
//  READ_LETTER ( c , stdin );
//  if ( c == 'y' )
//    {
//      printf ( "Enter the name of the file with the starting " );
//      printf ( "configuration.\n" );
//      infile = openfile ( in_fname , "r" );
//      /* Now, get the name of the file with the dGamma info, from infile. */
//      i = 0;
//      while ( ( dat->in_fname[i++] = getc(infile) ) != '\n' );
//      dat->in_fname[--i] = '\0';
//      printf ( "About to tell you the file name.\n" );
//      printf ( "I find that the file name to read is, %s\n" , dat->in_fname );
//      table_file = fopen ( dat->in_fname , "r" );
//      fclose ( table_file );
//      read_table ( dat , Gam );
//      /* Got the dGamma info.  Next, get other info about the starting
//	 configuration. */
//      fscanf ( infile , "%lf" , &(dat->dp) ); READ_TO_EOLN ( stdin );
//      fscanf ( infile , "%lf" , &(dat->p_max) ); READ_TO_EOLN ( stdin );
//      fscanf ( infile , "%lf" , &(dat->p_min) ); READ_TO_EOLN ( stdin );
//    }
//  else
//    {
//      prepare_table ( dat , Gam );
//      printf ( "Enter the following energies in GeV.\n" );
//      printf ( "First, the discretization spacing of momenta.\n" );
//      scanf ( "%lf" , &(dat->dp) ); READ_TO_EOLN ( stdin );
//      printf ( "Enter the maximum, and then the minimum, p to consider.\n" );
//      printf ( "Remember that probability which exceeds p_max or goes\n" );
//      printf ( "under p_min, gets thrown away; " );
//      printf ( "results are unreliable within\n" );
//      printf ( "a couple of T of each limit.  ENTER 1 PER LINE!!\n" );
//      scanf ( "%lf" , &(dat->p_max) ); READ_TO_EOLN ( stdin );
//      scanf ( "%lf" , &(dat->p_min) ); READ_TO_EOLN ( stdin );
//      if ( dat->p_max < dat->p_min )
//	{
//	  printf ( "You entered minimum then maximum! Fixing...\n" );
//	  junk = dat->p_max;
//	  dat->p_max = dat->p_min;
//	  dat->p_min = junk;
//	}
//    }
//  dat->n_p =  (int) (0.4 + dat->p_max / dat->dp );
//  dat->p_max = dat->dp * dat->n_p;
//  dat->n_pmin = (int) (0.4 + dat->p_min / dat->dp );
//  dat->n_p -= dat->n_pmin - 1;
//  dat->p_min = dat->dp * dat->n_pmin;
//
//  dat->n_kmin = 1 + 2 * ( (int) (2.0 / dat->dp ) );
//  dat->k_min = -dat->dp * dat->n_kmin;
//  dat->n_k = (int) ( ( 8 + dat->p_max ) / ( 2 * dat->dp ) );
//  dat->k_max = 2 * dat->dp * ( dat->n_k - 1 ) + dat->k_min;
//  /* Arranges the k's to be in the range from about -4 GeV
//     to about p_max + 4 GeV.  For T=400 MeV, this is enough that
//     population functions are only cut off at 10 T. */
//  *dGam1 = (double *) malloc ( sizeof(double) * dat->n_p * dat->n_k );
//  *dGam2 = (double *) malloc ( sizeof(double) * dat->n_p * dat->n_k );
//  *dGam3 = (double *) malloc ( sizeof(double) * dat->n_p * dat->n_k );
//  *dGam4 = (double *) malloc ( sizeof(double) * dat->n_p * dat->n_k );
//  *P = (double *) malloc ( sizeof(double) * dat->n_p );
//  *Pg = (double *) malloc ( sizeof(double) * dat->n_p );
//  *Pem = (double *) malloc ( sizeof(double) * dat->n_p );
//  /* Allocate space for the required records. */
//  if ( c == 'y' )
//    {
//      fscanf ( infile , "%lf" , &(dat->delta_x) );
//      fscanf ( infile , "%lf" , &(dat->alpha_s) );
//      fscanf ( infile , "%lf" , &(dat->alpha) );
//      fscanf ( infile , "%d" , &(dat->photon_22) );
//      fscanf ( infile , "%d" , &(dat->do_collisional) );
//      fscanf ( infile , "%d" , &(dat->collisional_only) );
//      fscanf ( infile , "%d" , &(dat->do_compton) );
//      fscanf ( infile , "%d" , &(dat->do_comptong) );
//      for ( i = 0 ; i < dat->n_p ; i++ )
//	{
//	  fscanf ( infile , "%lf" , *P + i );
//	  fscanf ( infile , "%lf" , *Pg + i );
//	  fscanf ( infile , "%lf" , *Pem + i );
//	}
//    }
//  else
//    {
//      printf ( "Enter the safety margin in the algorithm to update in x.\n" );
//      printf ( "(This should be a number in the neighborhood of 0.1.)\n" );
//      scanf ( "%lf" , &(dat->delta_x) ); READ_TO_EOLN ( stdin );
//      printf ( "Enter the desired value of alpha_strong.\n" );
//      scanf ( "%lf" , &(dat->alpha_s) ); READ_TO_EOLN ( stdin );
//      printf ( "Enter the desired value of alpha.\n" );
//      scanf ( "%lf" , &(dat->alpha) ); READ_TO_EOLN ( stdin );
//      printf ( "Enter 1 to include 2<->2 photon production, 0 to" );
//      printf ( "leave it out.\n" );
//      scanf ( "%d" , &(dat->photon_22) ); READ_TO_EOLN ( stdin );
//      printf ( "Enter 1 to include collisional E-loss, 0 to leave it out\n" );
//      scanf ( "%d" , &(dat->do_collisional) ); READ_TO_EOLN ( stdin );
//      if ( dat->do_collisional )
//	{
//	  printf ( "Enter 1 for collisional ONLY, 0 to do both.\n" );
//	  scanf ( "%d" , &(dat->collisional_only) ); READ_TO_EOLN ( stdin );
//	}
//      else
//	dat->collisional_only = 0;
//      printf ( "Enter 1 to include Compton q->g, 0 otherwise.\n" );
//      scanf ( "%d" , &(dat->do_compton) ); READ_TO_EOLN ( stdin );
//      printf ( "Same for Compton g->q processes.\n" );
//      scanf ( "%d" , &(dat->do_comptong) ); READ_TO_EOLN ( stdin );
//      printf ( "Enter p for polynomial, any other letter\n" );
//      printf ( "for a single initial value initial condition.\n" );
//      READ_LETTER ( c , stdin );
////      if ( c == 'p' )
////	initialize_power_law ( *P , *Pg , *Pem , *dat );
////      else
////	initialize_with_one_P ( *P , *Pg , *Pem , *dat );
//    }
//}
//
//void prep_dGamma ( Gamma_info *dat , dGammas *gam , 
//		    double T , double beta , double cos_phi , 
//		    double *dGamma , double *dGamma_gqq , 
//		   double *dGamma_ggg , double *dGamma_em )
//     /* Assumes that the Gamma_info and dGammas have been gotten already. */
//     /* Evaluates the rates of emission of hard collinear particles,
//	when in a medium of temperature T, moving with velocity beta
//	WRT the lab frame, for a hard parton moving at angle phi WRT
//	the boost direction of the medium.
//
//	For applications, the boost direction is the radial direction,
//	and phi is the azimuthal angle (cylindrical coordinates) of
//	the hard parton.
//
//	The basic formula is that the rest-frame energy is given by
//
//	p_rest = p_lab (1 - beta * cos_phi) / sqrt( 1 - beta * beta )
//	call 1/sqrt(1-beta*beta) = gamma.
//
//	so we want to know
//
//	dGamma[k_lab * gamma * (1-beta * cos_phi)] * (1-beta * cos_phi)
//
//	where the last bit is a Jacobian for dkdt between frames. 
//
//	When the dat.photon_22 is true, it adds the contribution
//	of Compton and pair processes under the approximation that
//	the photon momentum is the same as the particle's momentum
//	(which is valid at leading-log and receives O(T) corrections)
//	and that the particle emitting the photon is therefore
//	lost.  We do this by modifying dGamma(p,k=p) or dGamma(p,k=p+-Delta)
//	depending on whether p is an even or odd multiple of Delta.
//	We use the Next-to-Leading-Log calculation here.
//	The same procedure is used for purely QCD Compton/annihilation
//	processes.  These can be turned off by removing the
//	flag "DO_COMPTON".
//
//	when the do_collisional flag is true (1), it also adds a
//	"diffusion" verison of collisional energy loss.
//*/
//{
//  double e_scale , jacobian , gamma , tmp1 , tmp2;
//  double tot_dGamma , tot_dGammag , p , k , compton;
//  double dEdx_q , dEdx_g;
//  double mg2 , gs2 , nB , delta_kp;
//  int    i_p , i_k , here;
//
//  gamma = 1.0 / sqrt ( 1 - beta * beta );
//  jacobian = 1 - beta * cos_phi;
//  e_scale = gamma * jacobian / T;
//  delta_kp = dat->dp * e_scale;
//  nB = 1.0 / ( exp(delta_kp) - 1.0 );
//  gs2 = 4.0 * PI * dat->alpha_s;
//  mg2 = 0.5 * gs2 * ( 1.0 + dat->Nf / 6.0 );
//  
//  dat->dx_max = 0;
//  for ( i_p = 0 ; i_p < dat->n_p ; i_p++ )
//    {
//      tot_dGamma = 0;
//      tot_dGammag = 0;
//      for ( i_k = 0 ; i_k < dat->n_k ; i_k++ )
//	{
//	  p = dat->p_min + dat->dp * i_p;
//	  k = dat->k_min + 2 * dat->dp * i_k;
//	  p *= e_scale;
//	  k *= e_scale;
//	  here = i_p + i_k * dat->n_p;
//	  dGamma[here] = jacobian * 
//	    use_table ( p , k , gam->dGamma , 0 );
//	  dGamma_em[here] = jacobian *
//	    use_table ( p , k , gam->dGamma_em , 3 );
//	  if ( 4 * i_k - 2 * dat->n_kmin <= 3 + i_p + dat->n_pmin )
//	    { /* That is, if k < p/2 */
//	      dGamma_gqq[here] = dat->Nf * jacobian *
//		use_table ( p , k , gam->dGamma_gqq , 1 );
//	      dGamma_ggg[here] = jacobian *
//		use_table ( p , k , gam->dGamma_ggg , 2 );
//	    }
//	  else
//	    {
//	      dGamma_gqq[here] = 0;
//	      dGamma_ggg[here] = 0;
//	    }
//	  compton = T * jacobian * dat->cf / ( p * 2 * dat->dp * 32 * PI ) * 
//	    ( 0.5 * log(2*p/(dat->cf * dat->alpha_s * PI)) - 0.3614902 );
//	  /* The rate for quark->gluon by compton is,
//	     rate = [T^2 C_f^2 / 32 pi p] [ log[sqrt(pT)/gT] - C ]
//	     
//	     To turn that into a differential rate in k, divide
//	     by the k spacing, 1/2 dat->dp. 
//	     Use only when k=0 (gluons) or k=p (quarks).
//	     Only one power of T appears because p,k are scaled by T
//	     at this point.  We wrote one C_f to ease the use 
//	     in all three processes. */
//	  if ( dat->do_comptong )
//	    if ( ABS(k) < 1.5 * dat->dp * e_scale )
//	      { /* Then include Compton for gluons. */
//		/* We do this for k=0 rather than k=p because the 
//		   integrand only runs over k<p/2 and takes 2 final
//		   states with k,p-k energies.  So to get it right,
//		   we have to use k=0 and rely on the p-k to give
//		   the final state quark. */
//		dGamma_gqq[here] += dat->Nf * dat->df 
//		  * dat->cf * compton / ( (double) dat->da);
//		/* Should include an extra *2 but instead we do it
//		   for 2 rather than 1 value of k, namely just >and< 0. */
//	      }
//	  if ( ABS(k-p) < 1.5 * dat->dp * e_scale )
//	    { /* Then include Compton for quarks */
//	      compton *= ( (ABS(k-p)>0.4 * dat->dp * e_scale) ? 0.5 : 1 );
//	      if ( dat->do_compton )
//		dGamma[here] += dat->cf * compton; /* Rate quarks -> gluons */
//	    }
//	  if ( dat->photon_22 )
//	    dGamma_em[here] += compton;
//	  if ( dat->collisional_only )
//	    { /* Then undo all that work. */
//	      tot_dGamma -= dGamma[here];
//	      dGamma[here] = 0;
//	      dGamma_gqq[here] = 0;
//	      dGamma_ggg[here] = 0;
//	    }
//#ifdef NOISY
//	  if ( i_p == 100 && (p + 4 > k ) )
//	    {
//	      /* Then tell the world what is happening. */
//	      printf ( "p %e k %e qqg %e gqq %e ggg %e qqphot %e\n" , 
//		       p , k , dGamma[here] , dGamma_gqq[here] 
//		       , dGamma_ggg[here] , dGamma_em[here] );
//	    }
//#endif
//	  if ( dat->do_collisional )
//	    { /* Include "diffusion" style collisional energy loss */
//	      tmp1 = PI * (dat->alpha_s * dat->alpha_s) * T * T;
//	      tmp2 = log( p/mg2 ) - Gamma_E - .569961 - 1.662;
//	      dEdx_q = dat->Nf * 2/9.0 * tmp1 * ( tmp2 + log(2.) + 23/12.0 );
//	      dEdx_q += 4.0/3.0 * tmp1 * ( tmp2 + 13.0/6.0 );
//	      dEdx_g = dat->Nf * 0.5 * tmp1 * ( tmp2 + log(2.) + 13.0/6.0 );
//	      dEdx_g += 3 * tmp1 * ( tmp2 + 131.0/48.0 );
//	      if ( dEdx_q < 0 ) dEdx_q = 0; /* Can happen at small p */
//	      if ( dEdx_g < 0 ) dEdx_g = 0; /* due to log going negative */
//	      dEdx_q *= jacobian / ( gs2 * gs2 );
//	      dEdx_g *= jacobian / ( gs2 * gs2 );
//	      
//	      /*** diffrential rate ***/
//	      dEdx_q = dEdx_q / (2.0 * delta_kp * T);
//	      dEdx_g = dEdx_g / (2.0 * delta_kp * T);
//	      
//	      if (k < 1.5 * delta_kp && k > 0.5 * delta_kp)
//		{
//		  dGamma[here] += (1.0 + nB) / (delta_kp * T) * dEdx_q;
//		  dGamma_ggg[here] += (1.0 + nB) / (delta_kp * T) * dEdx_g;
//		}
//	      
//	      if (k > -1.5 * delta_kp && k < -0.5 * delta_kp)
//		{
//		  dGamma[here] += nB / (delta_kp * T) * dEdx_q;
//		  dGamma_ggg[here] += nB / (delta_kp * T) * dEdx_g;
//		}
//	      
//	    }
//	  tot_dGamma += dGamma[here];
//	  tot_dGammag += dGamma_gqq[here] + dGamma_ggg[here];
//	}
//      if ( tot_dGamma > dat->dx_max )
//	dat->dx_max = tot_dGamma;
//      if ( dat->include_gluons && tot_dGammag > dat->dx_max )
//	dat->dx_max = tot_dGammag;
//    }
//  dat->dx_max *= 2 * dat->dp;
//  dat->dx_max = 0.8 / dat->dx_max;
//  /* Value considered sufficient to avoid sign alternation problems. */
//}
//
////int main ()
////{
////  double     *dGamma , *dGamma_gqq , *dGamma_ggg , *dGamma_em;
////  double     x1 , x2 , x3;
////  Gamma_info dat;
////  dGammas    *Gam;
////  double     *P , *Pg , *Pem , *P_spare , x , T;
////  int        n_steps , i;
////  char       c , c_rescale , output_fname[100];
////  FILE       *output_file;
////
////  Gam = (dGammas *) malloc ( sizeof ( dGammas ) );
////  prep_equipment ( &dat , Gam , &dGamma , &dGamma_gqq , &dGamma_ggg , 
////		   &dGamma_em , &P , &Pg , &Pem);
//////  printf ( "Enter the temperature of the material in GeV.\n" );
//////  scanf ( "%lf" , &T ); READ_TO_EOLN ( stdin );
////}
