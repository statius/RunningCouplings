#include "smdr_internal.h"

/* This should appear in one source file: */
#include "smdr_pdg.h"

#define ZEROSAFE(a) (((a) > (SMDR_TOL)) ? (a) : (SMDR_TOL))

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Loads the MSbar input parameters from "static" variables
   X_in to the corresponding variables X. */

void SMDR_Load_Inputs ()
{
  Q = Q_in;
  g3 = g3_in; 
  g = g_in; 
  gp = gp_in;
  yt = yt_in; 
  yb = yb_in; 
  yc = yc_in; 
  ys = ys_in; 
  yu = yu_in; 
  yd = yd_in;
  ytau = ytau_in; 
  ymu = ymu_in; 
  ye = ye_in;
  k = k_in; 
  m2 = m2_in; 
  v = v_in; 
  Lambda = Lambda_in;
  Delta_alpha_had_5_MZ = Delta_alpha_had_5_MZ_in;

  return;
};

/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* Stores the current MSBar input parameters X to "static" inputs X_in. */
void SMDR_Save_Inputs ()
{
  Q_in = Q;
  g3_in = g3; 
  g_in = g; 
  gp_in = gp;
  yt_in = yt; 
  yb_in = yb; 
  yc_in = yc; 
  ys_in = ys; 
  yu_in = yu; 
  yd_in = yd;
  ytau_in = ytau; 
  ymu_in = ymu; 
  ye_in = ye;
  k_in = k; 
  m2_in = m2; 
  v_in = v; 
  Lambda_in = Lambda;
  Delta_alpha_had_5_MZ_in = Delta_alpha_had_5_MZ;

  return;
};

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Updates the minimal set of parameter combinations needed in betas.c */

void SMDR_update_for_betas (void)
{
  g32 = g3*g3;
  g34 = g32*g32;
  g36 = g34*g32;
  g38 = g36*g32;

  yt2 = yt*yt;
  yt4 = yt2*yt2;
  yt6 = yt4*yt2;
  yt8 = yt6*yt2;

  yb2 = yb*yb;
  yb4 = yb2*yb2;
  yb6 = yb4*yb2;
  yb8 = yb6*yb2;

  ytau2 = ytau*ytau; 
  ytau4 = ytau2*ytau2; 
  ytau6 = ytau4*ytau2; 
  ytau8 = ytau6*ytau2; 

  g2 = g*g;
  g4 = g2*g2;
  g6 = g4*g2;
  g8 = g6*g2;

  gp2 = gp*gp;
  gp4 = gp2*gp2;
  gp6 = gp4*gp2;
  gp8 = gp6*gp2;

  k2 = k*k;
  k3 = k2*k;
  k4 = k3*k;

  return;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Updates various other useful parameter combos. */

void SMDR_Update (void)
{
  SMDR_update_for_betas();

  v2 = v*v;
  v4 = v2*v2;
  Q2 = Q*Q;

  g10 = g8*g2;
  g12 = g10*g2;
  g14 = g12*g2;
  g16 = g14*g2;
  g18 = g16*g2;
  g20 = g18*g2;

  gp10 = gp8*gp2;
  gp12 = gp10*gp2;
  gp14 = gp12*gp2;
  gp16 = gp14*gp2;
  gp18 = gp16*gp2;

  yt10 = yt8*yt2;
  yt12 = yt10*yt2;

  k5 = k4*k;
  k6 = k5*k;
  k7 = k6*k;

  g2pgp2 = g2 + gp2;
  sqrtg2pgp2 = TVIL_SQRT(g2pgp2);
  g2pgp22 = g2pgp2 * g2pgp2;
  g2pgp23 = g2pgp22 * g2pgp2;
  g2pgp24 = g2pgp23 * g2pgp2;
  g2pgp25 = g2pgp24 * g2pgp2;

  g2m2k = g2 - 2*k;
  g2m2k2 = g2m2k * g2m2k;

  g2m4k = g2 - 4*k;
  g2m4k2 = g2m4k * g2m4k;

  g2m8k = g2 - 8*k;
  g2m8k2 = g2m8k * g2m8k;
  g2m8k3 = g2m8k2 * g2m8k;

  g2p8k = g2 + 8*k;
  g2p8k2 = g2p8k * g2p8k;

  g2mgp2 = g2 - gp2;
  g2mgp22 = g2mgp2 * g2mgp2;

  g2m2yt2 = g2 - 2*yt2;
  g2m2yt22 = g2m2yt2 * g2m2yt2;

  g2pgp2m2k = g2 + gp2 - 2*k;
  g2pgp2m2k2 = g2pgp2m2k * g2pgp2m2k;

  g2pgp2m4yt2 = g2 + gp2 - 4*yt2;
  g2pgp2m4yt22 = g2pgp2m4yt2 * g2pgp2m4yt2;

  g2pgp2m8k = g2 + gp2 - 8*k;
  g2pgp2m8k2 = g2pgp2m8k * g2pgp2m8k;
  g2pgp2m8k3 = g2pgp2m8k2 * g2pgp2m8k;

  g4pg2gp2mgp4 = g4 + g2*gp2 - gp4;
  g4pg2gp2mgp42 = g4pg2gp2mgp4 * g4pg2gp2mgp4;

  twog2pgp2 = 2*g2 + gp2;
  twog2pgp22 = twog2pgp2 * twog2pgp2;
  twog2pgp23 = twog2pgp22 * twog2pgp2;

  threeg2m5gp2 = 3*g2 - 5*gp2;
  threeg2m5gp22 = threeg2m5gp2 * threeg2m5gp2;

  threeg2mgp2 = 3*g2 - gp2;
  threeg2mgp22 = threeg2mgp2 * threeg2mgp2;

  threeg2pgp2 = 3*g2 + gp2;
  threeg2pgp22 = threeg2pgp2 * threeg2pgp2;

  auR = -(2./3.) * gp2/sqrtg2pgp2;
  adR =  (1./3.) * gp2/sqrtg2pgp2;
  aeR = gp2/sqrtg2pgp2;
  auL = 0.5 * sqrtg2pgp2 + auR;
  adL = -0.5 * sqrtg2pgp2 + adR;
  anL = 0.5 * sqrtg2pgp2;
  aeL = -0.5 * sqrtg2pgp2 + aeR;

  auL2 = auL*auL;
  auR2 = auR*auR;
  adL2 = adL*adL;
  adR2 = adR*adR;
  aeL2 = aeL*aeL;
  aeR2 = aeR*aeR;
  anL2 = anL*anL;

  W = (g2 * v2)/4;
  Z = (g2pgp2 * v2)/4;
  G = m2 + k * v2;
  H = m2 + 3 * k * v2;
  h = 2 * k * v2;
  T = yt2 * v2/2;
  b = yb2 * v2/2;
  tau = ytau2 * v2/2;

  h2 = h*h; 
  h3 = h2*h; 
  h4 = h3*h;
  h5 = h4*h;

  W2 = W*W; 
  W3 = W2*W; 
  W4 = W3*W; 
  W5 = W4*W; 
  W6 = W5*W;
  W7 = W6*W;
  
  Z2 = Z*Z; 
  Z3 = Z2*Z; 
  Z4 = Z3*Z; 
  Z5 = Z4*Z;
  Z6 = Z5*Z;
  
  T2 = T*T; 
  T3 = T2*T; 
  T4 = T3*T;

  b2 = b*b;
  tau2 = tau*tau;

  return;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Prints program name and version number. */

/* void SMDR_Version (void)  */
/* { */
/*   printf("(* %s Version %s *)\n", SMDR_NAME, SMDR_VERSION); */
/* } */

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_Ranges (void)
{
  SMDR_Check_g3_Range (0.5, 1.5);
  SMDR_Check_g_Range (0.4, 0.8);
  SMDR_Check_gp_Range (0.25, 0.5);
  SMDR_Check_lambda_Range (0.05, 0.2);
  SMDR_Check_yb_Range (0.05);
  SMDR_Check_ytau_Range (0.02);
  SMDR_Check_VEV_Range (199, 301);
  SMDR_Check_Q_Range (49, 501);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_Q_Range (SMDR_REAL Q_range_lo, SMDR_REAL Q_range_hi)
{
  if (Q < Q_range_lo)
    fprintf (stderr, 
    "WARNING: the input renormalization scale Q = %Lf seems too low.\n",Q);

  if (Q > Q_range_hi)
    fprintf (stderr, 
    "WARNING: the input renormalization scale Q = %Lf seems too high.\n",Q);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_VEV_Range (SMDR_REAL vev_range_lo, SMDR_REAL vev_range_hi)
{
  if (v < vev_range_lo)
    fprintf (stderr, 
    "WARNING: the input VEV, v = %Lf, seems too low.\n",v);

  if (v > vev_range_hi)
    fprintf (stderr, 
    "WARNING: the input VEV, v = %Lf, seems too high.\n",v);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_m2_Range (SMDR_REAL m2_range_lo, SMDR_REAL m2_range_hi)
{
  if (m2 < m2_range_lo)
    fprintf (stderr, 
    "WARNING: the input Higgs squared mass parameter, m2 = %Lf, seems too low.\n",m2);

  if (m2 > m2_range_hi)
    fprintf (stderr, 
    "WARNING: the input Higgs squared mass parameter, m2 = %Lf, seems too high.\n",m2);

  if (m2 > 0) fprintf (stderr, "Note that m2 should be negative.\n");
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_lambda_Range (SMDR_REAL k_range_lo, SMDR_REAL k_range_hi)
{
  if (k < k_range_lo)
    fprintf (stderr, 
    "WARNING: the input Higgs self coupling, lambda = %Lf, seems too low.\n",k);

  if (k > k_range_hi)
    fprintf (stderr, 
    "WARNING: the input Higgs self coupling, lambda = %Lf, seems too high.\n",k);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_yt_Range (SMDR_REAL yt_range_lo, SMDR_REAL yt_range_hi)
{
  if (yt < yt_range_lo)
    fprintf (stderr, 
    "WARNING: the input top-quark Yukawa coupling, yt = %Lf, seems too low.\n",yt);

  if (yt > yt_range_hi)
    fprintf (stderr, 
    "WARNING: the input top-quark Yukawa coupling, yt = %Lf, seems too high.\n",yt);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_yb_Range (SMDR_REAL yb_range_hi)
{
  if (yb < -TSIL_TOL)
    fprintf (stderr, 
    "WARNING: the input bottom-quark Yukawa coupling, yb = %Lf, is negative.\n",yb);

  if (yb > yb_range_hi)
    fprintf (stderr, 
    "WARNING: the input bottom-quark Yukawa coupling, yb = %Lf, seems too high.\n",yb);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_ytau_Range (SMDR_REAL ytau_range_hi)
{
  if (ytau < -TSIL_TOL)
    fprintf (stderr, 
    "WARNING: the input tau Yukawa coupling, ytau = %Lf, is negative.\n",ytau);

  if (ytau > ytau_range_hi)
    fprintf (stderr, 
    "WARNING: the input tau Yukawa coupling, ytau = %Lf, seems too high.\n",ytau);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_g3_Range (SMDR_REAL g3_range_lo, SMDR_REAL g3_range_hi)
{
  if (g3 < g3_range_lo)
    fprintf (stderr, 
    "WARNING: the input QCD gauge coupling, g3 = %Lf, seems too low.\n",g3);

  if (g3 > g3_range_hi)
    fprintf (stderr, 
    "WARNING: the input QCD gauge coupling, g3 = %Lf, seems too high.\n",g3);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_g_Range (SMDR_REAL g_range_lo, SMDR_REAL g_range_hi)
{
  if (g < g_range_lo)
    fprintf (stderr, 
    "WARNING: the input SU(2)_L gauge coupling, g = %Lf, seems too low.\n",g);

  if (g > g_range_hi)
    fprintf (stderr, 
    "WARNING: the input SU(2)_L gauge coupling, g = %Lf, seems too high.\n",g);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_gp_Range (SMDR_REAL gp_range_lo, SMDR_REAL gp_range_hi)
{
  if (gp < gp_range_lo)
    fprintf (stderr, 
    "WARNING: the input U(1)_Y gauge coupling, gp = %Lf, seems too low.\n",gp);

  if (gp > gp_range_hi)
    fprintf (stderr, 
    "WARNING: the input U(1)_Y gauge coupling, gp = %Lf, seems too high.\n",gp);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_Mhpole_Range (SMDR_REAL Mhpole, 
                             SMDR_REAL Mh_range_lo, SMDR_REAL Mh_range_hi)
{
  if (Mhpole < Mh_range_lo)
    fprintf (stderr, 
    "WARNING: the input Higgs pole mass, Mhpole = %Lf, seems too low.\n",Mhpole);

  if (Mhpole > Mh_range_hi)
    fprintf (stderr, 
    "WARNING: the input Higgs pole mass, Mhpole = %Lf, seems too high.\n",Mhpole);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Check_Mtpole_Range (SMDR_REAL Mtpole, 
                             SMDR_REAL Mt_range_lo, SMDR_REAL Mt_range_hi)
{
  if (Mtpole < Mt_range_lo)
    fprintf (stderr, 
    "WARNING: the input top-quark pole mass, Mtpole = %Lf, seems too low.\n",Mtpole);

  if (Mtpole > Mt_range_hi)
    fprintf (stderr, 
    "WARNING: the input top-quark pole mass, Mtpole = %Lf, seems too high.\n",Mtpole);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Returns TRUE (1) if the numbers are equal and FALSE (0) if they are
   not equal.  */

int SMDR_FPCompare (SMDR_REAL x, SMDR_REAL y)
{
  SMDR_REAL tmp;
  SMDR_REAL absx, absy;

  absx = TSIL_FABS(x); absy = TSIL_FABS(y);

  /* First check for 0 = 0? */
  if (absx < 1000*TSIL_TOL) {
    if (absy < 1000*TSIL_TOL) return TRUE;
    else return FALSE;
  }

  /* Make x the one with the larger abs value: */
  if (absx < absy) {
    tmp = y; y = x; x = tmp;
  }
  if (TSIL_FABS(x-y) < absx*TSIL_TOL) return TRUE;
  else return FALSE;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Gives the signed square root of a possibly negative real number    */

SMDR_REAL SMDR_SGNSQRT (SMDR_REAL x)
{
  if (x < 0.0L)
    return (-SMDR_SQRT(-x));
  else
    return (SMDR_SQRT(x));
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

void SMDR_Start_Timer (void)
{
  SMDR_Time_Start = clock();
  SMDR_Time_Total = 0;
  return;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

double SMDR_Timer (void)
{
  double timediff; 

  SMDR_Time_End = clock();
  timediff = difftime(SMDR_Time_End, SMDR_Time_Start)/CLOCKS_PER_SEC;
  /*
  printf("  Time used: %.2f seconds\n", timediff);
  */
  SMDR_Time_Start = clock();
  SMDR_Time_Total += timediff;
 
  return (timediff);
}

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

void SMDR_Eval_QCDQED_at_MZ (SMDR_REAL Q_MZ, 
                             SMDR_REAL Q_dec_thZW,
                             int loopOrder)
{
  SMDR_RGeval_QCDQED_53 (Q_MZ, Q_dec_thZW, loopOrder);

  SMDR_alphaS_5_MZ = alphaS_53;
  SMDR_alpha_5_MZ = alpha_53;
  SMDR_mb_MZ = mb_53;
  SMDR_mc_MZ = mc_53;
  SMDR_ms_MZ = ms_53;
  SMDR_mu_MZ = mu_53;
  SMDR_md_MZ = md_53;
  SMDR_mtau_MZ = mtau_53;
  SMDR_mmuon_MZ = mmuon_53;
  SMDR_melectron_MZ = melectron_53;

  return;
}
