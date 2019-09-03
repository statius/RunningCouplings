#ifndef _SMDR_H_
#define _SMDR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <time.h>

/* For convenience: */
typedef long double          SMDR_REAL;
typedef long double _Complex SMDR_COMPLEX;
#define SMDR_EXP             expl
#define SMDR_CEXP            cexpl
#define SMDR_LOG             logl
#define SMDR_CLOG            clogl
#define SMDR_FABS            fabsl
#define SMDR_CABS            cabsl
#define SMDR_SQRT            sqrtl
#define SMDR_CSQRT           csqrtl
#define SMDR_POW             powl
#define SMDR_CPOW            cpowl
#define SMDR_ATAN            atanl
#define SMDR_ATAN2           atan2l
#define SMDR_CREAL           creall
#define SMDR_CIMAG           cimagl
#define SMDR_CONJ            conjl
#define SMDR_EPSILON         LDBL_EPSILON
#define SMDR_TOL             1000.0L*LDBL_EPSILON

/* Standard Model MSbar input parameters, static: */
SMDR_REAL SMDR_Q_in;      /* MSbar renormalization scale */
SMDR_REAL SMDR_g3_in;     /* SU(3) gauge coupling */
SMDR_REAL SMDR_gp_in;     /* U(1)_Y gauge coupling */
SMDR_REAL SMDR_g_in;      /* SU(2) gauge coupling */
SMDR_REAL SMDR_yt_in;     /* Top Yukawa coupling */
SMDR_REAL SMDR_yb_in;     /* Bottom Yukawa coupling */
SMDR_REAL SMDR_yc_in;     /* Charm Yukawa coupling */
SMDR_REAL SMDR_ys_in;     /* Strange Yukawa coupling */
SMDR_REAL SMDR_yu_in;     /* Up Yukawa coupling */
SMDR_REAL SMDR_yd_in;     /* Down Yukawa coupling */
SMDR_REAL SMDR_ytau_in;   /* Tau Yukawa coupling */
SMDR_REAL SMDR_ymu_in;    /* Muon Yukawa coupling */
SMDR_REAL SMDR_ye_in;     /* Electron Yukawa coupling */
SMDR_REAL SMDR_lambda_in; /* Higgs self coupling */
SMDR_REAL SMDR_m2_in;     /* Higgs lagrangian mass^2 parameter */
SMDR_REAL SMDR_v_in;      /* Higgs VEV */
SMDR_REAL SMDR_Lambda_in; /* Field-independent vacuum energy */
SMDR_REAL SMDR_Delta_alpha_had_5_MZ_in; /* Non-perturbative hadronic 
                                           contribution to alpha */

/* Standard Model MSbar parameters, variable: */
SMDR_REAL SMDR_Q;      /* MSbar renormalization scale */
SMDR_REAL SMDR_g3;     /* SU(3) gauge coupling */
SMDR_REAL SMDR_gp;     /* U(1)_Y gauge coupling */
SMDR_REAL SMDR_g;      /* SU(2) gauge coupling */
SMDR_REAL SMDR_yt;     /* Top Yukawa coupling */
SMDR_REAL SMDR_yb;     /* Bottom Yukawa coupling */
SMDR_REAL SMDR_yc;     /* Charm Yukawa coupling */
SMDR_REAL SMDR_ys;     /* Strange Yukawa coupling */
SMDR_REAL SMDR_yu;     /* Up Yukawa coupling */
SMDR_REAL SMDR_yd;     /* Down Yukawa coupling */
SMDR_REAL SMDR_ytau;   /* Tau Yukawa coupling */
SMDR_REAL SMDR_ymu;    /* Muon Yukawa coupling */
SMDR_REAL SMDR_ye;     /* Electron Yukawa coupling */
SMDR_REAL SMDR_lambda; /* Higgs self coupling */
SMDR_REAL SMDR_m2;     /* Higgs lagrangian mass^2 parameter */
SMDR_REAL SMDR_v;      /* Higgs VEV */
SMDR_REAL SMDR_Lambda; /* Field-independent vacuum energy */
SMDR_REAL SMDR_Delta_alpha_had_5_MZ; /* Non-perturbative hadronic 
                                        contribution to alpha */

/* 5-quark, 3-lepton QCD+QED theory MSbar parameters */
SMDR_REAL SMDR_Q_53;       /* MSbar renormalization scale */
SMDR_REAL SMDR_alphaS_53;  /* QCD coupling */
SMDR_REAL SMDR_alpha_53;   /* QED coupling */
SMDR_REAL SMDR_mb_53;      /* bottom mass */
SMDR_REAL SMDR_mc_53;      /* charm mass */
SMDR_REAL SMDR_ms_53;      /* strange mass */
SMDR_REAL SMDR_mu_53;      /* up mass */
SMDR_REAL SMDR_md_53;      /* down mass */
SMDR_REAL SMDR_mtau_53;      /* tau mass */
SMDR_REAL SMDR_mmuon_53;      /* muon mass */
SMDR_REAL SMDR_melectron_53;      /* electron mass */

/* 4-quark, 3-lepton QCD+QED theory MSbar parameters */
SMDR_REAL SMDR_Q_43;       /* MSbar renormalization scale */
SMDR_REAL SMDR_alphaS_43;  /* QCD coupling */
SMDR_REAL SMDR_alpha_43;   /* QED coupling */
SMDR_REAL SMDR_mc_43;      /* charm mass */
SMDR_REAL SMDR_ms_43;      /* strange mass */
SMDR_REAL SMDR_mu_43;      /* up mass */
SMDR_REAL SMDR_md_43;      /* down mass */
SMDR_REAL SMDR_mtau_43;      /* tau mass */
SMDR_REAL SMDR_mmuon_43;      /* muon mass */
SMDR_REAL SMDR_melectron_43;      /* electron mass */

/* 4-quark, 2-lepton QCD+QED theory MSbar parameters */
SMDR_REAL SMDR_Q_42;       /* MSbar renormalization scale */
SMDR_REAL SMDR_alphaS_42;  /* QCD coupling */
SMDR_REAL SMDR_alpha_42;   /* QED coupling */
SMDR_REAL SMDR_mc_42;      /* charm mass */
SMDR_REAL SMDR_ms_42;      /* strange mass */
SMDR_REAL SMDR_mu_42;      /* up mass */
SMDR_REAL SMDR_md_42;      /* down mass */
SMDR_REAL SMDR_mmuon_42;      /* muon mass */
SMDR_REAL SMDR_melectron_42;      /* electron mass */

/* 3-quark, 2-lepton QCD+QED theory MSbar parameters */
SMDR_REAL SMDR_Q_32;       /* MSbar renormalization scale */
SMDR_REAL SMDR_alphaS_32;  /* QCD coupling */
SMDR_REAL SMDR_alpha_32;   /* QED coupling */
SMDR_REAL SMDR_ms_32;      /* strange mass */
SMDR_REAL SMDR_mu_32;      /* up mass */
SMDR_REAL SMDR_md_32;      /* down mass */
SMDR_REAL SMDR_mmuon_32;      /* muon mass */
SMDR_REAL SMDR_melectron_32;      /* electron mass */

/* Standard Model "on-shell" parameters and other outputs */
SMDR_REAL SMDR_GFermi;
SMDR_REAL SMDR_GFermi_DGG; /* Change or remove this eventually? */
SMDR_REAL SMDR_alpha;      /* Sommerfeld fine structure constant */
SMDR_REAL SMDR_alphaS_MZ;  /* In non-decoupled Standard Model at Q=MZ */
SMDR_REAL SMDR_alpha_MZ;   /* In non-decoupled Standard	Model at Q=MZ */
SMDR_REAL SMDR_s2W_MZ;     /* In non-decoupled Standard	Model at Q=MZ */
SMDR_REAL SMDR_alpha_MZ_PDG; /* t decoupled but not W */
SMDR_REAL SMDR_s2W_MZ_PDG;   /* t decoupled but not W */
SMDR_REAL SMDR_alphaS_5_MZ; /* t,h,Z,W all decoupled */
SMDR_REAL SMDR_alpha_5_MZ;  /* t,h,Z,W all decoupled */
SMDR_REAL SMDR_mb_MZ;
SMDR_REAL SMDR_mc_MZ;
SMDR_REAL SMDR_ms_MZ;
SMDR_REAL SMDR_mu_MZ;
SMDR_REAL SMDR_md_MZ;
SMDR_REAL SMDR_mtau_MZ;
SMDR_REAL SMDR_mmuon_MZ;
SMDR_REAL SMDR_melectron_MZ;
SMDR_REAL SMDR_Mt_pole;
SMDR_REAL SMDR_Gammat_pole;
SMDR_REAL SMDR_Mh_pole;
SMDR_REAL SMDR_Gammah_pole;
SMDR_REAL SMDR_MZ_pole;
SMDR_REAL SMDR_GammaZ_pole;
SMDR_REAL SMDR_MZ_BreitWigner;
SMDR_REAL SMDR_GammaZ_BreitWigner;
SMDR_REAL SMDR_MW_pole;
SMDR_REAL SMDR_GammaW_pole;
SMDR_REAL SMDR_MW_BreitWigner;
SMDR_REAL SMDR_GammaW_BreitWigner;
SMDR_REAL SMDR_Mb_pole;
SMDR_REAL SMDR_mbmb;
SMDR_REAL SMDR_mcmc;
SMDR_REAL SMDR_ms_2GeV;
SMDR_REAL SMDR_mu_2GeV;
SMDR_REAL SMDR_md_2GeV;
SMDR_REAL SMDR_Mtau_pole;
SMDR_REAL SMDR_Mmuon_pole;
SMDR_REAL SMDR_Melectron_pole;

/* Experimental reference values from Review of Particle Properties */
/* (These are declared and have their values set elsewhere; see
   smdr_pdg.h.) */
extern SMDR_REAL SMDR_GFermi_EXPT;
extern SMDR_REAL SMDR_GFermi_EXPT_UNC;
extern SMDR_REAL SMDR_alpha_EXPT;          /* Fine structure constant*/
extern SMDR_REAL SMDR_alpha_EXPT_UNC;
extern SMDR_REAL SMDR_s2W_MZ_EXPT;         /* MSbar at Q=MZ, top decoupled */
extern SMDR_REAL SMDR_s2W_MZ_EXPT_UNC;
extern SMDR_REAL SMDR_alphaS_MZ_EXPT;      /* MSbar at Q=MZ, top decoupled */
extern SMDR_REAL SMDR_alphaS_MZ_EXPT_UNC;
extern SMDR_REAL SMDR_alpha_MZ_EXPT;       /* MSbar at Q=MZ, top decoupled */
extern SMDR_REAL SMDR_alpha_MZ_EXPT_UNC;
extern SMDR_REAL SMDR_Delta_alpha_had_5_MZ_EXPT;
extern SMDR_REAL SMDR_Delta_alpha_had_5_MZ_EXPT_UNC;
extern SMDR_REAL SMDR_Mt_EXPT;             /* Pole mass, has renormalon ambiguity */
extern SMDR_REAL SMDR_Mt_EXPT_UNC;
extern SMDR_REAL SMDR_Mh_EXPT;
extern SMDR_REAL SMDR_Mh_EXPT_UNC;
extern SMDR_REAL SMDR_MZ_EXPT;             /* Experimental Breit-Wigner mass */
extern SMDR_REAL SMDR_MZ_EXPT_UNC;
extern SMDR_REAL SMDR_MW_EXPT;             /* Experimental Breit-Wigner mass */
extern SMDR_REAL SMDR_MW_EXPT_UNC;
extern SMDR_REAL SMDR_mbmb_EXPT;           /* MSbar mass evaluated at itself. */
extern SMDR_REAL SMDR_mbmb_EXPT_UNC_hi;
extern SMDR_REAL SMDR_mbmb_EXPT_UNC_lo;
extern SMDR_REAL SMDR_mcmc_EXPT;           /* MSbar mass evaluated at itself. */
extern SMDR_REAL SMDR_mcmc_EXPT_UNC_hi;
extern SMDR_REAL SMDR_mcmc_EXPT_UNC_lo;
extern SMDR_REAL SMDR_ms_2GeV_EXPT;        /* MSbar mass at Q = 2 GeV */
extern SMDR_REAL SMDR_ms_2GeV_EXPT_UNC_hi;
extern SMDR_REAL SMDR_ms_2GeV_EXPT_UNC_lo;
extern SMDR_REAL SMDR_mu_2GeV_EXPT;        /* MSbar mass at Q = 2 GeV */
extern SMDR_REAL SMDR_mu_2GeV_EXPT_UNC_hi;
extern SMDR_REAL SMDR_mu_2GeV_EXPT_UNC_lo;
extern SMDR_REAL SMDR_md_2GeV_EXPT;        /* MSbar mass at Q = 2 GeV */
extern SMDR_REAL SMDR_md_2GeV_EXPT_UNC_hi;
extern SMDR_REAL SMDR_md_2GeV_EXPT_UNC_lo;
extern SMDR_REAL SMDR_Mtau_EXPT;           /* Pole mass */
extern SMDR_REAL SMDR_Mtau_EXPT_UNC;
extern SMDR_REAL SMDR_Mmuon_EXPT;          /* Pole mass */
extern SMDR_REAL SMDR_Mmuon_EXPT_UNC;
extern SMDR_REAL SMDR_Melectron_EXPT;      /* Pole mass */
extern SMDR_REAL SMDR_Melectron_EXPT_UNC;

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ------------------------- User API ----------------------------------- */
/* NOTE: all global variables referred to below are understood to have 

   SMDR_ 

   prepended to them, unless one has included the file
   "smdr_internal.h".  If you do choose to include "smdr_internal.h",
   beware of variable name conflicts!
*/

/* ---------------------------------------------------------------------- */
/* ------------------------- In RGrun.c --------------------------------- */
/* For the functions in this file, loopOrder should be 1, 2, 3, 4, 5.
   In the case of loopOrder = 5, the known partial 5-loop loop order
   results are included.
*/

/* Copies the input MSbar parameter global variables:

   Q_in, g3_in, g_in, ..., m2_in, v_in, Lambda_in, Delta_alpha_had_5_MZ_in 

   to the corresponding working global variables:

   Q, g3, g, ..., m2, v, Lambda, Delta_alpha_had_5_MZ.

   Then, runs all of these variables (except the last, which is
   defined at MZ) from the scale Q=Q_in to the scale Q_final. Uses
   beta functions with loop order loopOrder. The results are stored in
   the global variables:

   Q=Q_final, g3, g, gp, yt, yb, ..., m2, v, Lambda.

   This function works simply by calling SMDR_Load_Inputs() followed
   by SMDR_RGrun_SM().
*/
int SMDR_RGeval_SM (SMDR_REAL Q_final, int loopOrder);

/* Runs the full Standard Model MSbar parameters in the non-decoupled
   theory:

   Q, g3, g, gp, yt, yb, ..., m2, v, Lambda

   from the current scale Q to the scale Q_final. Uses beta functions
   with loop order loopOrder. The results are stored in the global
   variables: Q=Q_final, g3, gp, g, yt, yb, ..., m2, v, Lambda.
*/
int SMDR_RGrun_SM (SMDR_REAL Q_final, int loopOrder);

/* ---------------------------------------------------------------------- */
/* --------------------- In RGrun_QCDQED.c ------------------------------ */
/* For the functions in this file, loopOrder should be 1, 2, 3, 4, 5.
   In the case of loopOrder = 5, the known partial 5-loop loop order
   results are included.
*/

/* Evaluates, using RG running and decoupling, the 5-quark,
   3-charged-lepton QCD+QED effective theory MSbar running parameter
   global variables:

   alphaS, alpha, mb, mc, ms, md, mu, mtau, mmuon, melectron 

   at the scale Qfinal, taking the decoupling scale for top, Higgs, Z,
   W to be Q_thZW_dec. The inputs are the global variables
   corresponding to the non-decoupled Standard Model input MSbar
   parameters:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., lambda_in, v_in

   The results are stored as the global variables:

   Q_53=Qfinal, alphaS_53, alpha_53, m<fermion>_53

   with <fermion> = b, c, s, d, u, tau, muon, electron.  Typically,
   Q_final is expected to be roughly between 4 GeV and MZ.
*/
int SMDR_RGeval_QCDQED_53 (SMDR_REAL Qfinal,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder);

/* Evaluates, using RG running and decoupling, the 4-quark,
   3-charged-lepton QCD+QED effective theory MSbar running parameters:

   alphaS, alpha, mc, ms, md, mu, mtau, mmuon, melectron

   at the scale Qfinal, taking decoupling scales of Q_thZW_dec for
   top, Higgs, Z, W, and taking Q_b_dec for the bottom quark
   decoupling scale.  The inputs are the global variables
   corresponding to the non-decoupled Standard Model input MSbar
   parameters:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., lambda_in, v_in.

   The results are stored as the global variables:

   Q_43=Qfinal, alphaS_43, alpha_43, m<fermion>_43,

   with <fermion> = c, s, d, u, tau, muon, electron.  Typically,
   Q_final is expected to be roughly between 1.5 GeV and 5 GeV.
*/
int SMDR_RGeval_QCDQED_43 (SMDR_REAL Qfinal,
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder);

/* Evaluates, using RG running and decoupling, the 4-quark,
   2-charged-lepton QCD+QED effective theory MSbar running parameters:

   alphaS, alpha, mc, ms, md, mu, mmuon, melectron 

   at the scale Qfinal, taking decoupling scales of Q_thZW_dec for
   top, Higgs, Z, W, and taking Q_b_dec for the bottom quark
   decoupling scale, and taking Q_tau_dec for the tau lepton
   decoupling scale. The inputs are the global variables corresponding
   to the non-decoupled Standard Model input MSbar parameters:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., lambda_in, v_in.

   The results are stored as the global variables:

   Q_42=Qfinal, alphaS_42, alpha_42, m<fermion>_42

   with <fermion> = c, s, d, u, muon, electron.  Typically, Q_final is
   expected to be roughly between 1 GeV and 2 GeV.
*/
int SMDR_RGeval_QCDQED_42 (SMDR_REAL Qfinal,
                           SMDR_REAL Q_tau_dec,
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder);

/* Evaluates, using RG running and decoupling, the 3-quark,
   2-charged-lepton QCD+QED effective theory MSbar running parameters:

   alphaS, alpha, ms, mu, md, mmuon, melectron

   at the scale Qfinal, taking decoupling scales of Q_thZW_dec for
   top, Higgs, Z, W, and taking Q_b_dec for the bottom quark
   decoupling scale, and Q_tau_dec for the tau lepton decoupling
   scale, and Q_c_dec for the charm quark decoupling scale. The inputs
   are the non-decoupled Standard Model input MSbar parameters:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., lambda_in, v_in.

   The results are stored as the global variables:

   Q_32=Qfinal, alphaS_32, alpha_32, m<fermion>_32, 

   with <fermion> = s, d, u, muon, electron. Typically, Q_final is
   expected to be less than about 1.5 GeV.
*/
int SMDR_RGeval_QCDQED_32 (SMDR_REAL Qfinal,
                           SMDR_REAL Q_c_dec,
                           SMDR_REAL Q_tau_dec,
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder);

/* Runs the MSbar parameters of the 5-quark, 3-charged lepton QCD+QED
   effective theory:

   Q_53, alphaS_53, alpha_53, m<fermion>_53, 

   with fermion = b, c, s, d, u, tau, muon, electron, from the current
   scale Q_53 to the scale Q_final. Uses beta functions with loop
   order loopOrder. The results are stored in the same global
   variables: Q_53=Qfinal, alphaS_53, alpha_53, m<fermion>_53.
*/
int SMDR_RGrun_QCDQED_53 (SMDR_REAL Qfinal, int loopOrder);

/* Runs the MSbar parameters of the 4-quark, 3-charged lepton QCD+QED
   effective theory:

   Q_43, alphaS_43, alpha_43, m<fermion>_43, 

   with <fermion> = c, s, d, u, tau, muon, electron, from the current
   scale Q_43 to the scale Q_final. Uses beta functions with loop
   order loopOrder. The results are stored in the same global
   variables: Q_43=Qfinal, alphaS_43, alpha_43, m<fermion>_43.
*/
int SMDR_RGrun_QCDQED_43 (SMDR_REAL Qfinal, int loopOrder);

/* Runs the MSbar parameters of the 4-quark, 2-charged lepton QCD+QED
   effective theory:

   Q_42, alphaS_42, alpha_42, m<fermion>_42, 

   with <fermion> = c, s, d, u, muon, electron from the current scale
   Q_42 to the scale Q_final. Uses beta functions with loop order
   loopOrder. The results are stored in the same global variables:
   Q_42=Qfinal, alphaS_42, alpha_42, m<fermion>_42.
*/
int SMDR_RGrun_QCDQED_42 (SMDR_REAL Qfinal, int loopOrder);

/* Runs the MSbar parameters of the 3-quark, 2-charged lepton QCD+QED
   effective theory:

   Q_32, alphaS_32, alpha_32, m<fermion>_32,

   with <fermion> = s, d, u, muon, electron, from the current scale
   Q_32 to the scale Q_final. Uses beta functions with loop order
   loopOrder. The results are stored in the same global variables:
   Q_32=Qfinal, alphaS_32, alpha_32, m<fermion>_32.
*/
int SMDR_RGrun_QCDQED_32 (SMDR_REAL Qfinal, int loopOrder);

/* This function is a tool used by the preceding functions. 

   Runs the MSbar parameters of the QCD+QED theory with nu up-type
   quarks, nd down-type quarks, and ne charged leptons, from the scale
   Q_init to the scale Q_final. The values of the QCD and QED
   couplings at Q_init are alphaS_init and alpha_init. The output
   results at Q_final are: alphaS_final, alpha_final, and the ratios
   cu=mu_final/mu_init, cd=md_final/md_init, ce=me_final/me_init.
*/
int SMDR_RGrun_QCDQED (SMDR_REAL Q_init, SMDR_REAL Q_final,
                       int loopOrder, int nu, int nd, int ne,
                       SMDR_REAL alphaS_init, SMDR_REAL alpha_init,
                       SMDR_REAL *alphaS_final, SMDR_REAL *alpha_final,
                       SMDR_REAL *cu, SMDR_REAL *cd, SMDR_REAL *ce);

/* ---------------------------------------------------------------------- */
/* --------------------- In decouple_thZW.c --------------------------- */

/* Simultaneously decouples the top quark and the Higgs, Z, and W
   bosons, at the current MSbar renormalization scale Q. The inputs
   are the global variables corresponding to the full Standard Model:

   Q, g3, g, gp, yt, yb, yc, ys, yd, yu, ytau, ymuon, yelectron,
   lambda, v

   The results are stored in the global variables: 

   alphaS_53, alpha_53, and m<fermion>_53, 

   where <fermion> = b, c, s, u, d, tau, muon, electron.  These are
   the MSbar QCD and EM couplings and the running masses in the
   QCD+QED effective theory that contains only the 5 quarks b,c,s,u,d
   and the charged leptons tau,mu,e, and the gluons and the photons.
   At present, the argument loopOrder is ignored, and all known
   effects are included.
*/
void SMDR_Decouple_thZW (int loopOrder);

/* ---------------------------------------------------------------------
   Evaluates:
     1) alphaS(MZ) and alphaS(MZ) and s^2(thetaW) in the MSbar scheme
        of the non-decoupled Standard Model theory, putting the
        results in the global variables: SMDR_alphaS_MZ,
        SMDR_alpha_MZ, SMDR_s2W_MZ.
     2) The Sommerfeld fine structure constant, using the results of
        1411.7040 Degrassi, Gambino, and Giardino. The result is
        stored in the global variable SMDR_alpha.
     3) alpha(MZ) and s^2(thetaW) in PDG MSbar scheme with only top
        decoupled but W not decoupled. The results are stored in the
        global variables: SMDR_alpha_MZ_PDG, SMDR_s2W_MZ_PDG.

   Uses as inputs the global variable Standard Model input values

   Q_in, g3_in, g_in, gp_in, yt_in, ...

   as well as the arguments Mtpole, Mhpole, MWpole.
*/
void SMDR_Eval_Gauge (SMDR_REAL Mtpole, SMDR_REAL Mhpole, SMDR_REAL MWpole);

/* ---------------------------------------------------------------------- */
/* ----------------------- In QCDQED_match.c ---------------------------- */
/* For the functions defined in this file:
   loopOrder = 0, 1, 2, 3, 4
   EM contributions beyond 2 loops are not included.
*/

/* Decouples the bottom quark at the scale given by the global
   variable Q_53 = the MSbar renormalization scale of the 5-quark,
   3-charged-lepton QCD+QED effective field theory. The inputs are the
   global variables:

   Q_53, alphaS_53, alpha_53, m<Fermion>_53

   with <Fermion> = b, c, s, d, u, tau, muon, electron.  The output
   MSbar parameters of the 4-quark, 3-charged-lepton QCD+QED effective
   theory are put in the global variables:

   Q_43=Q_53, alphaS_43, alpha_43, m<fermion>_43,

   with <fermion> = c, s, d, u, tau, muon, electron,
*/
void SMDR_Decouple_bottom (int loopOrder);

/* Decouples the tau lepton at the scale given by the global variable
   Q_43 = the MSbar renormalization scale of the 4-quark,
   3-charged-lepton QCD+QED effective field theory. The input are the
   global variables:

   Q_43, alphaS_43, alpha_43, m<Fermion>_43

   with <Fermion> = c, s, d, u, tau, muon, electron.  The output MSbar
   parameters of the 4-quark, 2-charged-lepton QCD+QED effective
   theory are put in the global variables:

   Q_42=Q_43, alphaS_42, alpha_42, m<fermion>_42,

   with <fermion> = c, s, d, u, muon, electron,
*/
void SMDR_Decouple_tau (int loopOrder);

/* Decouples the charm quark at the scale given by the global variable
   Q_42 = the MSbar renormalization scale of the 4-quark,
   2-charged-lepton QCD+QED effective field theory. The inputs are the
   global variables:

   Q_42, alphaS_42, alpha_42, m<Fermion>_42

   with <Fermion> = c, s, d, u, muon, electron.  The output MSbar
   parameters of the 3-quark, 2-charged-lepton QCD+QED effective
   theory are put in the global variables:

   Q_32=Q_42, alphaS_32, alpha_32, m<fermion>_32,

   with <fermion> = s, d, u, muon, electron,
*/
void SMDR_Decouple_charm (int loopOrder);

/* This function is a tool used by the preceding functions.
     The input arguments are:
     Fermion_type = "u", "d", or "e" = type of heavy fermion being decoupled
     Q_match   = MSbar decoupling scale.
     m_Fermion = MSbar mass of heavy fermion being decoupled, in full theory,
                 at Q_match
     alphaS_hi  = strong coupling in full theory, at scale Q_match.
     alpha_hi   = EM coupling in full theory, at scale Q_match
     nqlight = number of light quarks, not including the fermion being
               decoupled. Note that the number of light leptons doesn't matter
               in the approximation being used. (It would matter if 3-loop
               EM contributions were included.)

   The outputs are:
     alphaS_lo = strong coupling in decoupled theory, at scale Q_match.
     alpha_lo  = EM coupling in decoupled theory, at scale Q_match.
     zu = ratio of (MSbar light up-type mass in decoupled theory)/
                   (MSbar light up-type mass in full theory).
     zd = ratio of (MSbar light down-type mass in decoupled theory)/
                   (MSbar light down-type mass in full theory).
     ze = ratio of (MSbar light charged lepton mass in decoupled theory)/
                   (MSbar light charged lepton mass in full theory).

     Small mass corrections are not included here, but are added in
     the functions SMDR_Decouple_b and SMDR_Decouple_c.
*/
void SMDR_QCDQEDmatch (const char *Fermion_type,
                       SMDR_REAL Q_match,
                       SMDR_REAL m_Fermion,
                       SMDR_REAL alphaS_hi,
                       SMDR_REAL alpha_hi,
                       int nqlight,
                       int loopOrder,
                       SMDR_REAL *alphaS_lo,
                       SMDR_REAL *alpha_lo,
                       SMDR_REAL *zu,
                       SMDR_REAL *zd,
                       SMDR_REAL *ze);

/* ---------------------------------------------------------------------- */
/* -------------------------- In effpot.c: ------------------------------ */
/* For the functions defined in this file:   
   loopOrder = 0 = tree-level
   loopOrder = 1 = 1-loop
   loopOrder = 2 = 2-loop
   loopOrder = 2.5 = 2-loop plus 3-loop contributions in the large g3,yt limit
   loopOrder = 3 = full 3-loop
   loopOrder = 3.5 = full 3-loop plus 4-loop at leading order in QCD
*/

/* Returns the minimum value of the effective potential (in Landau
   gauge, in the MSbar scheme, with Goldstone boson resummation). 

    If the argument Q_eval is positive, then the inputs are obtained
    by first RG running the MSbar parameter global variables:

    Q_in, v_in, lambda_in, g3_in, g_in, gp_in, yt_in, yb_in, ytau_in

    to the scale Q_eval, using SMDR_RGeval_SM().

    If the argument Q_eval is negative, then the inputs are instead
    taken to be the current values of the MSbar parameter global
    working variables:

    Q, v, lambda, g3, g, gp, yt, yb, ytau

   The function adjusts the global working variable m2 to assure that 
   the resummed effective potential is minimized, and returns the 
   corresponding minimum value of the effective potential.
*/
SMDR_REAL SMDR_Eval_Veffmin (SMDR_REAL Q_eval, float loopOrder);


/* Returns the value of m2 that minimizes the effective potential (in
   Landau gauge, in the MSbar scheme, with Goldstone boson resummation). 
   
   If the argument Q_eval is positive, then the inputs are obtained
   by first RG running the MSbar parameter global variables:

   Q_in, v_in, lambda_in, g3_in, g_in, gp_in, yt_in, yb_in, ytau_in

   to the scale Q_eval, using SMDR_RGeval_SM().

   If the argument Q_eval is negative, then the inputs are instead
   taken to be the current values of the MSbar parameter global
   working variables:
   
   Q, v, lambda, g3, g, gp, yt, yb, ytau,

   which are assumed to be fixed. This function does not automatically
   change the global variable m2.
*/
SMDR_REAL SMDR_Eval_m2 (SMDR_REAL Q_eval, float loopOrder);


/* Returns the value of the VEV v at the minimum of the effective
   potential (in Landau gauge, in the MSbar scheme, with Goldstone
   boson resummation, at loop order loopOrder) 

   If the argument Q_eval is positive, then the inputs are obtained
   by first RG running the MSbar parameter global variables:

     Q_in, m2_in, lambda_in, g3_in, g_in, gp_in, yt_in, yb_in, ytau_in,

   to the scale Q_eval, using SMDR_RGeval_SM().

   If the argument Q_eval is negative, then the inputs are instead
   taken to be the current values of the MSbar parameter global
   working variables:

     Q, m2, lambda, g3, g, gp, yt, yb, ytau,

   which are assumed to be fixed. This function does not automatically
   change the global variable v.
*/
SMDR_REAL SMDR_Eval_vev (SMDR_REAL Q_eval, float loopOrder);


/* Returns the value of the quantity 

   Delta = -G = (vtree^2 - v^2)*lambda

   defined in section V.3 of 1709.02397. 

   If the argument Q_eval is positive, then the inputs are obtained
   by first RG running the MSbar parameter global variables:

   Q_in, v_in, lambda_in, g3_in, g_in, gp_in, yt_in, yb_in, ytau_in,

   to the scale Q_eval, using SMDR_RGeval_SM().

   If the argument Q_eval is negative, then the inputs are instead
   taken to be the current values of the MSbar parameter global
   working variables:

   Q, v, lambda, g3, g, gp, yt, yb, ytau,

   which are assumed to be fixed.
*/
SMDR_REAL SMDR_Eval_vevDelta (SMDR_REAL Q_eval, float loopOrder);

/* ---------------------------------------------------------------------- */
/* -------------------------- In Mh.c ----------------------------------- */
/* For the functions in this file, loopOrder may take the following
   values:
      0    tree level
      1    1-loop
      1.5  1-loop plus 2-loop QCD corrections
      2    full 2-loop
      2.3  full 2-loop plus leading 3-loop QCD
      2.5  full 2-loop plus leading 3-loop QCD and nonQCD terms.
    The computations use equations (2.43), (2.46)-(2.48), and
    (3.2)-(3.4) of 1407.4336.
*/

/*  Computes the Higgs pole mass. The results M_h and Gamma_h for the
    complex pole squared mass

    M_h^2 - i Gamma_h M_h 

    are returned as Mhpoleresult and Gammahpoleresult.

    If the argument Q_eval is positive, then the inputs are obtained
    by first RG running the MSbar parameter global variables:

    Q_in, g3_in, g_in, gp_in, lambda_in, yt_in, ..., v_in

    to the scale Q_eval, using SMDR_RGeval_SM().

    If the argument Q_eval is negative, then the inputs are instead
    given directly by the current values of the MSbar parameter global
    variables:

    Q, g3, g, gp, lambda, yt, ..., v.
*/
void SMDR_Eval_Mh_pole (SMDR_REAL Q_eval,
                        float loopOrder, 
                        SMDR_REAL *Mhpoleresult,
                        SMDR_REAL *Gammahpoleresult);

/*  Returns the Higgs self-coupling lambda, given the real part of the
    Higgs pole mass and the other MSbar Standard Model parameters,
    determined as follows:

    If the argument Q_eval is positive, then the inputs are obtained
    by first RG running the MSbar parameter global variables:

    Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

    to the scale Q_eval, using SMDR_RGeval_SM().

    If the argument Q_eval is negative, then the inputs are instead
    given by the current values of the MSbar parameter global
    variables:

    Q, g3, g, gp, yt, ..., v.
*/
SMDR_REAL SMDR_Eval_lambda (SMDR_REAL Q_eval, 
                            SMDR_REAL mhpole,
                            float loopOrder);

/* Controls the error in the iterated computation of the Higgs mass. */
SMDR_REAL SMDR_MHPOLE_TOLERANCE;

/* ---------------------------------------------------------------------- */
/* ------------------------- In MW.c ------------------------------------ */

/*  Computes the complex pole and Breit-Wigner squared masses of the W
    boson, using the calculation from 1503.03782.  The argument
    loopOrder may take the following values:

      0    tree level
      1    1-loop
      1.5  1-loop plus 2-loop QCD corrections
      2    full 2-loop

    If the argument Q_eval is positive, then the inputs are obtained
    by first RG running the MSbar parameter global variables:

    Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

    to the scale Q_eval, using SMDR_RGeval_SM().

    If the argument Q_eval is negative, then the inputs are instead
    given by the current values of the MSbar parameter global
    variables:

    Q, g3, g, gp, yt, ..., v.

    The results for the complex pole squared mass 

    M_W^2 - i Gamma_W M_W

    are returned in MWpoleresult, GammaWpoleresult.

    The results for the complex Breit-Wigner squared mass are also
    returned as MWBreitWignerresult, GammaWBreitWignerresult.
*/
void SMDR_Eval_MW_pole (SMDR_REAL Q_eval,
                        float loopOrder,
                        SMDR_REAL *MWpoleresult,
                        SMDR_REAL *GammaWpoleresult,
                        SMDR_REAL *MWBreitWignerresult,
                        SMDR_REAL *GammaWBreitWignerresult);

/* ---------------------------------------------------------------------- */
/* ------------------------- In MZ.c ------------------------------------ */

/*  Computes the complex pole and Breit-Wigner squared masses of the Z
    boson, using the calculation from 1505.04833. The argument
    loopOrder may take the following values:

      0    tree level
      1    1-loop
      1.5  1-loop plus 2-loop QCD corrections
      2    full 2-loop

    If the argument Q_eval is positive, then the inputs are obtained
    by first RG running the MSbar parameter global variables:

    Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

    to the scale Q_eval, using SMDR_RGeval_SM().

    If the argument Q_eval is negative, then the inputs are instead
    given by the current values of the MSbar parameter global
    variables:

    Q, g3, g, gp, yt, ..., v.

    The results for the complex pole squared mass

    M_Z^2 - i Gamma_Z M_Z

    are returned as MZpoleresult, GammaZpoleresult.

    The results for the complex Breit-Wigner squared mass are also
    returned as MZBreitWignerresult, GammaZBreitWignerresult.
*/
void SMDR_Eval_MZ_pole (SMDR_REAL Q_eval,
                        float loopOrder,
                        SMDR_REAL *MZpoleresult,
                        SMDR_REAL *GammaZpoleresult,
                        SMDR_REAL *MZBreitWignerresult,
                        SMDR_REAL *GammaZBreitWignerresult);

/* ---------------------------------------------------------------------- */
/* ------------------------------ In Mt.c: ------------------------------ */
/* For the functions in this file:
   method = 0 (expand around tree-level) or
            1 (expand around pole. Recommended; more accurate but slower).
   QCDLoopOrder = 0, 1, 2, 3, or 4.
   otherLoopOrder = 0 (tree-level) or
                    1 (1-loop) or
                    1.5 (2-loop mixed QCD/EW) or
                    2 (2-loop full).
   The sources for the relevant results are:
     4-loop pure QCD from 1502.01030, updated and perfected in 1606.06754.
     2-loop non-pure-QCD from 1604.01134.
*/

/* Computes the complex pole mass of the top quark. The results for
   the complex pole squared mass M_t^2 - i Gamma_t M_t are returned as
   Mtpoleresult, Gammatpoleresult.

   If the argument Q_eval is positive, then the inputs are obtained by
   first RG running the MSbar parameter global variables:

   Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

   to the scale Q_eval, using SMDR_RGeval_SM().

   If the argument Q_eval is negative, then the inputs are instead
   given by the current values of the MSbar parameter global
   variables:

   Q, g3, g, gp, yt, ..., v.
*/
void SMDR_Eval_Mt_pole (SMDR_REAL Q_eval,
                        int method, 
                        int QCDloopOrder, 
                        float otherLoopOrder,
                        SMDR_REAL *Mtpoleresult, 
                        SMDR_REAL *Gammatpoleresult);

/* Returns the top-quark Yukawa coupling, given the real part of the
   pole mass of the top quark, specified as MTpoletarget. The other
   arguments are as described for the previous function
   SMDR_Eval_Mt_pole().

   If the argument Q_eval is positive, then the inputs are obtained
   by first RG running the MSbar parameter global variables:

   Q_in, lambda_in, g3_in, g_in, gp_in, yb_in, ..., v_in

   to the scale Q_eval, using SMDR_RGeval_SM().

   If the argument Q_eval is negative, then the inputs are instead
   taken to be the current values of the MSbar parameter global
   variables:

   Q, lambda, g3, g, gp, yb, ..., v.
*/
SMDR_REAL SMDR_Eval_yt (SMDR_REAL Q_eval,
                        SMDR_REAL MTpoletarget,
                        int method,
                        int QCDLoopOrder,
                        float otherLoopOrder);

/* Controls the error in the iterated computation of the top pole mass. */
SMDR_REAL SMDR_MTPOLE_TOLERANCE;

/* ---------------------------------------------------------------------- */
/* ----------------------- In GFermi.c ---------------------------------- */

/* Returns GFermi at up to two loops. Based on Deltartilde given in
   section 3 of  arXiv:1907.02500 and the ancillary file Deltartilde.txt.
   See also Deltarbar in
   Kniehl and Veretin 1401.1844 eqs (37)-(40) and
   Kniehl, Pikelner, Veretin 1503.02138 eqs (60)-(62). Appendix A.2
   Kniehl, Pikelner, Veretin 1601.08143, file
   mr-1.3.2/mr/dr/drbar20.cpp in the computer code "mr".
 
  The argument loopOrder has allowed values:
    0    tree-level
    1    1-loop
    1.3  1-loop plus 2-loop leading order in QCD
    1.5  1-loop plus 2-loop leading order in QCD and yt
    2     Full 2-loop result

    If the argument Q_eval is positive, then the inputs are obtained
    by first RG running the MSbar parameters:

    Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

    to the scale Q_eval, using SMDR_RGeval_SM().

    If the argument Q_eval is negative, then the inputs are instead
    given by the current values of the MSbar parameters:

    Q, g3, g, gp, yt, ..., v.
*/
SMDR_REAL SMDR_Eval_GFermi (SMDR_REAL Q_eval, float loopOrder);

/* Returns the result for the Fermi decay constant G_Fermi as found by
   1411.7040 Degrassi, Gambino, and Giardino.  Uses as inputs the
   global variable Standard Model input values

   Q_in, g3_in, g_in, gp_in, yt_in, ...

   as well as the arguments Mtpole, Mhpole, MWpole.  This is an
   alternative to SMDR_Eval_GFermi().
*/
SMDR_REAL SMDR_Eval_GFermi_DGG (SMDR_REAL Mtpole,
                                SMDR_REAL Mhpole,
                                SMDR_REAL MWpole);

/* ---------------------------------------------------------------------- */
/* ----------------------- In Mlight.c ---------------------------------- */

/* Has the same effect as calling all of the following individual
   functions:

   SMDR_Eval_mbmb(); 
   SMDR_Eval_mcmc(); 
   SMDR_Eval_mquarks_2GeV();
   SMDR_Eval_Mtau_pole(); 
   SMDR_Eval_Mmuon_pole(); 
   SMDR_Eval_Melectron_pole();

   However, it is slightly more efficient than calling them
   individually.
*/   
void SMDR_Eval_Light_Masses ();

/* Returns the bottom-quark MSbar running mass evaluated at itself,
   mb(mb), in the 5-quark, 3-charged-lepton QCD+QED effective theory.
   Uses beta functions with loop order loopOrder = 1, 2, 3, 4, 5.  The
   inputs are the global variables corresponding to the MSbar
   parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., v_in.

   Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW.
*/
SMDR_REAL SMDR_Eval_mbmb (SMDR_REAL Q_dec_thZW, int loopOrder);

/* Returns the charm-quark MSbar running mass evaluated at itself,
   mc(mc), in the 4-quark, 2-charged-lepton QCD+QED effective theory.
   Uses beta functions with loop order loopOrder = 1, 2, 3, 4, 5.  The
   inputs are the global variables corresponding to the MSbar
   parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

   Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
   decoupling of the bottom and tau occur at Q_dec_bottom, Q_dec_tau.
*/
SMDR_REAL SMDR_Eval_mcmc (SMDR_REAL Q_dec_tau,
                          SMDR_REAL Q_dec_bottom,
                          SMDR_REAL Q_dec_thZW, int loopOrder);

/* Evaluates the strange, down, and up-quark MSbar running masss
   evaluated at Q=2 GeV, in the 4-quark, 3-charged-lepton QCD+QED
   effective theory.  Uses beta functions with loop order loopOrder =
   1, 2, 3, 4, 5.  The inputs are the global variables corresponding
   to the MSbar parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

   Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
   decoupling of the bottom at Q_dec_bottom.
*/
void SMDR_Eval_mquarks_2GeV (SMDR_REAL Q_dec_bottom,
                             SMDR_REAL Q_dec_thZW, 
                             int loopOrder,
                             SMDR_REAL *ms,
                             SMDR_REAL *mu,
                             SMDR_REAL *md);

/* Returns the pole mass of the bottom quark, computed in the 5-quark,
   3-charged-lepton QCD+QED effective field theory at the scale
   Q_eval, in the approximation of QCD loopOrder = 0, 1, 2, 3, or 4.
   (The RG running is done including all known effects up to 5-loop
   order.)  The starting inputs are the global variables corresponding
   to the MSbar parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

   Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW.

   NOTE: This just doesn't converge fast enough to be useful! 

   The Review of Particle Properties uses the 2-loop approximation,
   which is horrible. See 1606.06754 for discussion of the atrocious
   convergence properties of the expansion. Therefore, the bottom (and
   charm) pole masses are kind of useless as observables and are
   deprecated.
*/
SMDR_REAL SMDR_Eval_Mb_pole (SMDR_REAL Q_eval, 
                             SMDR_REAL Q_dec_thZW, 
                             int loopOrder);

/* Returns the pole mass of the tau lepton, computed in the 4-quark,
   3-charged-lepton QCD+QED effective field theory at the scale
   Q_eval, in the approximation of loopOrder = 0, 1, 2, 3.  (The RG
   running is done including all known effects up to 5-loop order.)
   The starting inputs are the global variables corresponding to the
   MSbar parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

   Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
   decoupling of the bottom quark occurs at the scale Q_dec_bottom.
*/
SMDR_REAL SMDR_Eval_Mtau_pole (SMDR_REAL Q_eval,
                               SMDR_REAL Q_dec_bottom,
                               SMDR_REAL Q_dec_thZW,
                               int loopOrder);

/* Returns the pole mass of the muon, computed in the 3-quark,
   2-charged-lepton QCD+QED effective field theory at the scale
   Q_eval, in the approximation of loopOrder = 0, 1, or 2.  (The RG
   running is done including all known effects up to 5-loop order.)
   The starting inputs are the global variables corresponding to the
   MSbar parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

   Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
   decoupling of the bottom, tau, and charm occur at the scales
   Q_dec_bottom, Q_dec_tau, and Q_dec_charm, respectively.
*/
SMDR_REAL SMDR_Eval_Mmuon_pole (SMDR_REAL Q_eval,
                                SMDR_REAL Q_dec_charm,
                                SMDR_REAL Q_dec_tau,
                                SMDR_REAL Q_dec_bottom,
                                SMDR_REAL Q_dec_thZW,
                                int loopOrder);

/* Returns the pole mass of the electron, computed in the 3-quark,
   2-charged-lepton QCD+QED effective field theory at the scale
   Q_eval, in the approximation of loopOrder = 0, 1, or 2.  (The RG
   running is done including all known effects up to 5-loop order.)
   The starting inputs are the global variables corresponding to the
   MSbar parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

   Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
   decoupling of the bottom, tau, and charm occur at the scales
   Q_dec_bottom, Q_dec_tau, and Q_dec_charm, respectively.
*/
SMDR_REAL SMDR_Eval_Melectron_pole (SMDR_REAL Q_eval,
                                    SMDR_REAL Q_dec_charm,
                                    SMDR_REAL Q_dec_tau,
                                    SMDR_REAL Q_dec_bottom,
                                    SMDR_REAL Q_dec_thZW,
                                    int loopOrder);

/* ---------------------------------------------------------------------- */
/* ------------------------- In arguments.c ----------------------------- */

/* This is a simple function to process a set of command line
   arguments. The arguments and their types are defined in the calling
   program, and passed in via the variables:

   nargs
   -----
   The total number of defined arguments, equal to the sum of the
   number of required arguments and the number of possible optional
   arguments.

   arglist 
   -------
   A list of strings (of dimension nargs) specifying the
   arguments. *Required* arguments must come first, and are given in a
   specified order that cannot be modified. A required argument is
   indicated with "req". After all required args are given, the
   arglist entries should define a set of flags that begin with
   '-'. Each of these specifies an *optional* argument. As an example,
   if there are three required arguments and two optional ones, then
   arglist might look like this:

   char arglist[] = {"req","req","req","-a","-i"};

   In this case the two optional argw would be specified on the
   command line as

   [-a <arg1> [-i <arg2>]]

   Optional arguments can be given in any order on the commad line,
   but the flags specified here must match up with the correspondng
   elements in the following arrays that define them.
   
   argtype
   -------
   An array of strings (of dimension nargs) defining the types of the
   arguments, in order. Allowed type strings are currently "real",
   "int", "string", and "toggle". As an example,

   char argtype[] = {"string","int","real","real","string"};

   Type "toggle" is used for an argument that acts as a toggle, i.e.,
   without an associated value.

   argvar
   ------
   This is an array of pointers-to-void (of dimension nargs), to data
   objects in the calling program to which the read values should be
   assigned. As an example,

   void argvar[] = {inFile, &i, &x, &Q, outFile};

   would assign the first argument in the variable inFile, the second
   in the variable i, and so on. Note that no ampersand is required
   for "string"-type variables.

   For arguments of type "toggle", the assigned variable should
   normally be of type int.  The effect of including the toggle will
   be to set this int to 1 (YES).

   As another example, for three optional arguments corresponding to
   an input filename, and output filename, and a real error tolerance,
   one could specify:

   int nargs = 3;
   char *arglist[] = {"-e","-i","-o"};
   char *argtype[] = {"real","string","string"};
   void *argvar[] = {&ERROR_TOLERANCE, inputFile, outputFile};

   Users should be sure to specify default values for all optional
   arguments in the calling program.

   Some simple tests for errors are included, but these are probably
   not bulletproof. Caveat emptor!
*/
int SMDR_Process_Arguments (int, char *[], int, char *[], char *[], void *[]);

/* ---------------------------------------------------------------------- */
/* --------------------------- In io.c ---------------------------------- */

/* Prints the Standard Model MSbar input parameters and other
   quantities at the current scale Q */
void SMDR_Display_MSbar_Parameters (void);
void SMDR_Display_v (void);
void SMDR_Display_m2 (void);
void SMDR_Display_Lambda (void);
void SMDR_Display_Delta_alpha_had5 (void);

/* Prints the Standard Model on-shell input parameters. */
void SMDR_Display_OS_Inputs (void);

/* The following versions of the above take a file pointer to write to
   and a string to prepend to each line of output. SMDR_Display_v ()
   is equivalent to SMDR_Write_v (stdout, ""). To put "# " at the
   beginning of each line, use SMDR_Write_v (stdout, "# ").
*/
void SMDR_Write_MSbar_Parameters (FILE *, char *);
void SMDR_Write_v (FILE *, char *);
void SMDR_Write_m2 (FILE *, char *);
void SMDR_Write_Lambda (FILE *, char *);
void SMDR_Write_Delta_alpha_had5 (FILE *, char *);

void SMDR_Write_OS_Inputs (FILE *, char *);

/* The following three functions can be used to read data from
   files. The list of pre-defined variables that can be read using
   these utilities can be found in the file "src/smdr_ioparams.h".
*/

/* The following function takes a file stream and a variable name,
   searches the stream for the variable, and sets the corresponding
   SMDR variable to the value found. If the value is not found, a
   warning is issued. Example:

   SMDR_Get_Value (fp, "SMDR_Q_in");
*/
int SMDR_Get_Value (FILE *, char *);

/* The following function takes a file name and a variable name,
   searches the file for the variable, and sets the corresponding SMDR
   variable to the value found. If the value is not found, an error
   results. Example:

   SMDR_Read_Value ("inputFileName", "SMDR_g3_in");

*/
int SMDR_Read_Value (char *, char *);

/* The following function takes a file name and an array of variable
   names, searches the file for the variables, and sets the
   corresponding SMDR variables to the values found. If any value is
   not found, an error results. Example:

   char *toRead[] = {"SMDR_Q_in", SMDR_g3_in"};
   SMDR_Read_Values ("inputFileName", 2, toRead);
*/
int SMDR_Read_Values (char *, int, char *[]);

/* Simple function that takes a number of variables and a list of
   their names, and steps through the list, printing the current value
   of each to stdout and prompting the user for a new value.
 */
int SMDR_Set_Values_Interactively (int, char *[]);

/* Writs a complete model file, including both MSbar and on-shell
   quantities, equivalent to ReferenceModel.dat.
*/
int SMDR_Write_Model_File (char *);

/* Print version information. SMDR_Display_Version () is equivalent to
   SMDR_Write_Version (stdout, "");
*/
void SMDR_Display_Version ();
void SMDR_Write_Version (FILE *, char *);

/* This utility prints a pre-defined array of strings to a file
   stream, with a prepended character.  See the files
   applications/fig_*.c for examples of its use.
*/
void SMDR_Write_Column_Data (FILE *, int, char *[], char *);

/* The following four functions read data from a file. Format should
   be as in ReferenceModel.dat, although the lines can be in any
   order, the trailing semicolons can be included or not, and any line
   starting with a '#' is treated as a comment. Blank lines are
   ignored.
*/
/* Reads the complete set of MSbar inputs: */   
int SMDR_Read_MSbar_Inputs (char *filename);

/* Reads the complete set of on-shell inputs: */   
int SMDR_Read_OS_Inputs (char *filename);

/* Reads all MSbar and on-shell input values: */
int SMDR_Read_Model_File (char *filename);


/* ---------------------------------------------------------------------- */
/* ----------------------- In utilities.c ------------------------------- */

/* Prints the version number, found in smdr_build.h, to stdout. */
void SMDR_Version (void);

/* Loads the Standard Model MSbar parameters from the "static" global
   input variables X_in to the corresponding global variables X, where

   X = Q, g3, g, gp, yt, yb, yc, ys, yu, yd, ytau, ymu, ye, lambda,
       m2, v, Lambda, and Delta_alpha_had_5_MZ.
*/
void SMDR_Load_Inputs (void);

/* Loads the Standard Model MSbar parameters from the current global
   input variables X to the corresponding global variables X_in, where

   X = Q, g3, g, gp, yt, yb, yc, ys, yu, yd, ytau, ymu, ye, lambda,
       m2, v, Lambda, and Delta_alpha_had_5_MZ.

   This function is the inverse of SMDR_Load_Inputs().
*/
void SMDR_Save_Inputs (void);

/* Updates various useful parameter combinations.*/
void SMDR_Update (void);

/* These functions performs sanity checks for the Standard Model MSbar
   parameters g3, g, gp, lambda, yb, ytau, v, Q.
*/
void SMDR_Check_Ranges (void);
void SMDR_Check_Q_Range (SMDR_REAL Q_range_lo, SMDR_REAL Q_range_hi);
void SMDR_Check_VEV_Range (SMDR_REAL vev_range_lo, SMDR_REAL vev_range_hi);
void SMDR_Check_m2_Range (SMDR_REAL m2_range_lo, SMDR_REAL m2_range_hi);
void SMDR_Check_lambda_Range (SMDR_REAL k_range_lo, SMDR_REAL k_range_hi);
void SMDR_Check_yt_Range (SMDR_REAL yt_range_lo, SMDR_REAL yt_range_hi);
void SMDR_Check_yb_Range (SMDR_REAL yb_range_hi);
void SMDR_Check_ytau_Range (SMDR_REAL ytau_range_hi);
void SMDR_Check_g3_Range (SMDR_REAL g3_range_lo, SMDR_REAL g3_range_hi);
void SMDR_Check_g_Range (SMDR_REAL g_range_lo, SMDR_REAL g_range_hi);
void SMDR_Check_gp_Range (SMDR_REAL gp_range_lo, SMDR_REAL gp_range_hi);
void SMDR_Check_Mhpole_Range (SMDR_REAL Mhpole, 
                              SMDR_REAL Mh_range_lo, SMDR_REAL Mh_range_hi);
void SMDR_Check_Mtpole_Range (SMDR_REAL Mtpole, 
                              SMDR_REAL Mt_range_lo, SMDR_REAL Mt_range_hi);

/* Computes the MSbar couplings and masses in the 5-quark, 3-charged-lepton
   QCD+QED effective theory, at the scale Q_MZ, taking t, h, Z, W to be
   decoupled at Q_dec_thZW. The RG running is evaluated at
   loopOrder = 1, 2, 3, 4, or 5.

   The input parameters are the Standard Model MSbar global variables inputs:

   Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in.

   The results are put in the global variables:

   alphaS_5_MZ, alpha_5_MZ, m<fermion>_53, with

   fermion = b, c, s, d, u, tau, muon, electron.
*/
void SMDR_Eval_QCDQED_at_MZ (SMDR_REAL Q_MZ, 
                             SMDR_REAL Q_dec_thZW,
                             int loopOrder);

/* Returns sgn(x) sqrt(|x|) */
SMDR_REAL SMDR_SGNSQRT (SMDR_REAL);

/* Timing stuff: */
void SMDR_Start_Timer (void);
double SMDR_Timer (void);
clock_t SMDR_Time_Start, SMDR_Time_End;
double SMDR_Time_Total;

/* ---------------------------------------------------------------------- */
/* ------------------------ In fit_obs.c -------------------------------- */

/* Finds the (non-decoupled) Standard Model MSbar parameters at the
   scale Q_target. The results are stored in the global variables:

   Q_in=Q_target, g3_in, g_in, gp_in, v_in, lambda_in, y<fermion>_in

   where fermion = t, b, c, s, d, u, tau, muon, electron.

   The inputs are the arguments:

   alphaS_target (in the 5-quark, 3-charged-lepton QCD+QED theory at Q=MZ),  
   alpha_target  (in the 5-quark, 3-charged-lepton theory at Q=MZ),  
   MZ_target (the Z-boson Breit-Wigner mass),
   GFermi_target (the Fermi constant),
   Mh_target (the Higgs boson pole mass),
   Mt_target (the top-quark pole mass),
   mbmb_target (the MSbar bottom mass evaluated at itself in the 5-quark,
               3-charged-lepton QCD+QED theory),
   mcmc_target (the MSbar charm mass evaluated at itself in the 4-quark,
               2-charged-lepton QCD+QED theory),
   ms_2GeV_target (the MSbar strange quark mass at Q=2 GeV, in the 3-quark,
                  2-charged-lepton QCD+QED theory),
   md_2GeV_target (the MSbar down quark mass at Q=2 GeV, in the 3-quark,
                  2-charged-lepton QCD+QED theory),
   mu_2GeV_target (the MSbar up quark mass at Q=2 GeV, in the 3-quark,
                  2-charged-lepton QCD+QED theory),
   Mtau_target (the tau lepton pole mass),
   Mmuon_target (the muon pole mass),
   Melectron_target (the electron pole mass),
   Delta_alpha_had (the non-perturbative hadronic contribution to the 
                   fine-structure constant at MZ) 

   The calculation proceeds by iteration, until the largest of the
   fractional errors in the target quantities is less than
   error_target. About 5 or fewer iterations are typically expected
   for error_target = 10^-8.
*/
int SMDR_Fit_Inputs (SMDR_REAL Q_target,
                     SMDR_REAL alphaS_target,
                     SMDR_REAL alpha_target,
                     SMDR_REAL GFermi_target,
                     SMDR_REAL MZ_target,
                     SMDR_REAL Mh_target,
                     SMDR_REAL Mt_target,
                     SMDR_REAL mbmb_target,
                     SMDR_REAL mcmc_target,
                     SMDR_REAL ms_2GeV_target,
                     SMDR_REAL md_2GeV_target,
                     SMDR_REAL mu_2GeV_target,
                     SMDR_REAL Mtau_target,    
                     SMDR_REAL Mmuon_target,
                     SMDR_REAL Melectron_target,
                     SMDR_REAL Delta_alpha_had,
                     SMDR_REAL error_target);

/* Evaluates coefficients that allow for accelerated convergence of
   the iteration process in SMDR_Fit_Inputs(), and stores them in the
   file accelcoeffs.h. This function is called by the command line
   program

   ./make_coeffs

   which should be run if/whenever the reference model in
   ReferenceModel.c changes.
*/
void SMDR_Make_Accelcoeffs ();

/* messages.c: */
void SMDR_Error (char *, char *, int);
void SMDR_Warn (char *, char *);

/* arguments.c: */
int SMDR_Process_Arguments (int, char *[], int, char *[], char *[], void *[]);
int SMDR_Write_Output (char *, char *);

#ifdef __cplusplus
}
#endif

#endif /* smdr.h */
