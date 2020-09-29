/* General header file for internal use.  Users should include smdr.h
   instead! */

#ifndef _SMDR_INTERNAL_H_
#define _SMDR_INTERNAL_H_

/* First, a bunch of C standard library stuff: */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

/* Then the local headers: */
#include "smdr.h"        /* This also includes tsil.h */
#include "smdr_build.h"  /* Version information */
#include "tsil.h"
#include "3vil.h"

/* MSbar input parameters (initial values) */
#define Q_in    SMDR_Q_in
#define g3_in   SMDR_g3_in
#define gp_in   SMDR_gp_in
#define g_in    SMDR_g_in
#define yt_in   SMDR_yt_in
#define yb_in   SMDR_yb_in
#define yc_in   SMDR_yc_in
#define ys_in   SMDR_ys_in
#define yu_in   SMDR_yu_in
#define yd_in   SMDR_yd_in
#define ytau_in SMDR_ytau_in
#define ymu_in  SMDR_ymu_in
#define ye_in   SMDR_ye_in
#define k_in    SMDR_lambda_in
#define m2_in   SMDR_m2_in
#define v_in    SMDR_v_in
#define Lambda_in    SMDR_Lambda_in
#define Delta_alpha_had_5_MZ_in SMDR_Delta_alpha_had_5_MZ_in

/* MSbar input parameters (variable values; obtained by running from initial
                           values) */
#define Q    SMDR_Q
#define g3   SMDR_g3
#define gp   SMDR_gp
#define g    SMDR_g
#define yt   SMDR_yt
#define yb   SMDR_yb
#define yc   SMDR_yc
#define ys   SMDR_ys
#define yu   SMDR_yu
#define yd   SMDR_yd
#define ytau SMDR_ytau
#define ymu  SMDR_ymu
#define ye   SMDR_ye
#define k    SMDR_lambda
#define m2   SMDR_m2
#define v    SMDR_v
#define Lambda    SMDR_Lambda
#define Delta_alpha_had_5_MZ SMDR_Delta_alpha_had_5_MZ

/* Useful definitions of various combinations of parameters. */
#define v2 SMDR_v2
#define v4 SMDR_v4
#define Q2 SMDR_Q2
#define g32 SMDR_g32
#define g34 SMDR_g34
#define g36 SMDR_g36
#define g38 SMDR_g38
#define g2 SMDR_g2
#define g4 SMDR_g4
#define g6 SMDR_g6
#define g8 SMDR_g8
#define g10 SMDR_g10
#define g12 SMDR_g12
#define g14 SMDR_g14
#define g16 SMDR_g16
#define g18 SMDR_g18
#define g20 SMDR_g20
#define gp2 SMDR_gp2
#define gp4 SMDR_gp4
#define gp6 SMDR_gp6
#define gp8 SMDR_gp8
#define gp10 SMDR_gp10
#define gp12 SMDR_gp12
#define gp14 SMDR_gp14
#define gp16 SMDR_gp16
#define gp18 SMDR_gp18
#define yt2 SMDR_yt2
#define yt4 SMDR_yt4
#define yt6 SMDR_yt6
#define yt8 SMDR_yt8
#define yt10 SMDR_yt10
#define yt12 SMDR_yt12
#define yb2 SMDR_yb2
#define yb4 SMDR_yb4
#define yb6 SMDR_yb6
#define yb8 SMDR_yb8
#define ytau2 SMDR_ytau2
#define ytau4 SMDR_ytau4
#define ytau6 SMDR_ytau6
#define ytau8 SMDR_ytau8
#define k2 SMDR_k2
#define k3 SMDR_k3
#define k4 SMDR_k4
#define k5 SMDR_k5
#define k6 SMDR_k6
#define k7 SMDR_k7
#define g2pgp2 SMDR_g2pgp2
#define auL SMDR_auL
#define auR SMDR_auR
#define adL SMDR_adL
#define adR SMDR_adR
#define aeL SMDR_aeL
#define aeR SMDR_aeR
#define anL SMDR_anL
#define auL2 SMDR_auL2
#define auR2 SMDR_auR2
#define adL2 SMDR_adL2
#define adR2 SMDR_adR2
#define aeL2 SMDR_aeL2
#define aeR2 SMDR_aeR2
#define anL2 SMDR_anL2
#define sqrtg2pgp2 SMDR_sqrtg2pgp2
#define g2pgp2 SMDR_g2pgp2
#define g2pgp22 SMDR_g2pgp22
#define g2pgp23 SMDR_g2pgp23
#define g2pgp24 SMDR_g2pgp24
#define g2pgp25 SMDR_g2pgp25
#define g2m2k SMDR_g2m2k
#define g2m2k2 SMDR_g2m2k2
#define g2m4k SMDR_g2m4k
#define g2m4k2 SMDR_g2m4k2
#define g2m8k SMDR_g2m8k
#define g2m8k2 SMDR_g2m8k2
#define g2m8k3 SMDR_g2m8k3
#define g2p8k SMDR_g2p8k
#define g2p8k2 SMDR_g2p8k2
#define g2mgp2 SMDR_g2mgp2
#define g2mgp22 SMDR_g2mgp22
#define g2m2yt2 SMDR_g2m2yt2
#define g2m2yt22 SMDR_g2m2yt22
#define g2pgp2m2k SMDR_g2pgp2m2k
#define g2pgp2m2k2 SMDR_g2pgp2m2k2
#define g2pgp2m4yt2 SMDR_g2pgp2m4yt2
#define g2pgp2m4yt22 SMDR_g2pgp2m4yt22
#define g2pgp2m8k SMDR_g2pgp2m8k
#define g2pgp2m8k2 SMDR_g2pgp2m8k2
#define g2pgp2m8k3 SMDR_g2pgp2m8k3
#define g4pg2gp2mgp4 SMDR_g4pg2gp2mgp4
#define g4pg2gp2mgp42 SMDR_g4pg2gp2mgp42
#define twog2pgp2 SMDR_twog2pgp2
#define twog2pgp22 SMDR_twog2pgp22
#define twog2pgp23 SMDR_twog2pgp23
#define threeg2m5gp2 SMDR_threeg2m5gp2
#define threeg2m5gp22 SMDR_threeg2m5gp22
#define threeg2mgp2 SMDR_threeg2mgp2
#define threeg2mgp22 SMDR_threeg2mgp22
#define threeg2pgp2 SMDR_threeg2pgp2
#define threeg2pgp22 SMDR_threeg2pgp22
#define T SMDR_T
#define W SMDR_W
#define Z SMDR_Z
#define G SMDR_G
#define H SMDR_H
#define h SMDR_h
#define b SMDR_b
#define tau SMDR_tau
#define h2 SMDR_h2
#define h3 SMDR_h3
#define h4 SMDR_h4
#define h5 SMDR_h5
#define W2 SMDR_W2
#define W3 SMDR_W3
#define W4 SMDR_W4
#define W5 SMDR_W5
#define W6 SMDR_W6
#define W7 SMDR_W7
#define Z2 SMDR_Z2
#define Z3 SMDR_Z3
#define Z4 SMDR_Z4
#define Z5 SMDR_Z5
#define Z6 SMDR_Z6
#define T2 SMDR_T2
#define T3 SMDR_T3
#define T4 SMDR_T4
#define b2 SMDR_b2
#define tau2 SMDR_tau2

#define beta_gp SMDR_beta_gp
#define beta_g SMDR_beta_g
#define beta_g3 SMDR_beta_g3
#define beta_yt SMDR_beta_yt
#define beta_yb SMDR_beta_yb
#define beta_ytau SMDR_beta_ytau
#define beta_yc SMDR_beta_yc
#define beta_yu SMDR_beta_yu
#define beta_ys SMDR_beta_ys
#define beta_yd SMDR_beta_yd
#define beta_ymu SMDR_beta_ymu
#define beta_ye SMDR_beta_ye
#define beta_k SMDR_beta_k
#define beta_m2 SMDR_beta_m2
#define beta_Lambda SMDR_beta_Lambda
#define beta_v SMDR_beta_v

#define Q_53 SMDR_Q_53
#define alphaS_53 SMDR_alphaS_53
#define alpha_53 SMDR_alpha_53
#define mb_53 SMDR_mb_53
#define mc_53 SMDR_mc_53
#define ms_53 SMDR_ms_53
#define mu_53 SMDR_mu_53
#define md_53 SMDR_md_53
#define mtau_53 SMDR_mtau_53
#define mmuon_53 SMDR_mmuon_53
#define melectron_53 SMDR_melectron_53

#define Q_43 SMDR_Q_43
#define alphaS_43 SMDR_alphaS_43
#define alpha_43 SMDR_alpha_43
#define mc_43 SMDR_mc_43
#define ms_43 SMDR_ms_43
#define mu_43 SMDR_mu_43
#define md_43 SMDR_md_43
#define mtau_43 SMDR_mtau_43
#define mmuon_43 SMDR_mmuon_43
#define melectron_43 SMDR_melectron_43

#define Q_42 SMDR_Q_42
#define alphaS_42 SMDR_alphaS_42
#define alpha_42 SMDR_alpha_42
#define mc_42 SMDR_mc_42
#define ms_42 SMDR_ms_42
#define mu_42 SMDR_mu_42
#define md_42 SMDR_md_42
#define mmuon_42 SMDR_mmuon_42
#define melectron_42 SMDR_melectron_42

#define Q_32 SMDR_Q_32
#define alphaS_32 SMDR_alphaS_32
#define alpha_32 SMDR_alpha_32
#define ms_32 SMDR_ms_32
#define mu_32 SMDR_mu_32
#define md_32 SMDR_md_32
#define mmuon_32 SMDR_mmuon_32
#define melectron_32 SMDR_melectron_32

/* Used by SMDR_RGrun_QCDQED and SMDR_RG_rk6_QCDQED.  */
#define Q_eff SMDR_Q_eff
#define alphaS_run SMDR_alphaS_run
#define alpha_run SMDR_alpha_run
#define lnmu SMDR_lnmu
#define lnmd SMDR_lnmd
#define lnme SMDR_lnme

#define beta_alphaS SMDR_beta_alphaS
#define beta_alpha SMDR_beta_alpha
#define beta_lnmu SMDR_beta_lnmu
#define beta_lnmd SMDR_beta_lnmd
#define beta_lnme SMDR_beta_lnme

/* Some useful combinations: */
SMDR_REAL v2, v4, Q2;
SMDR_REAL g32, g34, g36, g38;
SMDR_REAL g2, g4, g6, g8, g10, g12, g14, g16, g18, g20;
SMDR_REAL gp2, gp4, gp6, gp8, gp10, gp12, gp14, gp16, gp18;
SMDR_REAL yt2, yt4, yt6, yt8, yt10, yt12;
SMDR_REAL yb2, yb4, yb6, yb8;
SMDR_REAL ytau2, ytau4, ytau6, ytau8;
SMDR_REAL k2, k3, k4, k5, k6, k7;
SMDR_REAL g2pgp2;
SMDR_REAL auL,auR,adL,adR,aeL,aeR,anL;
SMDR_REAL auL2,auR2,adL2,adR2,aeL2,aeR2,anL2;

/* More useful combinations: */
SMDR_REAL sqrtg2pgp2, g2pgp2, g2pgp22, g2pgp23, g2pgp24, g2pgp25;
SMDR_REAL g2m2k, g2m2k2;
SMDR_REAL g2m4k, g2m4k2;
SMDR_REAL g2m8k, g2m8k2, g2m8k3;
SMDR_REAL g2p8k, g2p8k2;
SMDR_REAL g2mgp2, g2mgp22;
SMDR_REAL g2m2yt2, g2m2yt22;
SMDR_REAL g2pgp2m2k, g2pgp2m2k2;
SMDR_REAL g2pgp2m4yt2, g2pgp2m4yt22;
SMDR_REAL g2pgp2m8k, g2pgp2m8k2, g2pgp2m8k3;
SMDR_REAL g4pg2gp2mgp4, g4pg2gp2mgp42;
SMDR_REAL twog2pgp2, twog2pgp22, twog2pgp23;
SMDR_REAL threeg2m5gp2, threeg2m5gp22;
SMDR_REAL threeg2mgp2, threeg2mgp22;
SMDR_REAL threeg2pgp2, threeg2pgp22;

/* Tree-level running squared masses: */
SMDR_REAL T, W, Z, G, H, h, b, tau;

/* Useful running squared mass powers: */
SMDR_REAL h2, h3, h4, h5;
SMDR_REAL W2, W3, W4, W5, W6, W7;
SMDR_REAL Z2, Z3, Z4, Z5, Z6;
SMDR_REAL T2, T3, T4;
SMDR_REAL b2, tau2;

/* Beta functions for full Standard Model: */
SMDR_REAL beta_gp, beta_g, beta_g3, beta_yt, beta_yb, beta_ytau;
SMDR_REAL beta_yc, beta_yu, beta_ys, beta_yd, beta_ymu, beta_ye;
SMDR_REAL beta_k, beta_m2, beta_Lambda, beta_v;

/* Running quantities for QCDQED theory with t,h,Z,W decoupled: */
SMDR_REAL Q_eff, alphaS_run, alpha_run, lnmu, lnmd, lnme;

SMDR_REAL Q_53, alphaS_53, alpha_53;
SMDR_REAL mb_53, mc_53, ms_53, mu_53, md_53;
SMDR_REAL mtau_53, mmuon_53, melectron_53;

SMDR_REAL Q_43, alphaS_43, alpha_43;
SMDR_REAL mc_43, ms_43, mu_43, md_43;
SMDR_REAL mtau_43, mmuon_43, melectron_43;

SMDR_REAL Q_42, alphaS_42, alpha_42;
SMDR_REAL mc_42, ms_42, mu_42, md_42;
SMDR_REAL mmuon_42, melectron_42;

SMDR_REAL Q_32, alphaS_32, alpha_32;
SMDR_REAL ms_32, mu_32, md_32;
SMDR_REAL mmuon_32, melectron_32;

/* Beta functions for QCDQED theory with t,h,Z,W decoupled: */
SMDR_REAL beta_alphaS, beta_alpha, beta_lnmu, beta_lnmd, beta_lnme;

/* From betas.c */
int SMDR_Betas (int loopOrder);  

/* From utilities.c */
void SMDR_update_for_betas (void);   
SMDR_REAL SMDR_SGNSQRT (SMDR_REAL x); 
int SMDR_FPCompare (SMDR_REAL x, SMDR_REAL y); 

/* From rk6_RG.c */
int SMDR_RG_rk6 (SMDR_REAL *dt, int force_step, int loopOrder); 
int SMDR_RG_rk6_QCDQED (SMDR_REAL *dt, int force_step, int loopOrder,
                        int nu, int nd, int ne);

/* From betas_QCDQED.c */
void SMDR_Betas_QCDQED (int loopOrder, int nu, int nd, int ne);

/* From Mlight.c */
SMDR_REAL SMDR_f2lf (SMDR_REAL x);

/* io.c: */
int GetWord (FILE *, char *);
int SkipRestOfLine (FILE *);
int SkipWhiteSpace (FILE *);

/* Miscellaneous constants: */
#ifndef PI
#define PI     3.1415926535897932385L
#endif
#define PI2    9.8696044010893586188L
#define ln2    0.69314718055994530942L
#define Zeta2  1.6449340668482264365L
#define Zeta3  1.2020569031595942854L
#define SQRT2  1.41421356237309504880169L

/* 1/(16 pi^2) */
#define ONELOOPFACTOR (0.00633257397764611071524247L) 

/* 1/(16 pi^2)^2 */
#define TWOLOOPFACTOR (0.0000401014931823606843326281L)

/* 1/(16 pi^2)^3 */
#define THREELOOPFACTOR (0.0000002539456721913701894751L)

/* 1/(16 pi^2)^4 */
#define FOURLOOPFACTOR (0.0000000016081297554549204457L)

/* Possibly useful enums: */
enum {FALSE, TRUE};
enum {NO, YES};

#endif /* internal.h */
