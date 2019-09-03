/* 
   This file contains the list of parameters that may be read or
   written using the SMDR IO routines.

   Add any new options here to insure that they are available.
*/

#ifndef _IOPARAMS_H
#define _IOPARAMS_H

int nVals = 35;

char *varName[] = {
  "SMDR_Q_in",
  "SMDR_g3_in",
  "SMDR_gp_in",
  "SMDR_g_in",
  "SMDR_yt_in",
  "SMDR_yb_in",
  "SMDR_yc_in",
  "SMDR_ys_in",
  "SMDR_yu_in",
  "SMDR_yd_in",
  "SMDR_ytau_in",
  "SMDR_ymu_in",
  "SMDR_ye_in",
  "SMDR_lambda_in",
  "SMDR_m2_in",
  "SMDR_v_in",
  "SMDR_Lambda_in",
  "SMDR_Delta_alpha_had_5_MZ_in",
  "SMDR_Mt_pole",
  "SMDR_Mh_pole",
  "SMDR_MZ_BreitWigner",
  "SMDR_MZ_pole",
  "SMDR_MW_BreitWigner",
  "SMDR_MW_pole",
  "SMDR_GFermi",
  "SMDR_alpha",
  "SMDR_alphaS_5_MZ",
  "SMDR_mbmb",
  "SMDR_mcmc",
  "SMDR_ms_2GeV",
  "SMDR_md_2GeV",
  "SMDR_mu_2GeV",
  "SMDR_Mtau_pole",
  "SMDR_Mmuon_pole",
  "SMDR_Melectron_pole"};

SMDR_REAL *varValue[] = {
  &SMDR_Q_in,
  &SMDR_g3_in,
  &SMDR_gp_in,
  &SMDR_g_in,
  &SMDR_yt_in,
  &SMDR_yb_in,
  &SMDR_yc_in,
  &SMDR_ys_in,
  &SMDR_yu_in,
  &SMDR_yd_in,
  &SMDR_ytau_in,
  &SMDR_ymu_in,
  &SMDR_ye_in,
  &SMDR_lambda_in,
  &SMDR_m2_in,
  &SMDR_v_in,
  &SMDR_Lambda_in,
  &SMDR_Delta_alpha_had_5_MZ_in,
  &SMDR_Mt_pole,
  &SMDR_Mh_pole,
  &SMDR_MZ_BreitWigner,
  &SMDR_MZ_pole,
  &SMDR_MW_BreitWigner,
  &SMDR_MW_pole,
  &SMDR_GFermi,
  &SMDR_alpha,
  &SMDR_alphaS_5_MZ,
  &SMDR_mbmb,
  &SMDR_mcmc,
  &SMDR_ms_2GeV,
  &SMDR_md_2GeV,
  &SMDR_mu_2GeV,
  &SMDR_Mtau_pole,
  &SMDR_Mmuon_pole,
  &SMDR_Melectron_pole};

#endif
