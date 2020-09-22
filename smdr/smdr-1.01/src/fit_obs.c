#include "smdr_internal.h"

#define ZEROSAFE(a) (((a) > (SMDR_TOL)) ? (a) : (SMDR_TOL))

/* --------------------------------------------------------------------- */
/* Takes the on-shell observables as arguments, and finds the best fit
   MSbar parameters at the MSbar renormalization scale Q_target. 
   The resulting MSbar quantities are set to be SMDR_Q = SMDR_Q_in,
   SMDR_g3 = SMDR_g3_in, etc.
*/
int  SMDR_Fit_Inputs (SMDR_REAL Q_target,
                      SMDR_REAL alphaS_5_MZ_target,
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
                      SMDR_REAL Delta_alpha_had_target,
                      SMDR_REAL error_target)
{
  char funcname[] = "SMDR_Fit_Inputs";
  int j;
  SMDR_REAL vratio, v2ratio;
  SMDR_REAL KZ, Ka;
  SMDR_REAL err_total, err_alphaS, err_alpha, err_MZ, err_GFermi, err_Mh; 
  SMDR_REAL err_Mt, err_mbmb, err_mcmc, err_ms_2GeV, err_md_2GeV;
  SMDR_REAL err_mu_2GeV, err_Mtau, err_Mmuon, err_Melectron;  
  SMDR_REAL g3rat, lambdarat, ytrat, ybrat, ycrat, ysrat, ydrat, yurat;
  SMDR_REAL ytaurat, ymurat, yerat;
  SMDR_REAL del_v, del_g3, del_g, del_gp, del_yt, del_yb, del_yc, del_lambda;
  SMDR_REAL del_ys, del_yu, del_yd, del_ytau, del_ymu, del_ye;
  int result_niters;
  #include "accelcoeffs.h"

  if (mbmb_target < 3.0) {
    SMDR_Warn (funcname, "Target mb(mb) is too small.");
    SMDR_Warn (funcname, "Setting it and all lighter fermion masses to 0.\n");
    mbmb_target = 0;
    mcmc_target = 0;
    ms_2GeV_target = 0;
    md_2GeV_target = 0;
    mu_2GeV_target = 0;
    Mtau_target = 0;
    Mmuon_target = 0;
    Melectron_target = 0;
  } else if (mcmc_target < 1.0) {
    SMDR_Warn (funcname, "Target mc(mc) is too small for evaluation in QCD.");
    SMDR_Warn (funcname, "Setting it and all lighter fermion masses to 0.\n");
    mcmc_target = 0;
    ms_2GeV_target = 0;
    md_2GeV_target = 0;
    mu_2GeV_target = 0;
    Mmuon_target = 0;
    Melectron_target = 0;
  }

  /* Don't compute Mt, Mh with needlessly high tolerances... */ 
  SMDR_MTPOLE_TOLERANCE = 1.0e-8;
  SMDR_MHPOLE_TOLERANCE = 1.0e-8;
  if (error_target < 1.0e-7) {
    SMDR_MTPOLE_TOLERANCE = error_target/10.;
    SMDR_MHPOLE_TOLERANCE = error_target/10.;
  }
  if (SMDR_MTPOLE_TOLERANCE < 1.01e-15) SMDR_MTPOLE_TOLERANCE = 1.01e-15;
  if (SMDR_MHPOLE_TOLERANCE < 1.01e-15) SMDR_MHPOLE_TOLERANCE = 1.01e-15;

  /* Use ReferenceModel.dat as our initial guess for the MSbar
     parameters, and their corresponding OS predictions.
  */
  /* SMDR_Read_MSbar_Inputs ("ReferenceModel.dat"); */
  /* SMDR_Read_Value ("ReferenceModel.dat", "SMDR_v_in"); */
  /* SMDR_Read_OS_Inputs ("ReferenceModel.dat"); */
  #include "AUTOrefmodel.h"

  SMDR_RGeval_SM (Mt_target, 5);

  SMDR_Save_Inputs ();
  SMDR_Delta_alpha_had_5_MZ_in = Delta_alpha_had_target;

  for (j=0; j<16; j++) {

    /* Find the new input parameter estimates by comparing OS observables to the 
       targets. On the 0th iteration, the OS observables and MSbar parameters
       are both as specified in "ReferenceModel.dat".
    */
    v2ratio = SMDR_GFermi/GFermi_target;
    vratio = SMDR_SQRT (v2ratio);

    KZ = (SMDR_g_in * SMDR_g_in + SMDR_gp_in * SMDR_gp_in) * 
         (MZ_target * MZ_target)/
         (SMDR_MZ_BreitWigner * SMDR_MZ_BreitWigner * v2ratio);

    Ka = ((SMDR_g_in * SMDR_g_in * SMDR_gp_in * SMDR_gp_in)/
          (SMDR_g_in * SMDR_g_in + SMDR_gp_in * SMDR_gp_in)) * 
          (alpha_target/SMDR_alpha);

    del_g  = (SMDR_SQRT(KZ * (1. + SMDR_SQRT(1. - 4. * Ka/KZ))/2.)/g_in) - 1.;
    del_gp = (SMDR_SQRT(KZ * (1. - SMDR_SQRT(1. - 4. * Ka/KZ))/2.)/gp_in) - 1.;

    del_g3 = SMDR_SQRT (alphaS_5_MZ_target/SMDR_alphaS_5_MZ) - 1.;
    del_v = vratio - 1.;
    del_lambda = ((Mh_target * Mh_target)/
                 (SMDR_Mh_pole * SMDR_Mh_pole * v2ratio)) - 1.;
    del_yt = (Mt_target/(SMDR_Mt_pole * vratio)) - 1.;
    del_yb = (mbmb_target/(ZEROSAFE(SMDR_mbmb) * vratio)) - 1.;
    del_yc = (mcmc_target/(ZEROSAFE(SMDR_mcmc) * vratio)) - 1.;
    del_ys = (ms_2GeV_target/(ZEROSAFE(SMDR_ms_2GeV) * vratio)) - 1.;
    del_yd = (md_2GeV_target/(ZEROSAFE(SMDR_md_2GeV) * vratio)) - 1.;
    del_yu = (mu_2GeV_target/(ZEROSAFE(SMDR_mu_2GeV) * vratio)) - 1.;

    del_ytau = (Mtau_target/(ZEROSAFE(SMDR_Mtau_pole) * vratio)) - 1.;
    del_ymu = (Mmuon_target/(ZEROSAFE(SMDR_Mmuon_pole) * vratio)) - 1.;
    del_ye = (Melectron_target/(ZEROSAFE(SMDR_Melectron_pole) * vratio)) - 1.;

/*
    printf("del_v  = %.16Lf\n",del_v);
    printf("del_k  = %.16Lf\n",del_lambda);
    printf("del_g  = %.16Lf\n",del_g);
    printf("del_gp = %.16Lf\n",del_gp);
    printf("del_g3 = %.16Lf\n",del_g3);
    printf("del_yt = %.16Lf\n",del_yt);
    printf("del_yb = %.16Lf\n",del_yb);
    printf("del_yc = %.16Lf\n",del_yc);
    printf("del_ys = %.16Lf\n",del_ys);
    printf("del_yd = %.16Lf\n",del_yd);
    printf("del_yu = %.16Lf\n",del_yu);
    printf("del_ytau = %.16Lf\n",del_ytau);
    printf("del_ymu  = %.16Lf\n",del_ymu);
    printf("del_ye   = %.16Lf\n",del_ye); 
*/

    g3rat = 1. + cog3_g3 * del_g3 
               + cog3_g32 * del_g3 * del_g3 
               + cog3_yt * del_yt;

    lambdarat = 1. + del_lambda 
                   + colambda_yt * del_yt 
                   + colambda_g3 * del_g3;

    ytrat = 1. + coyt_yt * del_yt 
               + coyt_yt2 * del_yt * del_yt
               + coyt_g3 * del_g3;

    ybrat = 1. + coyb_yb * del_yb 
               + coyb_yb2 * del_yb * del_yb +
               + coyb_g3 * del_g3
               + coyb_yt * del_yt;

    ycrat = 1. + coyc_yc * del_yc 
               + coyc_yc2 * del_yc * del_yc 
               + coyc_g3 * del_g3 
               + coyc_yt * del_yt
               + coyc_yb * del_yb;

    ysrat = 1. + del_ys 
               + coys_g3 * del_g3 
               + coys_yt * del_yt 
               + coys_yb * del_yb;

    ydrat = 1. + del_yd 
               + coyd_g3 * del_g3 
               + coyd_yt * del_yt 
               + coyd_yb * del_yb;

    yurat = 1. + del_yu 
               + coyu_g3 * del_g3 
               + coyu_yt * del_yt 
               + coyu_yb * del_yb;

    ytaurat = 1. + del_ytau;
    ymurat = 1. + del_ymu;
    yerat = 1. + del_ye;

    if ((ybrat < 0) || (mbmb_target < SMDR_TOL)) ybrat = 0;
    if ((ycrat < 0) || (mcmc_target < SMDR_TOL)) ycrat = 0;
    if ((ysrat < 0) || (ms_2GeV_target < SMDR_TOL)) ysrat = 0;
    if ((ydrat < 0) || (md_2GeV_target < SMDR_TOL)) ydrat = 0;
    if ((yurat < 0) || (mu_2GeV_target < SMDR_TOL)) yurat = 0;
    if ((ytaurat < 0) || (Mtau_target < SMDR_TOL)) ytaurat = 0;
    if ((ymurat < 0) || (Mmuon_target < SMDR_TOL)) ymurat = 0;
    if ((yerat < 0) || (Melectron_target < SMDR_TOL)) yerat = 0;

    SMDR_v_in *= vratio;
    SMDR_g3_in *= g3rat;
    SMDR_g_in  *= 1. + del_g;
    SMDR_gp_in *= 1. + del_gp;
    SMDR_lambda_in *= lambdarat;
    SMDR_yt_in *= ytrat;
    SMDR_yb_in *= ybrat;
    SMDR_yc_in *= ycrat;
    SMDR_ys_in *= ysrat;
    SMDR_yd_in *= ydrat;
    SMDR_yu_in *= yurat;
    SMDR_ytau_in *= ytaurat;
    SMDR_ymu_in *= ymurat;
    SMDR_ye_in *= yerat;

    err_alphaS = SMDR_alphaS_5_MZ/alphaS_5_MZ_target - 1.;
    err_alpha = SMDR_alpha/alpha_target - 1.;
    err_MZ = SMDR_MZ_BreitWigner/MZ_target - 1.;
    err_GFermi = SMDR_GFermi/GFermi_target - 1.;
    err_Mh = SMDR_Mh_pole/Mh_target - 1.;
    err_Mt = SMDR_Mt_pole/Mt_target - 1.;
    err_mbmb = (SMDR_mbmb - mbmb_target)/SMDR_mbmb_EXPT;
    err_mcmc = (SMDR_mcmc - mcmc_target)/SMDR_mcmc_EXPT;
    err_ms_2GeV = (SMDR_ms_2GeV - ms_2GeV_target)/SMDR_ms_2GeV_EXPT;
    err_md_2GeV = (SMDR_md_2GeV - md_2GeV_target)/SMDR_md_2GeV_EXPT;
    err_mu_2GeV = (SMDR_mu_2GeV - mu_2GeV_target)/SMDR_mu_2GeV_EXPT;
    err_Mtau = (SMDR_Mtau_pole - Mtau_target)/SMDR_Mtau_EXPT;
    err_Mmuon = (SMDR_Mmuon_pole - Mmuon_target)/SMDR_Mmuon_EXPT;
    err_Melectron = (SMDR_Melectron_pole - Melectron_target)/
                    SMDR_Melectron_EXPT;

/*
    printf("err_mbmb = %.16Lf\n",err_mbmb);
    printf("err_mcmc = %.16Lf\n",err_mcmc);
    printf("err_ms_2GeV = %.16Lf\n",err_ms_2GeV);
    printf("err_md_2GeV = %.16Lf\n",err_md_2GeV);
    printf("err_mu_2GeV = %.16Lf\n",err_mu_2GeV);
    printf("err_Mtau = %.16Lf\n",err_Mtau);
    printf("err_Mmuon  = %.16Lf\n",err_Mmuon);
    printf("err_Melectron   = %.16Lf\n",err_Melectron); 
*/

    err_total = SMDR_FABS(err_alphaS) +
                SMDR_FABS(err_alpha) + 
                SMDR_FABS(err_MZ) + 
                SMDR_FABS(err_GFermi) + 
                SMDR_FABS(err_Mh) + 
                SMDR_FABS(err_Mt) + 
                SMDR_FABS(err_mbmb) + 
                SMDR_FABS(err_mcmc) + 
                SMDR_FABS(err_ms_2GeV) + 
                SMDR_FABS(err_md_2GeV) + 
                SMDR_FABS(err_mu_2GeV) + 
                SMDR_FABS(err_Mtau) + 
                SMDR_FABS(err_Mmuon) + 
                SMDR_FABS(err_Melectron);

/*
    printf("SMDR_Q_in  = %.8Lf;\n", SMDR_Q_in);
    printf("SMDR_g3_in = %.16Lf;\n", SMDR_g3_in);
    printf("SMDR_g_in  = %.16Lf;\n", SMDR_g_in);
    printf("SMDR_gp_in = %.16Lf;\n", SMDR_gp_in);
    printf("SMDR_v_in  = %.14Lf;\n", SMDR_v_in);
    printf("SMDR_lambda_in  = %.16Lf;\n", SMDR_lambda_in);
    printf("SMDR_yt_in = %.16Lf;\n", SMDR_yt_in);
    printf("SMDR_yb_in = %.16Lf;\n", SMDR_yb_in);
    printf("SMDR_yc_in = %.18Lf;\n", SMDR_yc_in);
    printf("SMDR_ys_in = %.18Lf;\n", SMDR_ys_in);
    printf("SMDR_yd_in = %.20Lf;\n", SMDR_yd_in);
    printf("SMDR_yu_in = %.20Lf;\n", SMDR_yu_in);
    printf("SMDR_ytau_in = %.16Lf;\n", SMDR_ytau_in);
    printf("SMDR_ymu_in  = %.17Lf;\n", SMDR_ymu_in);
    printf("SMDR_ye_in   = %.20Lf;\n", SMDR_ye_in);
    printf("SMDR_Delta_alpha_had_5_MZ_in = %.12Lf;\n", SMDR_Delta_alpha_had_5_MZ_in);
*/

    if (j>0) {
      printf("Iteration %d: ",j);
      printf("Total fractional error = %.15Lf\n", err_total); 
      if (err_total < error_target) break;
    }
  
    /* Now evaluate all of the on-shell observables. */
    SMDR_Eval_Mt_pole (SMDR_Mt_EXPT, 1, 4, 2, &SMDR_Mt_pole, &SMDR_Gammat_pole);
    SMDR_Eval_Mh_pole (160., 2.5, &SMDR_Mh_pole, &SMDR_Gammah_pole);
    SMDR_Eval_MZ_pole (160., 2, &SMDR_MZ_pole, &SMDR_GammaZ_pole,
                       &SMDR_MZ_BreitWigner, &SMDR_GammaZ_BreitWigner);

    /* This one is computed only because it is needed by SMDR_Eval_Gauge() */
    SMDR_Eval_MW_pole (160., 2, &SMDR_MW_pole, &SMDR_GammaW_pole,
                       &SMDR_MW_BreitWigner, &SMDR_GammaW_BreitWigner);

    SMDR_GFermi = SMDR_Eval_GFermi (SMDR_Mt_EXPT, 2);
    SMDR_Eval_Gauge (SMDR_Mt_pole, SMDR_Mh_pole, SMDR_MW_BreitWigner);
    SMDR_Eval_QCDQED_at_MZ (SMDR_MZ_EXPT, SMDR_MZ_EXPT, 5);

    if (mbmb_target > SMDR_TOL) 
      SMDR_mbmb = SMDR_Eval_mbmb(SMDR_MZ_EXPT, 5);
    else SMDR_mbmb = 0;

    if (mcmc_target > SMDR_TOL) 
      SMDR_mcmc = SMDR_Eval_mcmc(SMDR_Mtau_EXPT, SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 
                                 5);
    else SMDR_mcmc = 0;
   
    SMDR_Eval_mquarks_2GeV (SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 5,
                            &SMDR_ms_2GeV, &SMDR_mu_2GeV, &SMDR_md_2GeV);

    SMDR_Mtau_pole = SMDR_Eval_Mtau_pole (SMDR_Mtau_EXPT,
                                          SMDR_mbmb_EXPT,
                                          SMDR_MZ_EXPT, 3);

    SMDR_Mmuon_pole = SMDR_Eval_Mmuon_pole (SMDR_mcmc_EXPT,
                                            SMDR_mcmc_EXPT,
                                            SMDR_Mtau_EXPT,
                                            SMDR_mbmb_EXPT,
                                            SMDR_MZ_EXPT, 2);
    SMDR_Melectron_pole = SMDR_Eval_Melectron_pole (SMDR_mcmc_EXPT,
                                                    SMDR_mcmc_EXPT,
                                                    SMDR_Mtau_EXPT,
                                                    SMDR_mbmb_EXPT,
                                                    SMDR_MZ_EXPT, 2);

    /*
    printf("\n");
    printf("SMDR_Mt_pole = %.14Lf;         frac. err. = %.10Lf\n", SMDR_Mt_pole, err_Mt);
    printf("SMDR_Mh_pole = %.14Lf;         frac. err. = %.10Lf\n", SMDR_Mh_pole, err_Mh);
    printf("SMDR_MZ_BreitWigner = %.14Lf;   frac. err. = %.10Lf\n", SMDR_MZ_BreitWigner, err_MZ);
    printf("SMDR_GFermi = %.20Lf;      frac. err. = %.10Lf\n", SMDR_GFermi, err_GFermi);
    printf("SMDR_alpha = 1/%.16Lf;         frac. err. = %.10Lf\n",1/SMDR_alpha, err_alpha);
    printf("SMDR_alphaS_5_MZ = %.16Lf;       frac. err. = %.10Lf\n",SMDR_alphaS_5_MZ, err_alphaS);
    printf("SMDR_mbmb = %.16Lf;              frac. err. = %.10Lf\n",SMDR_mbmb, err_mbmb);
    printf("SMDR_mcmc = %.16Lf;              frac. err. = %.10Lf\n",SMDR_mcmc, err_mcmc);
    printf("SMDR_ms_2GeV = %.16Lf;           frac. err. = %.10Lf\n",SMDR_ms_2GeV, err_ms_2GeV);
    printf("SMDR_md_2GeV = %.16Lf;           frac. err. = %.10Lf\n",SMDR_md_2GeV, err_md_2GeV);
    printf("SMDR_mu_2GeV = %.16Lf;           frac. err. = %.10Lf\n",SMDR_mu_2GeV, err_mu_2GeV);
    printf("SMDR_Mtau_pole = %.16Lf;         frac. err. = %.10Lf\n",SMDR_Mtau_pole, err_Mtau);
    printf("SMDR_Mmuon_pole = %.17Lf;        frac. err. = %.10Lf\n",SMDR_Mmuon_pole, err_Mmuon);
    printf("SMDR_Melectron_pole = %.19Lf;  frac. err. = %.10Lf\n",SMDR_Melectron_pole, err_Melectron);
    */

    result_niters = j;
  }

  SMDR_Load_Inputs();

  SMDR_Lambda = 0;
  SMDR_Lambda_in = SMDR_Lambda = -SMDR_Eval_Veffmin (-1, 3.5);
  SMDR_m2_in = SMDR_m2;

  SMDR_RGeval_SM (Q_target, 5);
  SMDR_Save_Inputs();

  return(result_niters);
};

/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
/* This function finds the coefficients used by SMDR_Fit_Inputs() to
   hopefully accelerate convergence in the iteration.  The results are
   put into the file accelcoeffs.h.  This function should be run once
   when the SMDR is first installed, and again if/when
   ReferenceModel.dat changes.
*/

void SMDR_Make_Accelcoeffs () 
{
  SMDR_REAL cog3_g3, cog3_g32;
  SMDR_REAL colambda_g3, coyt_g3, coyb_g3, coyc_g3, coys_g3, coyd_g3, coyu_g3;
  SMDR_REAL coyt_yt, coyt_yt2, colambda_yt, cog3_yt, coyb_yt, coyc_yt;
  SMDR_REAL coys_yt, coyd_yt, coyu_yt;
  SMDR_REAL coyb_yb, coyb_yb2, coyc_yb;
  SMDR_REAL coys_yb, coyd_yb, coyu_yb;
  SMDR_REAL coyc_yc, coyc_yc2;
  SMDR_REAL r1_g3, r2_g3, g3_p0;
  SMDR_REAL r1_lambda, lambda_p0;
  SMDR_REAL r1_yt, r2_yt, yt_p0;
  SMDR_REAL r1_yb, r2_yb, yb_p0;
  SMDR_REAL r1_yc, r2_yc, yc_p0;
  SMDR_REAL r1_ys, ys_p0;
  SMDR_REAL r1_yd, yd_p0;
  SMDR_REAL r1_yu, yu_p0;
  SMDR_REAL delt = 0.005;
  FILE *outfile;
  char outFileName[] = "accelcoeffs.h";

  /* SMDR_Read_Model_File ("ReferenceModel.dat"); */
#include "AUTOrefmodel.h"

  g3_p0 = SMDR_g3_in;
  lambda_p0 = SMDR_lambda_in;
  yt_p0 = SMDR_yt_in;
  yb_p0 = SMDR_yb_in;
  yc_p0 = SMDR_yc_in;
  ys_p0 = SMDR_ys_in;
  yd_p0 = SMDR_yd_in;
  yu_p0 = SMDR_yu_in;

  /* First consider variations in alphaS_5_MZ */
  SMDR_Fit_Inputs (SMDR_Q_in,
                   (1. + delt) * (1. + delt) * SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   SMDR_Mt_pole,
                   SMDR_mbmb,
                   SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   1.0e-15);

  r1_g3 = SMDR_g3_in/g3_p0;
  r1_lambda = SMDR_lambda_in/lambda_p0;
  r1_yt = SMDR_yt_in/yt_p0;
  r1_yb = SMDR_yb_in/yb_p0;
  r1_yc = SMDR_yc_in/yc_p0;
  r1_ys = SMDR_ys_in/ys_p0;
  r1_yd = SMDR_yd_in/yd_p0;
  r1_yu = SMDR_yu_in/yu_p0;


  /* SMDR_Read_Model_File ("ReferenceModel.dat"); */
#include "AUTOrefmodel.h"

  SMDR_Fit_Inputs (SMDR_Q_in,
                   (1. - delt) * (1. - delt) * SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   SMDR_Mt_pole,
                   SMDR_mbmb,
                   SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   1.0e-15);

  r2_g3 = SMDR_g3_in/g3_p0;

  cog3_g3 = (r1_g3 - r2_g3)/(2. * delt);
  cog3_g32 = (r1_g3 + r2_g3 - 2.)/(2. * delt * delt);

  colambda_g3 = (r1_lambda-1.)/delt;
  coyt_g3 = (r1_yt-1.)/delt;
  coyb_g3 = (r1_yb-1.)/delt;
  coyc_g3 = (r1_yc-1.)/delt;
  coys_g3 = (r1_ys-1.)/delt;
  coyd_g3 = (r1_yd-1.)/delt;
  coyu_g3 = (r1_yu-1.)/delt;

  /* Now consider variations in Mt_pole. */

  /* SMDR_Read_Model_File ("ReferenceModel.dat"); */
#include "AUTOrefmodel.h"

  SMDR_Fit_Inputs (SMDR_Q_in,
                   SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   (1. + delt) * SMDR_Mt_pole,
                   SMDR_mbmb,
                   SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   1.0e-15);

  r1_g3 = SMDR_g3_in/g3_p0;
  r1_lambda = SMDR_lambda_in/lambda_p0;
  r1_yt = SMDR_yt_in/yt_p0;
  r1_yb = SMDR_yb_in/yb_p0;
  r1_yc = SMDR_yc_in/yc_p0;
  r1_ys = SMDR_ys_in/ys_p0;
  r1_yd = SMDR_yd_in/yd_p0;
  r1_yu = SMDR_yu_in/yu_p0;


  /* SMDR_Read_Model_File ("ReferenceModel.dat"); */
#include "AUTOrefmodel.h"

  SMDR_Fit_Inputs (SMDR_Q_in,
                   SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   (1. - delt) * SMDR_Mt_pole,
                   SMDR_mbmb,
                   SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   1.0e-15);

  r2_yt = SMDR_yt_in/yt_p0;

  coyt_yt = (r1_yt - r2_yt)/(2. * delt);
  coyt_yt2 = (r1_yt + r2_yt - 2.)/(2. * delt * delt);

  colambda_yt = (r1_lambda-1.)/delt;
  cog3_yt = (r1_g3-1.)/delt;
  coyb_yt = (r1_yb-1.)/delt;
  coyc_yt = (r1_yc-1.)/delt;
  coys_yt = (r1_ys-1.)/delt;
  coyd_yt = (r1_yd-1.)/delt;
  coyu_yt = (r1_yu-1.)/delt;

  /* Now consider variations in mbmb. */

  /* SMDR_Read_Model_File ("ReferenceModel.dat"); */
#include "AUTOrefmodel.h"

  SMDR_Fit_Inputs (SMDR_Q_in,
                   SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   SMDR_Mt_pole,
                   (1. + delt) * SMDR_mbmb,
                   SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   1.0e-15);

  r1_yb = SMDR_yb_in/yb_p0;
  r1_yc = SMDR_yc_in/yc_p0;
  r1_ys = SMDR_ys_in/ys_p0;
  r1_yd = SMDR_yd_in/yd_p0;
  r1_yu = SMDR_yu_in/yu_p0;

  /* SMDR_Read_Model_File ("ReferenceModel.dat"); */
#include "AUTOrefmodel.h"

  SMDR_Fit_Inputs (SMDR_Q_in,
                   SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   SMDR_Mt_pole,
                   (1. - delt) * SMDR_mbmb,
                   SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   1.0e-15);

  r2_yb = SMDR_yb_in/yb_p0;

  coyb_yb = (r1_yb - r2_yb)/(2. * delt);
  coyb_yb2 = (r1_yb + r2_yb - 2.)/(2. * delt * delt);
  coyc_yb = (r1_yc - 1.)/delt;
  coys_yb = (r1_ys - 1.)/delt;
  coyd_yb = (r1_yd - 1.)/delt;
  coyu_yb = (r1_yu - 1.)/delt;

  /* Now consider variations in mcmc. */

  /* SMDR_Read_Model_File ("ReferenceModel.dat"); */
#include "AUTOrefmodel.h"

  SMDR_Fit_Inputs (SMDR_Q_in,
                   SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   SMDR_Mt_pole,
                   SMDR_mbmb,
                   (1. + delt) * SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   1.0e-15);

  r1_yc = SMDR_yc_in/yc_p0;

  /* SMDR_Read_Model_File ("ReferenceModel.dat"); */
#include "AUTOrefmodel.h"

  SMDR_Fit_Inputs (SMDR_Q_in,
                   SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   SMDR_Mt_pole,
                   SMDR_mbmb,
                   (1. - delt) * SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   1.0e-15);

  r2_yc = SMDR_yc_in/yc_p0;

  coyc_yc = (r1_yc - r2_yc)/(2. * delt);
  coyc_yc2 = (r1_yc + r2_yc - 2.)/(2. * delt * delt);

  /* Print these to the file accelcoeffs.h */
  outfile = fopen (outFileName, "w");

  fprintf(outfile,"SMDR_REAL cog3_g3     = %.10Lf;\n",cog3_g3);
  fprintf(outfile,"SMDR_REAL cog3_g32    = %.10Lf;\n",cog3_g32);
  fprintf(outfile,"SMDR_REAL colambda_g3 = %.10Lf;\n",colambda_g3);
  fprintf(outfile,"SMDR_REAL coyt_g3     = %.10Lf;\n",coyt_g3);
  fprintf(outfile,"SMDR_REAL coyb_g3     = %.10Lf;\n",coyb_g3);
  fprintf(outfile,"SMDR_REAL coyc_g3     = %.10Lf;\n",coyc_g3);
  fprintf(outfile,"SMDR_REAL coys_g3     = %.10Lf;\n",coys_g3);
  fprintf(outfile,"SMDR_REAL coyd_g3     = %.10Lf;\n",coyd_g3);
  fprintf(outfile,"SMDR_REAL coyu_g3     = %.10Lf;\n",coyu_g3);

  fprintf(outfile,"SMDR_REAL coyt_yt     = %.10Lf;\n",coyt_yt);
  fprintf(outfile,"SMDR_REAL coyt_yt2    = %.10Lf;\n",coyt_yt2);
  fprintf(outfile,"SMDR_REAL colambda_yt = %.10Lf;\n",colambda_yt);
  fprintf(outfile,"SMDR_REAL cog3_yt     = %.10Lf;\n",cog3_yt);
  fprintf(outfile,"SMDR_REAL coyb_yt     = %.10Lf;\n",coyb_yt);
  fprintf(outfile,"SMDR_REAL coyc_yt     = %.10Lf;\n",coyc_yt);
  fprintf(outfile,"SMDR_REAL coys_yt     = %.10Lf;\n",coys_yt);
  fprintf(outfile,"SMDR_REAL coyd_yt     = %.10Lf;\n",coyd_yt);
  fprintf(outfile,"SMDR_REAL coyu_yt     = %.10Lf;\n",coyu_yt);

  fprintf(outfile,"SMDR_REAL coyb_yb     = %.10Lf;\n",coyb_yb);
  fprintf(outfile,"SMDR_REAL coyb_yb2    = %.10Lf;\n",coyb_yb2);
  fprintf(outfile,"SMDR_REAL coyc_yb     = %.10Lf;\n",coyc_yb);
  fprintf(outfile,"SMDR_REAL coys_yb     = %.10Lf;\n",coys_yb);
  fprintf(outfile,"SMDR_REAL coyd_yb     = %.10Lf;\n",coyd_yb);
  fprintf(outfile,"SMDR_REAL coyu_yb     = %.10Lf;\n",coyu_yb);

  fprintf(outfile,"SMDR_REAL coyc_yc     = %.10Lf;\n",coyc_yc);
  fprintf(outfile,"SMDR_REAL coyc_yc2    = %.10Lf;\n",coyc_yc2);

  fclose (outfile);
  return;
}
