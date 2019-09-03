/* RG running of the basic model paramters. */

#include "smdr_internal.h"

/* ---------------------------------------------------------------- */
int SMDR_RGeval_QCDQED_53 (SMDR_REAL Qfinal, 
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder) 
{
  int loopOrderDec = loopOrder;
  if (loopOrderDec > 4) loopOrderDec = 4;

  SMDR_RGeval_SM (Q_thZW_dec, loopOrder);
  SMDR_Decouple_thZW (loopOrderDec);

  return (SMDR_RGrun_QCDQED_53 (Qfinal, loopOrder));
}

/* ---------------------------------------------------------------- */
int SMDR_RGeval_QCDQED_43 (SMDR_REAL Qfinal, 
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder) 
{
  int loopOrderDec = loopOrder;
  if (loopOrderDec > 4) loopOrderDec = 4;

  SMDR_RGeval_QCDQED_53 (Q_b_dec, Q_thZW_dec, loopOrder);
  SMDR_Decouple_bottom (loopOrderDec);
  return (SMDR_RGrun_QCDQED_43 (Qfinal, loopOrder));
}

/* ---------------------------------------------------------------- */
int SMDR_RGeval_QCDQED_42 (SMDR_REAL Qfinal, 
                           SMDR_REAL Q_tau_dec,
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder) 
{
  int loopOrderDec = loopOrder;
  if (loopOrderDec > 4) loopOrderDec = 4;

  SMDR_RGeval_QCDQED_43 (Q_tau_dec, Q_b_dec, Q_thZW_dec, loopOrder);
  SMDR_Decouple_tau (loopOrderDec);
  return (SMDR_RGrun_QCDQED_42 (Qfinal, loopOrder));
}

/* ---------------------------------------------------------------- */
int SMDR_RGeval_QCDQED_32 (SMDR_REAL Qfinal, 
                           SMDR_REAL Q_c_dec,
                           SMDR_REAL Q_tau_dec,
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder) 
{
  int loopOrderDec = loopOrder;
  if (loopOrderDec > 4) loopOrderDec = 4;

  SMDR_RGeval_QCDQED_42 (Q_c_dec, Q_tau_dec, Q_b_dec, Q_thZW_dec, loopOrder);
  SMDR_Decouple_charm (loopOrderDec);
  return (SMDR_RGrun_QCDQED_32 (Qfinal, loopOrder));
}

/* ---------------------------------------------------------------- */
int SMDR_RGrun_QCDQED_53 (SMDR_REAL Qfinal, int loopOrder) 
{
  SMDR_REAL alphaS_final, alpha_final;
  SMDR_REAL cu, cd, ce;

  int success = SMDR_RGrun_QCDQED (Q_53, Qfinal, loopOrder, 2, 3, 3,
                                   alphaS_53, alpha_53, 
                                   &alphaS_final, &alpha_final,
                                   &cu, &cd, &ce);

  Q_53 = Qfinal;
  alphaS_53 = alphaS_final;
  alpha_53 = alpha_final;
  mb_53 *= cd;
  mc_53 *= cu;
  ms_53 *= cd;
  mu_53 *= cu; 
  md_53 *= cd;
  mtau_53 *= ce; 
  mmuon_53 *= ce; 
  melectron_53 *= ce; 

  return(success);
}

/* ---------------------------------------------------------------- */
int SMDR_RGrun_QCDQED_43 (SMDR_REAL Qfinal, int loopOrder) 
{
  SMDR_REAL alphaS_final, alpha_final;
  SMDR_REAL cu, cd, ce;
  int success = SMDR_RGrun_QCDQED (Q_43, Qfinal, loopOrder, 2, 2, 3,
                                   alphaS_43, alpha_43, 
                                   &alphaS_final, &alpha_final,
                                   &cu, &cd, &ce);
  Q_43 = Qfinal;
  alphaS_43 = alphaS_final;
  alpha_43 = alpha_final;
  mc_43 *= cu;
  ms_43 *= cd;
  mu_43 *= cu; 
  md_43 *= cd;
  mtau_43 *= ce; 
  mmuon_43 *= ce; 
  melectron_43 *= ce; 

  return(success);
}

/* ---------------------------------------------------------------- */
int SMDR_RGrun_QCDQED_42 (SMDR_REAL Qfinal, int loopOrder) 
{
  SMDR_REAL alphaS_final, alpha_final;
  SMDR_REAL cu, cd, ce;

  int success = SMDR_RGrun_QCDQED (Q_42, Qfinal, loopOrder, 2, 2, 2,
                                   alphaS_42, alpha_42, 
                                   &alphaS_final, &alpha_final,
                                   &cu, &cd, &ce);
  Q_42 = Qfinal;
  alphaS_42 = alphaS_final;
  alpha_42 = alpha_final;
  mc_42 *= cu;
  ms_42 *= cd;
  mu_42 *= cu; 
  md_42 *= cd;
  mmuon_42 *= ce; 
  melectron_42 *= ce; 

  return(success);
}

/* ---------------------------------------------------------------- */
int SMDR_RGrun_QCDQED_32 (SMDR_REAL Qfinal, int loopOrder) 
{
  SMDR_REAL alphaS_final, alpha_final;
  SMDR_REAL cu, cd, ce;
  int success = SMDR_RGrun_QCDQED (Q_32, Qfinal, loopOrder, 1, 2, 2,
                                   alphaS_32, alpha_32, 
                                   &alphaS_final, &alpha_final,
                                   &cu, &cd, &ce);
  Q_32 = Qfinal;
  alphaS_32 = alphaS_final;
  alpha_32 = alpha_final;
  ms_32 *= cd;
  mu_32 *= cu; 
  md_32 *= cd;
  mmuon_32 *= ce; 
  melectron_32 *= ce; 

  return(success);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Runs the parameters to Q_final at loop order specified.          */

int SMDR_RGrun_QCDQED (SMDR_REAL Q_init, SMDR_REAL Q_final, 
                       int loopOrder, int nu, int nd, int ne,
                       SMDR_REAL alphaS_init, SMDR_REAL alpha_init,
                       SMDR_REAL *alphaS_final, SMDR_REAL *alpha_final,
                       SMDR_REAL *cu, SMDR_REAL *cd, SMDR_REAL *ce)
{
  SMDR_REAL t_final, t, dt;
  int min_steps, max_steps;
  int force_step;
  int goodsteps, badsteps; 
  int rk6status; /* 1 for success or forced; 0 for need retry. */
  char funcname[] = "QCDQEDrun";

  if ( (loopOrder < 1) || (loopOrder > 5) )
    SMDR_Error (funcname,
    "Invalid loop order specified, should be 0, 1, 2, 3, 4, or 5", 3);
  
  Q_eff = Q_init;
  alphaS_run = alphaS_init;
  alpha_run = alpha_init;
  lnmu = 0;
  lnmd = 0;
  lnme = 0;

  t_final = TSIL_LOG(Q_final/Q_init);

  if (TSIL_FABS(t_final) < TSIL_TOL) {
    *alphaS_final = alphaS_init;
    *alpha_final = alpha_init;
    *cu = 1;
    *cd = 1;
    *ce = 1;    
    return 0;
  }

  min_steps = 1 + ((int) (5.0L * TSIL_FABS(t_final)));
  max_steps = 1 + ((int) (125.0L * TSIL_FABS(t_final)));

  dt = t_final/(1 + ((int) (124.0L * TSIL_FABS(t_final))));
  t = 0.0L;

  SMDR_Betas_QCDQED (loopOrder, nu, nd, ne);

  goodsteps = badsteps = 0;

  while ( TSIL_FABS(dt) < 0.5 * TSIL_FABS(t_final - t) )
    {
      if ( TSIL_FABS(dt) < TSIL_FABS(t_final/max_steps) )
	{
	  force_step = 1;
	  dt = t_final/((SMDR_REAL) max_steps);
	}
      else force_step = 0;
      
      if ( TSIL_FABS(dt) > TSIL_FABS(t_final)/min_steps )
	dt = (t_final)/min_steps;

      rk6status = SMDR_RG_rk6_QCDQED (&dt, force_step, loopOrder, nu, nd, ne); 
      t = TSIL_LOG (Q_eff/Q_init);
      
      if (1 == rk6status) {
	goodsteps += (1 - force_step);
	badsteps += force_step;
      }
    }

  /* The remaining distance is less than 2.0 times the step size.  So,
     take exactly two more steps. */
  dt = 0.5L*(t_final - t);
  SMDR_RG_rk6_QCDQED (&dt, 1, loopOrder, nu, nd, ne);
  t = TSIL_LOG (Q_eff/Q_init);

  /* Arrange final step to land exactly on t_final, and force it. */
  dt = t_final - t;
  SMDR_RG_rk6_QCDQED (&dt, 1, loopOrder, nu, nd, ne);

  /*
  printf("goodsteps = %d    badsteps  = %d\n",goodsteps, badsteps);
  */

  *alphaS_final = alphaS_run;
  *alpha_final = alpha_run;
  *cu = SMDR_EXP(lnmu);
  *cd = SMDR_EXP(lnmd);
  *ce = SMDR_EXP(lnme);

  /* Return a status code eventually */
  return 0;
}
