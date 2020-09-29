#include "smdr_internal.h"

/* 
  Sources: 
  9912391 and 0608026 for 3-loop QCD+QED effects (b, tau, mu, e)  
  1502.01030 and 1606.06754 for 4-loop QCD effects (b only). 
  0708.1729 for lighter fermion mass effects.
  See Mathematica files QCDQED and makeMbot4loop for derivations.
*/

/* --------------------------------------------------------------------- */
SMDR_REAL SMDR_Eval_mbmb (SMDR_REAL Q_dec_thZW, int loopOrder)
{
  SMDR_REAL mb_mb, mb_mb_previous;
  SMDR_REAL alphaS_mbmb, alpha_mbmb, cu, cd, ce;
  int max_iters = 20;
  int i;

  /* Use PDG experimental value as our initial guess. */
  SMDR_RGeval_QCDQED_53 (SMDR_mbmb_EXPT, Q_dec_thZW, loopOrder);
  mb_mb = mb_53;
  if (mb_mb < SMDR_TOL) return (0);

  for (i=0; i<max_iters; i++) {
    /*
    printf("  Iteration %d:  mb(mb) = %.10Lf;\n", i, mb_mb);
    */
    mb_mb_previous = mb_mb;
    SMDR_RGrun_QCDQED (Q_53, mb_mb, loopOrder, 2, 3, 3,
                       alphaS_53, alpha_53,
                       &alphaS_mbmb, &alpha_mbmb, &cu, &cd, &ce);

    /* This seems to minimize iterations needed for convergence: */
    mb_mb = (mb_mb + 5.*mb_53 * cd)/6.;
    if ( SMDR_FABS(1.0 - mb_mb/mb_mb_previous) < 1.0e-14) break;
  }

  Q_53 = mb_53 = mb_mb;
  alphaS_53 = alphaS_mbmb;
  alpha_53 = alpha_mbmb;
  mc_53 *= cu;
  ms_53 *= cd;
  md_53 *= cd;
  mu_53 *= cu;
  mtau_53 *= ce;
  mmuon_53 *= ce;
  melectron_53 *= ce;

  return (mb_53);
}

/* --------------------------------------------------------------------- */
SMDR_REAL SMDR_Eval_mcmc (SMDR_REAL Q_dec_tau,
                          SMDR_REAL Q_dec_bottom,
                          SMDR_REAL Q_dec_thZW, int loopOrder)
{
  SMDR_REAL mc_mc, mc_mc_previous;
  SMDR_REAL alphaS_mcmc, alpha_mcmc, cu, cd, ce;
  int max_iters = 20;
  int i;

  /* Use PDG experimental value as our initial guess. */
  SMDR_RGeval_QCDQED_42 (SMDR_mcmc_EXPT,
                         Q_dec_tau, 
                         Q_dec_bottom, 
                         Q_dec_thZW, loopOrder);
  mc_mc = mc_42;
  if (mc_mc < SMDR_TOL) return (0);

  for (i=0; i<max_iters; i++) {
    /*
    printf("  Iteration %d:  mc(mc) = %.10Lf;\n", i, mc_mc);
    */
    mc_mc_previous = mc_mc;
    SMDR_RGrun_QCDQED (Q_42, mc_mc, loopOrder, 2, 2, 2,
                       alphaS_42, alpha_42,
                       &alphaS_mcmc, &alpha_mcmc, &cu, &cd, &ce);
    /* This seems to minimize iterations needed for convergence: */
    mc_mc = (mc_mc + 2.*mc_42 * cu)/3.;
    if ( SMDR_FABS(1.0 - mc_mc/mc_mc_previous) < 1.0e-14) break;
  }

  Q_42 = mc_42 = mc_mc;
  alphaS_42 = alphaS_mcmc;
  alpha_42 = alpha_mcmc;
  ms_42 *= cd;
  md_42 *= cd;
  mu_42 *= cu;
  mmuon_42 *= ce;
  melectron_42 *= ce;

  return (mc_mc);
}

/* --------------------------------------------------------------------- */
/* 
  Returns the pole mass of the bottom quark, computed in the 5-quark,
  3-charged-lepton QCD+QED effective field theory at the scale Q_eval,
  in the approximation of loopOrder loops.
  The starting inputs are the global variables corresponding to the MSbar
  parameters of the full Standard Model:
  Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.
  Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW.     

  NOTE: This just doesn't converge fast enough to be useful! 
  The Review of Particle Properties uses the 2-loop approximation, which 
  is horrible. See 1606.06754 for discussion of the atrocious convergence 
  properties of the expansion. Therefore, the bottom (and charm) pole masses 
  are kind of useless as observables and are deprecated.
*/
SMDR_REAL SMDR_Eval_Mb_pole (SMDR_REAL Q_eval, 
                             SMDR_REAL Q_dec_thZW, 
                             int loopOrder)
{
  SMDR_REAL lnbarB, lnbarB2, lnbarB3, lnbarB4;
  SMDR_REAL zmb;
  SMDR_REAL Mbpole = 5;
  SMDR_REAL aSoPi, aoPi;
  SMDR_REAL aSoPi2, aSoPi3, aSoPi4, aoPi2, aoPi3;
  SMDR_REAL xf, Lxf;
  SMDR_REAL zmb2charm, zmb3charm;
  SMDR_REAL oldresult;  
  int i;
  char funcname[] = "SMDR_Eval_Mb_pole";

  if ((loopOrder < 0) || (loopOrder > 4)) SMDR_Error (funcname,
    "Invalid loopOrder, should be 0, 1, 2, 3, or 4.", 1);

  SMDR_RGeval_QCDQED_53 (Q_eval, Q_dec_thZW, 5);

  aSoPi = alphaS_53/PI;
  aoPi = alpha_53/PI;
  aSoPi2 = aSoPi * aSoPi;
  aSoPi3 = aSoPi2 * aSoPi;
  aSoPi4 = aSoPi3 * aSoPi;
  aoPi2 = aoPi * aoPi;
  aoPi3 = aoPi2 * aoPi;

  for (i = 0; i < 20; i++) {
    oldresult = Mbpole;
    lnbarB = 2. * SMDR_LOG (Mbpole/Q_53);
    lnbarB2 = lnbarB * lnbarB;
    lnbarB3 = lnbarB2 * lnbarB;
    lnbarB4 = lnbarB3 * lnbarB;

    zmb = 1;

    if (loopOrder > 0) 
      zmb += (3.*lnbarB - 4.) * (aSoPi/3. + aoPi/36.);

    /* Charm quark mass effects from 0708.1729 */
    if ((mc_53 > 0.000001) && (loopOrder > 1)) {
      xf = mc_53/Mbpole;
      Lxf = SMDR_LOG(xf);

      /* From eqs. (20), (21) of 0708.1729 */
      /*
      zmb2charm = xf * (-1.6449340668482262 * (1 + xf*xf) + xf + 
                   xf*xf*xf*(1.24739 - 0.722222 * Lxf + 0.333333 * Lxf * Lxf)); 
      */

      /* Alternate version, checked that this is the same. */
      zmb2charm = -SMDR_f2lf(xf)/48.;

      /* From eqs. (20), (21) of 0708.1729 */
      zmb3charm = xf * (-17.9766 + 13.7079 * Lxf + 
                    xf * (12.0355 + 10.687 * Lxf) +
                    xf*xf* (12.286 + 13.2205 * Lxf) +
                    xf*xf*xf * (-3.0196 - 0.6528*Lxf + 0.0698*Lxf*Lxf - 
                      0.0370*Lxf*Lxf*Lxf));

      /* From eqs. (21), (22) of 0708.1729 */
      zmb3charm += (2. + 1.5 * SMDR_LOG(Q_53*Q_53/(mc_53 * mc_53))) * xf * (
         -1.096623 + 1.333333*xf - 3.289868*xf*xf + 
         xf*xf*xf*(2.8448797 - 1.481481*Lxf + 0.8888889*Lxf*Lxf) +
         xf*xf*xf*xf*xf*(-0.2785185 + 0.3555556*Lxf));
    } else {
      zmb2charm = 0;
      zmb3charm = 0;
    }

    if (loopOrder > 1) {
      zmb += aSoPi2*(-10.16681750113412 + 4.736111111111111*lnbarB - 
               0.4583333333333333*lnbarB2) + 
             aoPi*aSoPi*(-0.1512777246760184 - 0.19444444444444444*lnbarB + 
               0.08333333333333333*lnbarB2) + 
             aoPi2*(1.1295686244986776 - 0.40933641975308643*lnbarB + 
               0.09606481481481481*lnbarB2);

      zmb += aSoPi2*zmb2charm;
    }

    if (loopOrder > 2) {
      zmb += aSoPi3*(-101.45441694139859 + 34.762883169060714*lnbarB - 
               6.16087962962963*lnbarB2 + 0.43287037037037035*lnbarB3) +
             aoPi*aSoPi2*(0.9294775699448863 - 0.8105858687538399*lnbarB + 
               0.5231481481481481*lnbarB2 - 0.03819444444444445*lnbarB3) + 
             aoPi2*aSoPi*(0.6438042387488245 + 0.9882375058281452*lnbarB - 
               0.4533179012345679*lnbarB2 + 0.09606481481481481*lnbarB3) + 
             aoPi3*(-8.999958085162628 + 4.706835233020573*lnbarB - 
               0.888320901920439*lnbarB2 + 0.14498671124828533*lnbarB3); 

      zmb += aSoPi3*zmb3charm;
    }

    if (loopOrder > 3)
      /* From 1606.06754 equation (20) and (43), used ancillary file. */
      zmb += aSoPi4 * (-1278.7 + 500.23310595472594*lnbarB - 
             83.38977278773129*lnbarB2 + 9.95630787037037*lnbarB3 - 
             0.5140335648148148*lnbarB4);

    Mbpole = mb_53/zmb;
    if (SMDR_FABS(oldresult - Mbpole) < 1.0e-14) break;
  }

  return (Mbpole);
}

/* ----------------------------------------------------------------------- */
SMDR_REAL twoloopQEDpart(SMDR_REAL nu, 
                         SMDR_REAL nd, 
                         SMDR_REAL ne,
                         SMDR_REAL lnbar) 
{
  SMDR_REAL lnbar2 = lnbar * lnbar;

  return( 103./128. - (21.*lnbar)/32. + (9.*lnbar2)/32.  
          - (27.*Zeta2)/8. +  3.*ln2*Zeta2 - (3.*Zeta3)/4. 
          + nd*(71./288. - (13.*lnbar)/72. + lnbar2/24. + Zeta2/6.)   
          + ne*(71./96. - (13.*lnbar)/24. + lnbar2/8. + Zeta2/2.) 
          + nu*(71./72. - (13.*lnbar)/18. + lnbar2/6. + (2.*Zeta2)/3.));
}

/* ----------------------------------------------------------------------- */
/*
   Tau lepton pole mass, working in the effective QCD+QED theory with
   u,d,s,c active quarks and tau,mu,e active leptons.
*/
SMDR_REAL SMDR_Eval_Mtau_pole (SMDR_REAL Q_eval, 
                               SMDR_REAL Q_dec_bottom,
                               SMDR_REAL Q_dec_thZW,
                               int loopOrder)
{
  SMDR_REAL Mtaupole = SMDR_Mtau_EXPT; /* First guess is experimental value. */
  SMDR_REAL aSoPi, aoPi, aoPi2;
  SMDR_REAL lnbarTau, lnbarTau2, lnbarTau3;
  SMDR_REAL zmtau;
  SMDR_REAL oldresult;  
  int max_iters = 20;
  int i;
  char funcname[] = "SMDR_Mtaupole";

  if ((loopOrder < 0) || (loopOrder > 3)) SMDR_Error (funcname,
    "Invalid loopOrder, should be 0, 1, 2, or 3.", 1);

  SMDR_RGeval_QCDQED_43 (Q_eval, Q_dec_bottom, Q_dec_thZW, 5);

  if (mtau_43 < SMDR_TOL) return(0);

  aSoPi = alphaS_43/PI;
  aoPi = alpha_43/PI;
  aoPi2 = aoPi * aoPi;

  for (i = 0; i < max_iters; i++) {
    oldresult = Mtaupole;
    lnbarTau = 2. * SMDR_LOG (Mtaupole/Q_43);
    lnbarTau2 = lnbarTau * lnbarTau;
    lnbarTau3 = lnbarTau2 * lnbarTau;
    zmtau = 1;

    if (loopOrder > 0) 
      zmtau += aoPi * ((3./4.)*lnbarTau - 1.);

    if (loopOrder > 1) { 
      zmtau += aoPi2 * twoloopQEDpart(2,2,3,lnbarTau);

      /* Charm and muon mass effect corrections are very small. */
      zmtau += -(aoPi2/24.) * SMDR_f2lf(SMDR_mcmc_EXPT/SMDR_Mtau_EXPT);
      zmtau += -(aoPi2/32.) * SMDR_f2lf(SMDR_Mmuon_EXPT/SMDR_Mtau_EXPT); 

      /* Strange-quark mass correction uses mquark = 400 MeV, very small. */
      zmtau += -(aoPi2/96.) * SMDR_f2lf(0.4/SMDR_Mtau_EXPT); 

      /* Down and up mass corrections both use mquark = 300 MeV, very small. */
      zmtau += -(5.*aoPi2/96.) * SMDR_f2lf(0.3/SMDR_Mtau_EXPT); 
    }

    if (loopOrder > 2) 
      zmtau += aoPi2 * (
        aSoPi*(3.1635406840679976 - 0.22925476724579677*lnbarTau + 
          0.41666666666666667*lnbarTau2) + 
        aoPi*(-48.67926994511686 + 36.77729721855527*lnbarTau - 
          10.220582561728397*lnbarTau2 + 1.7782600308641976*lnbarTau3));

    Mtaupole = mtau_43/zmtau;

    if (SMDR_FABS(oldresult - Mtaupole) < 1.0e-14) break;
  }

  return (Mtaupole);
}

/* ----------------------------------------------------------------------- */
/*
   Lepton pole mass, working in the effective QCD+QED theory with
   s, d, u active quarks and e, mu active leptons.
*/
SMDR_REAL eval_Mlepton_pole (SMDR_REAL mleptonMS,
                             SMDR_REAL alpha_,
                             SMDR_REAL Q_,
                             SMDR_REAL m_other_lepton,
                             int loopOrder)
{
  SMDR_REAL lnbar;
  SMDR_REAL zm;
  SMDR_REAL Mpole = mleptonMS;
  SMDR_REAL aoPi = alpha_/PI;
  SMDR_REAL aoPi2 = aoPi * aoPi;
  SMDR_REAL oldresult;  
  int max_iters = 20;
  int i;

  if (mleptonMS < SMDR_TOL) return (0);

  for (i = 0; i < max_iters; i++) {
    oldresult = Mpole;
    lnbar = 2. * SMDR_LOG (Mpole/Q_);

    zm = 1;

    if (loopOrder > 0) 
      zm += aoPi * ((3./4.)*lnbar - 1.);

    /* To avoid numerical disasters, only include the 2-loop parts if
       both the electron and muon masses are semi-realistic. */
    if ((loopOrder > 1) && (mleptonMS > 0.0001) && (m_other_lepton > 0.0001))
    { 
      zm += aoPi2 * twoloopQEDpart(1,2,2,lnbar);

      /* Mass correction from other active lepton */
      zm += -(aoPi2/32.) * SMDR_f2lf(m_other_lepton/Mpole); 

      /* Strange-quark mass correction uses mquark = 400 MeV. */
      zm += -(aoPi2/96.) * SMDR_f2lf(0.4/Mpole); 

      /* Down and up mass corrections both use mquark = 300 MeV. */
      zm += -(5.*aoPi2/96.) * SMDR_f2lf(0.3/Mpole); 
    }

    Mpole = mleptonMS/zm;

    if (SMDR_FABS(oldresult - Mpole) < 1.0e-14) break;
  }

  return (Mpole);
}

/* ----------------------------------------------------------------------- */
/*
   Muon pole mass, working in the effective QCD+QED theory with
   s, d, u active quarks and e, mu active leptons.
*/
SMDR_REAL SMDR_Eval_Mmuon_pole (SMDR_REAL Q_eval, 
                                SMDR_REAL Q_dec_charm,
                                SMDR_REAL Q_dec_tau,
                                SMDR_REAL Q_dec_bottom,
                                SMDR_REAL Q_dec_thZW,
                                int loopOrder)
{
  char funcname[] = "SMDR_Mmuon_pole";

  if ((loopOrder < 0) || (loopOrder > 2)) SMDR_Error (funcname,
    "Invalid loopOrder, should be 0, 1, or 2.", 1);

  SMDR_RGeval_QCDQED_32 (Q_eval, 
                         Q_dec_charm, 
                         Q_dec_tau,
                         Q_dec_bottom, 
                         Q_dec_thZW, 
                         5);

  return (eval_Mlepton_pole(mmuon_32, alpha_32, Q_32, SMDR_Melectron_EXPT, 
                            loopOrder));
}

/* ----------------------------------------------------------------------- */
/*
   Electron pole mass, working in the effective QCD+QED theory with
   s, d, u active quarks and e, mu active leptons.
*/
SMDR_REAL SMDR_Eval_Melectron_pole (SMDR_REAL Q_eval, 
                                    SMDR_REAL Q_dec_charm,
                                    SMDR_REAL Q_dec_tau,
                                    SMDR_REAL Q_dec_bottom,
                                    SMDR_REAL Q_dec_thZW,
                                    int loopOrder)
{
  char funcname[] = "SMDR_Melectron_pole";

  if ((loopOrder < 0) || (loopOrder > 2)) SMDR_Error (funcname,
    "Invalid loopOrder, should be 0, 1, or 2.", 1);

  SMDR_RGeval_QCDQED_32 (Q_eval, 
                         Q_dec_charm, 
                         Q_dec_tau,
                         Q_dec_bottom, 
                         Q_dec_thZW, 
                         5);

  return (eval_Mlepton_pole(melectron_32, alpha_32, Q_32, SMDR_Mmuon_EXPT, 
                            loopOrder));
}

/* --------------------------------------------------------------------- */
/* 
   This function is used to include the 2-loop effects of non-zero fermion
   mass mf to the pole mass MF of a heavier fermion. x = mf/MF. The 
   function vanishes when x=0. 

   For small x, this asymptotes to 8 Pi^2 x - 48 x^2 +... .

   For x = 1, it equals 8 Pi^2 - 24.

   For large x it asymptotes to: 
     302/9 + (8*Pi^2)/3 + (104*Log[x])/3 + 16*Log[x]^2
*/
SMDR_REAL SMDR_f2lf (SMDR_REAL x) 
{
  SMDR_REAL logx;
  SMDR_REAL result;
  SMDR_REAL x2, x3, x4;

  if (x < TSIL_TOL) return(0);

  logx = SMDR_LOG(x);
  x2 = x*x;
  x3 = x2*x;
  x4 = x3*x;

  result = SMDR_CREAL( -24.*x2 + 8.*(1. + x2)*(1. + 3.*x - x2)*Zeta2 
    - 16.*x2*logx + 16.*logx*logx + 8.*(1. + x4)*TSIL_Dilog(1. - 1./x2)  
    + 16.*x*(1. + x2)*(TSIL_Dilog((1. - x)/(1. + x)) - 
    TSIL_Dilog((-1. + x)/(1. + x))) );

  return (result);
}

/* --------------------------------------------------------------------- */

void SMDR_Eval_Light_Masses () {

  SMDR_Eval_QCDQED_at_MZ (SMDR_MZ_EXPT, SMDR_MZ_EXPT, 5);

  SMDR_mbmb = SMDR_Eval_mbmb (SMDR_MZ_EXPT, 5);
  SMDR_mcmc = SMDR_Eval_mcmc (SMDR_Mtau_EXPT, SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 5);
  
  /* This is deprecated! Pole mass of bottom doesn't converge fast enough.
  SMDR_Mb_pole = SMDR_Eval_Mb_pole (SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 4);
  */

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

  return;
}

/* ----------------------------------------------------------- */
void SMDR_Eval_mquarks_2GeV (SMDR_REAL Q_dec_bottom,
                             SMDR_REAL Q_dec_thZW, 
                             int loopOrder,
                             SMDR_REAL *ms, 
                             SMDR_REAL *mu, 
                             SMDR_REAL *md)
{
  SMDR_RGeval_QCDQED_43 (2.0, Q_dec_bottom, Q_dec_thZW, loopOrder);
  *ms = SMDR_ms_43;
  *mu = SMDR_mu_43;
  *md = SMDR_md_43;
  return;
}

/* ----------------------------------------------------------- */

