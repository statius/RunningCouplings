/* 
   Higgs self-energy functions and complex pole mass calculation. 
   The computation uses equations (2.43), (2.46)-(2.48), and (3.2)-(3.4) of
   arXiv:1407.4336. The argument loopOrder has allowed values:

   0    tree level
   1    1-loop level
   1.5  1-loop plus 2-loop QCD
   2    2-loop level
   2.3  2-loop plus leading 3-loop QCD
   2.5  2-loop plus leading 3-loop QCD and yt
*/

#include "smdr_internal.h"

/* Parameters conrolling iteration: */
#define MAXITERS 20

/* Local functions; not needed or useful elsewhere: */
SMDR_COMPLEX SMDR_Mh_Pi1h (SMDR_REAL);
SMDR_COMPLEX SMDR_Mh_Pi2QCDh (SMDR_REAL);
SMDR_COMPLEX SMDR_Mh_Pi2nonQCDh (SMDR_REAL);
SMDR_COMPLEX SMDR_Mh_Pi3QCDh (void);
SMDR_COMPLEX SMDR_Mh_Pi3nonQCDh (void);
int          SMDR_Mh_DoTSILh (float, SMDR_REAL);

/* #define global variables used in this file, for convenience with safety. */
#define AW SMDR_Mh_AW
#define AZ SMDR_Mh_AZ
#define Ah SMDR_Mh_Ah
#define At SMDR_Mh_At
#define Ab SMDR_Mh_Ab
#define BWW SMDR_Mh_BWW
#define BZZ SMDR_Mh_BZZ
#define Bhh SMDR_Mh_Bhh
#define Btt SMDR_Mh_Btt
#define Bbb SMDR_Mh_Bbb
#define Btautau SMDR_Mh_Btautau
#define Ihhh SMDR_Mh_Ihhh
#define IttZ SMDR_Mh_IttZ
#define Ihtt SMDR_Mh_Ihtt
#define I0tW SMDR_Mh_I0tW
#define IWWZ SMDR_Mh_IWWZ
#define IhWW SMDR_Mh_IhWW
#define IhZZ SMDR_Mh_IhZZ
#define I00W SMDR_Mh_I00W
#define I00Z SMDR_Mh_I00Z
#define I00h SMDR_Mh_I00h
#define I00t SMDR_Mh_I00t
#define I0hW SMDR_Mh_I0hW
#define I0hZ SMDR_Mh_I0hZ
#define I0WZ SMDR_Mh_I0WZ
#define UWW00 SMDR_Mh_UWW00
#define UZZ00 SMDR_Mh_UZZ00
#define S00W SMDR_Mh_S00W
#define S00Z SMDR_Mh_S00Z
#define TW00 SMDR_Mh_TW00
#define TZ00 SMDR_Mh_TZ00
#define MWWWW0 SMDR_Mh_MWWWW0
#define TW0W SMDR_Mh_TW0W
#define Tbar0WW SMDR_Mh_Tbar0WW
#define Mtttt0 SMDR_Mh_Mtttt0
#define Tt0t SMDR_Mh_Tt0t
#define Tbar0tt SMDR_Mh_Tbar0tt
#define Mbbbb0 SMDR_Mh_Mbbbb0
#define Tb0b SMDR_Mh_Tb0b
#define Tbar0bb SMDR_Mh_Tbar0bb
#define Mhhhhh SMDR_Mh_Mhhhhh
#define Uhhhh SMDR_Mh_Uhhhh
#define Shhh SMDR_Mh_Shhh
#define MttttZ SMDR_Mh_MttttZ
#define UtttZ SMDR_Mh_UtttZ
#define SttZ SMDR_Mh_SttZ
#define TttZ SMDR_Mh_TttZ
#define TZtt SMDR_Mh_TZtt
#define Mtttth SMDR_Mh_Mtttth
#define Uttth SMDR_Mh_Uttth
#define Shtt SMDR_Mh_Shtt
#define Ttth SMDR_Mh_Ttth
#define Thtt SMDR_Mh_Thtt
#define MtZtZt SMDR_Mh_MtZtZt
#define UZZtt SMDR_Mh_UZZtt
#define Mththt SMDR_Mh_Mththt
#define Uhhtt SMDR_Mh_Uhhtt
#define MtWtW0 SMDR_Mh_MtWtW0
#define UWW0t SMDR_Mh_UWW0t
#define Utt0W SMDR_Mh_Utt0W
#define S0tW SMDR_Mh_S0tW
#define TW0t SMDR_Mh_TW0t
#define Tt0W SMDR_Mh_Tt0W
#define Tbar0tW SMDR_Mh_Tbar0tW
#define MWWWWZ SMDR_Mh_MWWWWZ
#define UWWWZ SMDR_Mh_UWWWZ
#define SWWZ SMDR_Mh_SWWZ
#define TWWZ SMDR_Mh_TWWZ
#define TZWW SMDR_Mh_TZWW
#define MWWWWh SMDR_Mh_MWWWWh
#define UWWWh SMDR_Mh_UWWWh
#define ShWW SMDR_Mh_ShWW
#define TWWh SMDR_Mh_TWWh
#define ThWW SMDR_Mh_ThWW
#define MZZZZh SMDR_Mh_MZZZZh
#define UZZZh SMDR_Mh_UZZZh
#define ShZZ SMDR_Mh_ShZZ
#define TZZh SMDR_Mh_TZZh
#define ThZZ SMDR_Mh_ThZZ
#define MWZWZW SMDR_Mh_MWZWZW
#define UZZWW SMDR_Mh_UZZWW
#define MhWhWW SMDR_Mh_MhWhWW
#define UhhWW SMDR_Mh_UhhWW
#define MhZhZZ SMDR_Mh_MhZhZZ
#define UhhZZ SMDR_Mh_UhhZZ
#define lnbart SMDR_Mh_lnbart
#define lnbart2 SMDR_Mh_lnbart2
#define lnbart3 SMDR_Mh_lnbart3
#define lnbarh SMDR_Mh_lnbarh

/* Analytic cases: */
SMDR_COMPLEX AW, AZ, Ah, At, Ab,
  BWW, BZZ, Bhh, Btt, Bbb, Btautau,
  Ihhh, IttZ, Ihtt, I0tW, IWWZ, IhWW, IhZZ, 
  I00W, I00Z, I00h, I00t, I0hW, I0hZ, I0WZ, 
  UWW00, UZZ00, S00W, S00Z, TW00, TZ00, MWWWW0, 
  TW0W, Tbar0WW, Mtttt0, Tt0t, Tbar0tt, Tb0b, Tbar0bb, Mbbbb0;

/* These require calls to TSIL_Evaluate: */
SMDR_COMPLEX
  Mhhhhh, Uhhhh, Shhh,
  MttttZ, UtttZ, SttZ, TttZ, TZtt,
  Mtttth, Uttth, Shtt, Ttth, Thtt,
  MtZtZt, UZZtt,
  Mththt, Uhhtt,
  MtWtW0, UWW0t, Utt0W, S0tW, TW0t, Tt0W, Tbar0tW,
  MWWWWZ, UWWWZ, SWWZ, TWWZ, TZWW,
  MWWWWh, UWWWh, ShWW, TWWh, ThWW,
  MZZZZh, UZZZh, ShZZ, TZZh, ThZZ,
  MWZWZW, UZZWW,
  MhWhWW, UhhWW,
  MhZhZZ, UhhZZ;

/* These are needed for the 3-loop part: */
SMDR_REAL lnbart, lnbart2, lnbart3, lnbarh;

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/*  Computes Higgs pole mass at up to two loops. The argument
    loopOrder may take the following values:

    0    tree level
    1    1-loop
    1.5  1-loop plus 2-loop QCD corrections
    2    full 2-loop
    2.3  full 2-loop plus leading 3-loop QCD
    2.5  full 2-loop plus leading 3-loop QCD and nonQCD terms
*/

void SMDR_Eval_Mh_pole (SMDR_REAL Q_eval,
                        float loopOrder,
                        SMDR_REAL *Mhpoleresult,
                        SMDR_REAL *Gammahpoleresult)
{
  SMDR_COMPLEX se, CM2, last;
  SMDR_REAL PM_TOL;
  int i;
  char funcname[] = "SMDR_Eval_Mh_pole";

  PM_TOL = 1.0e-8;
  if ((SMDR_MHPOLE_TOLERANCE < PM_TOL) && (SMDR_MHPOLE_TOLERANCE > 1.0e-15))
     PM_TOL = SMDR_MHPOLE_TOLERANCE;

  /* If Q_eval is negative, then we just use the current Q.
     Otherwise, we run all the parameters from Q_in to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

  if ( (TSIL_FABS(loopOrder) > 0.0001) &&
       (TSIL_FABS(loopOrder-1) > 0.0001) &&
       (TSIL_FABS(loopOrder-1.5) > 0.0001) &&
       (TSIL_FABS(loopOrder-2) > 0.0001) &&
       (TSIL_FABS(loopOrder-2.3) > 0.0001) &&
       (TSIL_FABS(loopOrder-2.5) > 0.0001) )
    SMDR_Error (funcname, 
    "Invalid loop order specified, should be 0, 1, 1.5, 2, 2.3, or 2.5.", 3);

  /* Check input parameters for sanity: */
  SMDR_Check_lambda_Range (0.05, 0.2);
  SMDR_Check_VEV_Range (199, 301);
  SMDR_Check_g3_Range (0.5, 1.5);
  SMDR_Check_g_Range (0.4, 0.8);
  SMDR_Check_gp_Range (0.25, 0.5);
  SMDR_Check_yt_Range (0.5, 1.5);
  SMDR_Check_yb_Range (0.05);
  SMDR_Check_ytau_Range (0.02);
  SMDR_Check_Q_Range (49, 501);

  SMDR_Update ();

  /* If higher order, we iterate... */
  /* Starting value: */
  CM2 = h;

  for (i=0; i<MAXITERS; i++) {

    if (loopOrder < 0.999) break;

    last = CM2;

    /* Perform needed TSIL evaluations: */
    SMDR_Mh_DoTSILh (loopOrder, SMDR_CREAL(CM2));

    se = ONELOOPFACTOR * SMDR_Mh_Pi1h (SMDR_CREAL(CM2));

    if (loopOrder > 1.499)
      se += TWOLOOPFACTOR * SMDR_Mh_Pi2QCDh (SMDR_CREAL(CM2));

    if (loopOrder > 1.999)
      se += TWOLOOPFACTOR * SMDR_Mh_Pi2nonQCDh (SMDR_CREAL(CM2));

    if (loopOrder > 2.299)
      se += THREELOOPFACTOR * SMDR_Mh_Pi3QCDh ();

    if (loopOrder > 2.499)
      se += THREELOOPFACTOR * SMDR_Mh_Pi3nonQCDh ();

    /* New estimate is: */
    CM2 = h + se;

    if (TSIL_CABS((CM2 - last)/last) < PM_TOL) break;
  }

  if (i == MAXITERS)
    SMDR_Warn (funcname, "Too many iterations! Convergence not assured.");

  *Mhpoleresult = SMDR_SQRT(SMDR_CREAL(CM2));
  *Gammahpoleresult = -SMDR_CIMAG(CM2)/(*Mhpoleresult);

  return;
}


/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Higgs self coupling determination, given the Higgs pole mass. Does
   not set the model value of lambda to the determined value
   automatically. 

   The argument loopOrder has allowed values:

   0    tree level
   1    1-loop level
   1.5  1-loop plus 2-loop QCD
   2    2-loop level
   2.3  2-loop plus leading 3-loop QCD
   2.5  2-loop plus leading 3-loop QCD and yt
*/

SMDR_REAL SMDR_Eval_lambda (SMDR_REAL Q_eval,
                            SMDR_REAL mhpole, 
                            float loopOrder)
{
  SMDR_REAL m2hpole = mhpole*mhpole, lam, orig, last;
  SMDR_COMPLEX se;
  SMDR_REAL PM_TOL;
  int i;
  char funcname[] = "SMDR_Eval_lambda";

  PM_TOL = 1.0e-8;
  if ((SMDR_MHPOLE_TOLERANCE < PM_TOL) && (SMDR_MHPOLE_TOLERANCE > 1.0e-15))
     PM_TOL = SMDR_MHPOLE_TOLERANCE;

  /* If Q_eval is negative, then we just use the current Q.
     Otherwise, we run all the parameters from Q_in to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

  if ( (TSIL_FABS(loopOrder) > 0.0001) &&
       (TSIL_FABS(loopOrder-1) > 0.0001) &&
       (TSIL_FABS(loopOrder-1.5) > 0.0001) &&
       (TSIL_FABS(loopOrder-2) > 0.0001) &&
       (TSIL_FABS(loopOrder-2.3) > 0.0001) &&
       (TSIL_FABS(loopOrder-2.5) > 0.0001) )
    SMDR_Error (funcname, 
    "Invalid loopOrder specified, should be 0, 1, 1.5, 2, 2.3, or 2.5.", 3);

  /* Check input parameters for sanity: */
  SMDR_Check_Mhpole_Range (mhpole, 120, 130);
  SMDR_Check_VEV_Range (199, 301);
  SMDR_Check_g3_Range (0.5, 1.5);
  SMDR_Check_g_Range (0.4, 0.8);
  SMDR_Check_gp_Range (0.25, 0.5);
  SMDR_Check_yt_Range (0.5, 1.5);
  SMDR_Check_yb_Range (0.05);
  SMDR_Check_ytau_Range (0.02);
  SMDR_Check_Q_Range (49, 501);

  SMDR_Update ();

  /* Tree-level result: */
  lam = 0.5*m2hpole/v2;

  /* Save original value for later restoration: */
  orig = k;

  /* If higher order, we iterate... */
  for (i=0; i<MAXITERS; i++) {

    if (loopOrder < 0.001) break;
    last = lam;
    k = lam;
    SMDR_Update ();

    /* Perform needed TSIL evaluations: */
    SMDR_Mh_DoTSILh (loopOrder, m2hpole);

    se = ONELOOPFACTOR * SMDR_Mh_Pi1h (m2hpole);

    if (loopOrder > 1.499)
      se += TWOLOOPFACTOR * SMDR_Mh_Pi2QCDh (m2hpole);

    if (loopOrder > 1.999)
      se += TWOLOOPFACTOR * SMDR_Mh_Pi2nonQCDh (m2hpole);

    if (loopOrder > 2.299)
      se += THREELOOPFACTOR * SMDR_Mh_Pi3QCDh ();

    if (loopOrder > 2.499)
      se += THREELOOPFACTOR * SMDR_Mh_Pi3nonQCDh ();

    /* New estimate is: */
    lam = 0.5*(m2hpole - SMDR_CREAL(se))/v2;

    if (TSIL_FABS((lam - last)/last) < PM_TOL) break;
  }

  if (i == MAXITERS)
    SMDR_Warn (funcname, "Too many iterations! Convergence not assured.");

  /* Reset original value: */
  k = orig;
  SMDR_Update ();
  return(lam);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* One-loop self energy function for SM Higgs.                      */
/* 1407.4336 equation (2.46), plus bottom and tau contributions.    */

SMDR_COMPLEX SMDR_Mh_Pi1h (SMDR_REAL s)
{
  SMDR_REAL s2 = s*s;
  SMDR_REAL charm = yc * yc * v2/2.;
  SMDR_REAL strange = ys * ys * v2/2.;
  SMDR_REAL muon = ymu * ymu * v2/2.;
  SMDR_COMPLEX Bcc, Bss, Bmumu;

  Bcc = TSIL_B (charm, charm, s, Q2);
  Bss = TSIL_B (strange, strange, s, Q2);
  Bmumu = TSIL_B (muon, muon, s, Q2);

  return 3*yt2*(4*T - s)*Btt 
    + 3*yb2*(4*b - s)*Bbb 
    + ytau2*(4*tau - s)*Btautau 
    + 3*yc*yc*(4*charm - s)*Bcc
    + 3*ys*ys*(4*strange - s)*Bss
    + ymu*ymu*(4*muon - s)*Bmumu
    - 18*k2*v2*Bhh
    + 0.5*(g2 + gp2)*((s - 3*Z - s2/(4*Z))*BZZ - s*AZ/(2*Z) + 2*Z)
    + g2*((s - 3*W - s2/(4*W))*BWW - s*AW/(2*W) + 2*W);
}

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
/* Two-loop self energy function for SM Higgs (QCD contribution).       */
/* 1407.4336 equation (2.47), plus bottom, charm, strange contributions.*/

SMDR_COMPLEX SMDR_Mh_Pi2QCDh (SMDR_REAL s)
{
  SMDR_COMPLEX result;
  SMDR_COMPLEX lnbarms = SMDR_LOG(s/Q2) - I * PI;
  SMDR_COMPLEX smallmasspart;

  smallmasspart = s*(-131. + 68.*lnbarms - 12.*lnbarms*lnbarms + 48.*Zeta3);

  result = g32 * yt2 * (-17*s + 80*T - 28*At*At/T + (80*T - 36*s)*Btt + 
    (-176 + 36*s/T)*At*Btt + 14*s*Btt*Btt + (36*s - 168*T)*Tt0t + 
    16*(s - 4*T)*Tbar0tt + 8*(4*T - s)*(s - 2*T)*Mtttt0);

  if (b > 0.001) {
    result += g32 * yb2 * 
      (-17*s + 80*b - 28*Ab*Ab/b + (80*b - 36*s)*Bbb + 
      (-176 + 36*s/b)*Ab*Bbb + 14*s*Bbb*Bbb + (36*s - 168*b)*Tb0b + 
      16*(s - 4*b)*Tbar0bb + 8*(4*b - s)*(s - 2*b)*Mbbbb0);
  } else {
    result += g32 * yb2 * smallmasspart;
  }

  result += g32 * (yc * yc + ys * ys) * smallmasspart;

  return (result);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Two-loop self energy function for SM Higgs evaluated at the        */
/* specified squared momentum; non-QCD contributions.                 */
/* 1407.4336 ancillary file coefficients.txt and equation (2.48).     */

SMDR_COMPLEX SMDR_Mh_Pi2nonQCDh (SMDR_REAL s)
{
  SMDR_COMPLEX M2h2NONQCD,coMhhhhh,coUhhhh,coShhh,coIhhh,coMtttt0,coTt0t,
    coTbar0tt,coMttttZ,coUtttZ,coSttZ,coTttZ,coTZtt,coIttZ,coMtttth,coUttth,
    coShtt,coTtth,coThtt,coIhtt,coMtZtZt,coUZZtt,coMththt,coUhhtt,coMtWtW0,
    coUWW0t,coUtt0W,coS0tW,coTW0t,coTt0W,coI0tW,coMWWWWZ,coUWWWZ,coSWWZ,
    coTWWZ,coTZWW,coIWWZ,coMWWWW0,coTW0W,coTbar0WW,coMWWWWh,coUWWWh,coShWW,
    coTWWh,coIhWW,coMZZZZh,coUZZZh,coShZZ,coThZZ,coTZZh,coIhZZ,coMWZWZW,
    coUZZWW,coMhWhWW,coUhhWW,coMhZhZZ,coUhhZZ,coUWW00,coUZZ00,coS00W,coTW00,
    coS00Z,coTZ00,coI00W,coI00Z,coI00h,coI00t,coI0hW,coI0hZ,coI0WZ,coBWWBWW,
    coBZZBZZ,coBhhBhh,coBttBtt,coBZZBtt,coBWWBtt,coBhhBtt,coBWWBZZ,coBWWBhh,
    coBZZBhh,coAWBWW,coAZBWW,coAhBWW,coAtBWW,coAWBZZ,coAZBZZ,coAhBZZ,coAtBZZ,
    coAWBhh,coAZBhh,coAhBhh,coAtBhh,coAWBtt,coAZBtt,coAhBtt,coAtBtt,coAWAW,
    coAZAZ,coAhAh,coAtAt,coAZAt,coAWAt,coAhAt,coAWAZ,coAWAh,coAZAh,coBWW,
    coBZZ,coBhh,coBtt,coAW,coAZ,coAh,coAt,co1;

  coMhhhhh = -648*k4*v4;
 
  coUhhhh = 108*k3*v2;
 
  coShhh = 102*k2;
 
  coIhhh = -48*k2;
 
  coMtttt0 = (-16*g2*gp2*v4*yt2*(-2*k + yt2)*(-k + yt2))/(3*(g2 + gp2));
 
  coTt0t = v2*(4*g2*(6*k - 7*yt2)*yt2 - (4*g4*(6*k - 7*yt2)*yt2)/(g2 + gp2));
 
  coTbar0tt = (-32*g2*gp2*v2*yt2*(-k + yt2))/(3*(g2 + gp2));
 
  coMttttZ = v4*yt2*((-16*g4*(k - yt2)*(2*k - yt2))/(3*(g2 + gp2)) + 
      (256*g2*k2 - 64*gp2*k2 - 9*g4*yt2 + 6*g2*gp2*yt2 - 17*gp4*yt2 - 
        402*g2*k*yt2 + 78*gp2*k*yt2 + 146*g2*yt4 - 14*gp2*yt4)/24);
 
  coUtttZ = v2*(-((g2 + gp2)*(9*g4 - 30*g2*gp2 + 25*gp4)*(g2 + gp2 - 8*k))/
       (256*(g2 + gp2 - 8*yt2)) - (16*g4*(k - 2*yt2)*yt2)/(3*(g2 + gp2)) + 
      (27*g6 - 63*g4*gp2 - 15*g2*gp4 + 75*gp6 + 1512*g4*k - 432*g2*gp2*k + 
        2664*gp4*k + 504*g4*yt2 - 912*g2*gp2*yt2 + 1144*gp4*yt2 - 
        1088*g2*k*yt2 - 6208*gp2*k*yt2 - 9344*g2*yt4 + 896*gp2*yt4 + 
        9216*k*yt4)/768);
 
  coSttZ = (256*g4*k*yt2)/(3*(g4 + 2*g2*gp2 + gp4)) - 
     (9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 - 18*g4*yt2 - 132*g2*gp2*yt2 + 
       14*gp4*yt2)/(24*(g2 + gp2 - 2*k)) + 
     (8*g2*(8*g2*k - 37*g2*yt2 - 120*k*yt2))/(9*(g2 + gp2)) + 
     (32*g4*yt4)/(3*(g2 + gp2)*k) + (-207*g4 + 138*g2*gp2 - 391*gp4 - 
       368*g2*k + 272*gp2*k + 2674*g2*yt2 - 286*gp2*yt2 + 2400*k*yt2 + 
       108*yt4)/72 + (9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 + 18*g4*yt2 - 
       156*g2*gp2*yt2 + 82*gp4*yt2 - 584*g2*yt4 + 56*gp2*yt4)/(48*k);
 
  coTttZ = v2*(((9*g4 - 30*g2*gp2 + 25*gp4)*(15*g4 + 30*g2*gp2 + 15*gp4 - 
         56*g2*k - 56*gp2*k + 64*k2))/(768*(g2 + gp2 - 8*yt2)) - 
      (128*g4*k*(k - yt2)*yt2)/(3*(g4 + 2*g2*gp2 + gp4)) - 
      (yt2*(9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 - 18*g4*yt2 - 
         132*g2*gp2*yt2 + 14*gp4*yt2))/(48*(g2 + gp2 - 2*k)) + 
      (yt2*(9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 + 18*g4*yt2 - 
         156*g2*gp2*yt2 + 82*gp4*yt2 - 584*g2*yt4 + 56*gp2*yt4))/(96*k) - 
      (16*g2*(3*g2*k2 - 10*g2*k*yt2 - 30*k2*yt2 + 13*g2*yt4 + 30*k*yt4))/
       (9*(g2 + gp2)) + (16*g4*yt6)/(3*(g2 + gp2)*k) + 
      (-405*g6 + 945*g4*gp2 + 225*g2*gp4 - 1125*gp6 + 1512*g4*k - 
        5040*g2*gp2*k + 4200*gp4*k + 10560*g2*k2 - 4800*gp2*k2 - 
        4536*g4*yt2 + 11664*g2*gp2*yt2 - 11448*gp4*yt2 - 41536*g2*k*yt2 + 
        9664*gp2*k*yt2 - 38400*k2*yt2 + 59872*g2*yt4 - 6688*gp2*yt4 + 
        31488*k*yt4 + 6912*yt6)/2304);
 
  coTZtt = v2*(-((g2 + gp2)*(9*g4 - 30*g2*gp2 + 25*gp4)*(g2 + gp2 - 8*k))/
       (128*(g2 + gp2 - 8*yt2)) + (32*g4*(3*k - yt2)*yt2)/(3*(g2 + gp2)) - 
      ((g2 + gp2)*(9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 - 18*g4*yt2 - 
         132*g2*gp2*yt2 + 14*gp4*yt2))/(64*(g2 + gp2 - 2*k)) + 
      ((g2 + gp2 + 4*yt2)*(9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 - 18*g4*yt2 - 
         132*g2*gp2*yt2 + 14*gp4*yt2))/(384*k) + 
      (-315*g6 - 321*g4*gp2 - 529*g2*gp4 - 523*gp6 + 2520*g4*k + 
        48*g2*gp2*k + 4184*gp4*k - 1296*g4*yt2 + 3168*g2*gp2*yt2 - 
        3216*gp4*yt2 - 33408*g2*k*yt2 + 12672*gp2*k*yt2 + 14016*g2*yt4 - 
        1344*gp2*yt4)/1152);
 
  coIttZ = ((9*g4 - 30*g2*gp2 + 25*gp4)*k)/(24*(g2 + gp2 - 8*yt2)) - 
     (64*g4*k*yt2)/(3*(g4 + 2*g2*gp2 + gp4)) + 
     (9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 - 18*g4*yt2 - 132*g2*gp2*yt2 + 
       14*gp4*yt2)/(48*(g2 + gp2 - 2*k)) + 
     (36*g4 - 24*g2*gp2 + 68*gp4 + 55*g2*k - 25*gp2*k - 474*g2*yt2 + 
       6*gp2*yt2 - 56*k*yt2)/24 - (8*g2*(g2*k - 6*g2*yt2 - 10*k*yt2))/
      (3*(g2 + gp2)) - (16*g4*yt4)/(3*(g2 + gp2)*k) - 
     (9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 + 18*g4*yt2 - 156*g2*gp2*yt2 + 
       82*gp4*yt2 - 584*g2*yt4 + 56*gp2*yt4)/(96*k);
 
  coMtttth = 6*v4*yt4*(-k2 - 2*k*yt2 + yt4);
 
  coUttth = -12*v2*yt2*(-3*k2 + yt4);

  coShtt = (15*yt4)/2 - (12*yt6)/k;

  coTtth = v2*(3*(12*k + 5*yt2)*yt4 - (6*yt8)/k);

  coThtt = 0;
 
  coIhtt = -6*yt4 + (6*yt6)/k;
 
  coMtZtZt = (v4*yt2*(-9*g6 - 3*g4*gp2 - 11*g2*gp4 - 17*gp6 + 36*g4*k - 
       120*g2*gp2*k + 100*gp4*k + 18*g4*yt2 + 132*g2*gp2*yt2 - 14*gp4*yt2 - 
       144*g2*k*yt2 - 144*gp2*k*yt2 + 288*k2*yt2))/24;
 
  coUZZtt = v2*(((9*g4 - 30*g2*gp2 + 25*gp4)*(3*g4 + 6*g2*gp2 + 3*gp4 - 
         8*g2*k - 8*gp2*k + 16*k2))/(192*(g2 + gp2 - 8*yt2)) - 
      (128*g4*k2*yt2)/(3*(g4 + 2*g2*gp2 + gp4)) + 
      ((g2 + gp2)*(9*g6 + 3*g4*gp2 + 11*g2*gp4 + 17*gp6 - 18*g4*yt2 - 
         132*g2*gp2*yt2 + 14*gp4*yt2))/(96*(g2 + gp2 - 2*k)) - 
      (16*g2*k*(g2*k - 2*g2*yt2 - 10*k*yt2))/(3*(g2 + gp2)) + 
      (45*g6 + 87*g4*gp2 + 103*g2*gp4 + 61*gp6 - 72*g4*k - 144*g2*gp2*k - 
        72*gp4*k + 880*g2*k2 - 400*gp2*k2 - 72*g4*yt2 - 528*g2*gp2*yt2 + 
        56*gp4*yt2 - 1760*g2*k*yt2 + 800*gp2*k*yt2 - 896*k2*yt2)/192);
 
  coMththt = -36*k*v4*(3*k - 2*yt2)*yt4;
 
  coUhhtt = 18*k*v2*(8*k - 3*yt2)*yt2;
 
  coMtWtW0 = (-3*v4*(g2 - 2*yt2)*yt2*(g4 - 4*g2*k + g2*yt2 + 4*k*yt2))/4;
 
  coUWW0t = 3*v2*((g4*(g2 - 2*yt2)*(g2 + yt2))/(8*(g2 - 2*k)) + 
      (16*k2*yt4)/g2 + (g6 - 2*g4*k + g4*yt2 - 2*g2*k*yt2 + 16*k2*yt2 - 
        12*k*yt4)/2);
 
  coUtt0W = (3*v2*(6*g4*k + g4*yt2 - 6*g2*k*yt2 - 9*g2*yt4 + 20*k*yt4 - 
       10*yt6))/4;
 
  coS0tW = (-3*g2*(g2 - 2*yt2)*(g2 + yt2))/(2*(g2 - 2*k)) + 
     (3*(g2 - 2*yt2)*(g2 + yt2)*(g2 + 2*yt2))/(4*k) - (48*k*yt4)/g2 + 
     (-14*g4 + 16*g2*k + 8*g2*yt2 + 48*k*yt2 + 15*yt4)/2;
 
  coTW0t = v2*((-9*g4*(g2 - 2*yt2)*(g2 + yt2))/(16*(g2 - 2*k)) + 
      (3*g2*(g2 - 2*yt2)*(g2 + yt2)*(g2 + 2*yt2))/(32*k) + 
      (-13*g6 + 104*g4*k - 63*g4*yt2 + 168*g2*k*yt2 - 66*g2*yt4 - 192*k*yt4)/
       16);
 
  coTt0W = v2*yt2*((12*k*(4*k - yt2)*yt2)/g2 - (3*g2*(g2 - 2*yt2)*(g2 + yt2))/
       (8*(g2 - 2*k)) + (3*(g2 - 2*yt2)*(g2 + yt2)*(g2 + 2*yt2))/(16*k) + 
      (-8*g2*k + 5*g2*yt2 - 24*k*yt2 + 12*yt4)/4);
 
  coI0tW = (3*g2*(g2 - 2*yt2)*(g2 + yt2))/(4*(g2 - 2*k)) - 
     (3*(g2 - 2*yt2)*(g2 + yt2)*(g2 + 2*yt2))/(8*k) + 
     (3*(5*g4 - 7*g2*yt2 + 16*k*yt2 - 6*yt4))/4 + (24*k*yt4)/g2;
 
  coMWWWWZ = ((g4*(g2 - 4*k)*(3*g4 - 8*g2*k + 16*k2))/(4*(g2 + gp2)) + 
      (51*g8 - 18*g6*gp2 - g4*gp4 + 8*g6*k - 96*g4*gp2*k - 8*g2*gp4*k - 
        592*g4*k2 + 224*g2*gp2*k2 - 16*gp4*k2 + 768*g2*k3 - 256*gp2*k3)/64)*v4;
 
  coUWWWZ = (-(g4*(g2 - 4*k)*(3*g2 - 4*k))/(2*(g2 + gp2)) - 
      (3*g10)/(2*(g2 + gp2)*(g2 - 2*k)) - (g4*(51*g4 - 18*g2*gp2 - gp4))/
       (32*(g2 - 2*k)) - (2*gp4*k2)/g2 + (-21*g6 - g4*gp2 - 122*g4*k + 
        68*g2*gp2*k + 6*gp4*k + 480*g2*k2 - 240*gp2*k2)/8)*v2;
 
  coSWWZ = (4*g4*(14*g2 - 41*k))/(3*(g2 + gp2)) + 
     (6*g8)/((g2 + gp2)*(g2 - 2*k)) + (g2*(51*g4 - 18*g2*gp2 - gp4))/
      (8*(g2 - 2*k)) + ((3*g2 - gp2)*(33*g4 + 22*g2*gp2 + gp4))/
      (16*(g2 + gp2 - 2*k)) - (3*g8)/((g2 + gp2)*k) - 
     (201*g6 - 3*g4*gp2 - 21*g2*gp4 - gp6)/(32*k) + (2*gp4*k)/g2 - 
     (48*g6*k)/(g4 + 2*g2*gp2 + gp4) + (751*g4 - 202*g2*gp2 - 5*gp4 - 
       3280*g2*k + 976*gp2*k)/48;
 
  coTWWZ = ((3*g10)/((g2 + gp2)*(g2 - 2*k)) + (g4*(51*g4 - 18*g2*gp2 - gp4))/
       (16*(g2 - 2*k)) + (g2*(3*g2 - gp2)*(33*g4 + 22*g2*gp2 + gp4))/
       (64*(g2 + gp2 - 2*k)) - (3*g10)/(4*(g2 + gp2)*k) - 
      (g2*(201*g6 - 3*g4*gp2 - 21*g2*gp4 - gp6))/(128*k) - 
      (12*g6*(g2 - 2*k)*k)/(g4 + 2*g2*gp2 + gp4) + 
      (g4*(49*g4 - 194*g2*k + 96*k2))/(6*(g2 + gp2)) + 
      (463*g6 + 128*g4*gp2 + 9*g2*gp4 - 1952*g4*k + 308*g2*gp2*k + 36*gp4*k + 
        1344*g2*k2 + 96*gp2*k2)/48)*v2;
 
  coTZWW = ((g2*(3*g2 - gp2)*(33*g4 + 22*g2*gp2 + gp4))/(64*(g2 - 2*k)) + 
      (3*(3*g2 - gp2)*(g2 + gp2)*(33*g4 + 22*g2*gp2 + gp4))/
       (128*(g2 + gp2 - 2*k)) - ((3*g2 - gp2)*(3*g2 + gp2)*
        (33*g4 + 22*g2*gp2 + gp4))/(256*k) + (gp4*(gp2 - 8*k)*k)/(4*g2) + 
      (g4*(3*g4 - 20*g2*k + 16*k2))/(g2 + gp2) + 
      (567*g6 - 39*g4*gp2 - 31*g2*gp4 - gp6 - 2160*g4*k + 216*g2*gp2*k + 
        320*gp4*k + 1728*g2*k2 - 2304*gp2*k2)/96)*v2;
 
  coIWWZ = -((g4*(7*g2 - 8*k))/(g2 + gp2)) - (3*g8)/((g2 + gp2)*(g2 - 2*k)) - 
     (g2*(51*g4 - 18*g2*gp2 - gp4))/(16*(g2 - 2*k)) - 
     ((3*g2 - gp2)*(33*g4 + 22*g2*gp2 + gp4))/(32*(g2 + gp2 - 2*k)) + 
     (3*g8)/(2*(g2 + gp2)*k) + (201*g6 - 3*g4*gp2 - 21*g2*gp4 - gp6)/(64*k) - 
     (gp4*k)/g2 + (12*g6*k)/(g4 + 2*g2*gp2 + gp4) + 
     (-131*g4 + 50*g2*gp2 + gp4 + 400*g2*k - 208*gp2*k)/16;
 
  coMWWWW0 = (g2*gp2*(g2 - 4*k)*(3*g4 - 8*g2*k + 16*k2)*v4)/(4*(g2 + gp2));
 
  coTW0W = (11*g2*gp2*(3*g4 - 8*g2*k + 16*k2)*v2)/(3*(g2 + gp2));
 
  coTbar0WW = (g2*gp2*(3*g4 - 8*g2*k + 16*k2)*v2)/(g2 + gp2);
 
  coMWWWWh = -((3*g8 - 16*g6*k - 32*g4*k2 + 256*g2*k3 + 256*k4)*v4)/16;
 
  coUWWWh = ((3*g8)/(4*(g2 - 2*k)) + 2*k*(-g4 - 16*g2*k + 48*k2) - (128*k4)/g2)*
     v2;
 
  coShWW = (-3*g6)/(2*(g2 - 2*k)) + (3*g6)/(4*k) + (g4 + 96*g2*k - 528*k2)/12 + 
     (128*k3)/g2;
 
  coTWWh = ((-9*g8)/(2*(g2 - 2*k)) + (3*g8)/(16*k) + (23*g6 + 26*g4*k + 384*k3)/
       6)*v2;
 
  coIhWW = (5*g6)/(4*(g2 - 2*k)) - (3*g6)/(8*k) - 2*(3*g2 - 4*k)*k - (64*k3)/g2;
 
  coMZZZZh = -((3*g8 + 12*g6*gp2 + 18*g4*gp4 + 12*g2*gp6 + 3*gp8 - 16*g6*k - 
        48*g4*gp2*k - 48*g2*gp4*k - 16*gp6*k - 32*g4*k2 - 64*g2*gp2*k2 - 
        32*gp4*k2 + 256*g2*k3 + 256*gp2*k3 + 256*k4)*v4)/32;
 
  coUZZZh = ((3*(g8 + 4*g6*gp2 + 6*g4*gp4 + 4*g2*gp6 + gp8))/
       (8*(g2 + gp2 - 2*k)) + k*(-g4 - 2*g2*gp2 - gp4 - 16*g2*k - 16*gp2*k + 
        48*k2) - (64*k4)/(g2 + gp2))*v2;
 
  coShZZ = (-3*(g6 + 3*g4*gp2 + 3*g2*gp4 + gp6))/(4*(g2 + gp2 - 2*k)) + 
     (3*(g6 + 3*g4*gp2 + 3*g2*gp4 + gp6))/(8*k) + 
     (g4 + 2*g2*gp2 + gp4 + 96*g2*k + 96*gp2*k - 528*k2)/24 + 
     (64*k3)/(g2 + gp2);
 
  coThZZ = 0;
 
  coTZZh = ((-9*(g8 + 4*g6*gp2 + 6*g4*gp4 + 4*g2*gp6 + gp8))/
       (4*(g2 + gp2 - 2*k)) + (3*(g8 + 4*g6*gp2 + 6*g4*gp4 + 4*g2*gp6 + gp8))/
       (32*k) + (23*g6 + 69*g4*gp2 + 69*g2*gp4 + 23*gp6 + 26*g4*k + 
        52*g2*gp2*k + 26*gp4*k + 384*k3)/12)*v2;
 
  coIhZZ = (5*(g6 + 3*g4*gp2 + 3*g2*gp4 + gp6))/(8*(g2 + gp2 - 2*k)) - 
     (3*(g6 + 3*g4*gp2 + 3*g2*gp4 + gp6))/(16*k) - (3*g2 + 3*gp2 - 4*k)*k - 
     (32*k3)/(g2 + gp2);
 
  coMWZWZW = ((99*g8 + 33*g6*gp2 - 19*g4*gp4 - g2*gp6 - 312*g6*k + 
       44*g4*gp2*k - 16*g2*gp4*k - 4*gp6*k + 176*g4*k2 + 96*g2*gp2*k2 - 
       16*gp4*k2 - 256*g2*k3)*v4)/32;
 
  coUZZWW = (-((3*g2 - gp2)*(g2 + gp2)*(33*g4 + 22*g2*gp2 + gp4))/
       (64*(g2 + gp2 - 2*k)) - (2*g4*(3*g2 - 4*k)*k)/(g2 + gp2) + 
      (24*g6*k2)/(g4 + 2*g2*gp2 + gp4) + (-33*g6 + 11*g4*gp2 + 21*g2*gp4 + 
        gp6 + 38*g4*k + 12*g2*gp2*k - 18*gp4*k - 96*g2*k2)/16)*v2;
 
  coMhWhWW = (-3*k*(3*g6 - 12*g4*k + 64*k3)*v4)/2;
 
  coUhhWW = ((-27*g8)/(8*(g2 - 2*k)) + 
      (3*(9*g6 + 36*g4*k - 112*g2*k2 + 192*k3))/8)*v2;
 
  coMhZhZZ = (-3*k*(3*g6 + 9*g4*gp2 + 9*g2*gp4 + 3*gp6 - 12*g4*k - 
       24*g2*gp2*k - 12*gp4*k + 64*k3)*v4)/4;
 
  coUhhZZ = ((-27*(g8 + 4*g6*gp2 + 6*g4*gp4 + 4*g2*gp6 + gp8))/
       (16*(g2 + gp2 - 2*k)) + (3*(9*g6 + 27*g4*gp2 + 27*g2*gp4 + 9*gp6 + 
         36*g4*k + 72*g2*gp2*k + 36*gp4*k - 112*g2*k2 - 112*gp2*k2 + 192*k3))/
       16)*v2;
 
  coUWW00 = g4*((9*g4)/(8*(g2 - 2*k)) + (9*(g2 - 2*k))/2)*v2;
 
  coUZZ00 = (63*g4 + 6*g2*gp2 + 103*gp4)*
     ((g4 + 2*g2*gp2 + gp4)/(96*(g2 + gp2 - 2*k)) + (g2 + gp2 - 2*k)/24)*v2;
 
  coS00W = g2*((-3*(5*g2 - 16*k))/2 - (9*g4)/(2*(g2 - 2*k)) + (9*g4)/(4*k));
 
  coTW00 = ((-3*g4*(g2 - 8*k))/4 - (27*g8)/(16*(g2 - 2*k)) + (9*g8)/(32*k))*v2;
 
  coS00Z = -((g2 + gp2)*(63*g4 + 6*g2*gp2 + 103*gp4))/(24*(g2 + gp2 - 2*k)) + 
     ((g2 + gp2)*(63*g4 + 6*g2*gp2 + 103*gp4))/(48*k) + 
     (320*g4*k)/(9*(g2 + gp2)) + (-315*g4 - 30*g2*gp2 - 515*gp4 - 1552*g2*k + 
       1648*gp2*k)/72;
 
  coTZ00 = (63*g4 + 6*g2*gp2 + 103*gp4)*(-(g4 + 2*g2*gp2 + gp4)/
       (64*(g2 + gp2 - 2*k)) + (g4 + 2*g2*gp2 + gp4)/(384*k) + 
      (-g2 - gp2 + 8*k)/144)*v2;
 
  coI00W = (9*g6)/(4*(g2 - 2*k)) - (9*g6)/(8*k) - (6*g4*k)/(g2 + gp2) - 
     (3*g6*k)/(g4 + 2*g2*gp2 + gp4) + (3*g2*(3*g2 + 10*k))/2;
 
  coI00Z = ((g2 + gp2)*(63*g4 + 6*g2*gp2 + 103*gp4))/(48*(g2 + gp2 - 2*k)) - 
     ((g2 + gp2)*(63*g4 + 6*g2*gp2 + 103*gp4))/(96*k) + (9*g4*k)/(g2 + gp2) + 
     (63*g4 + 6*g2*gp2 + 103*gp4 - 144*g2*k)/24;
 
  coI00h = -12*k2;
 
  coI00t = -24*k*yt2;
 
  coI0hW = 3*k*(g2 + 8*k);
 
  coI0hZ = (3*k*(g2 + gp2 + 8*k))/2;
 
  coI0WZ = (9*g4*k)/(g2 + gp2) + 3*(4*g2 + gp2)*k + 
     (3*g6*k)/(g4 + 2*g2*gp2 + gp4);
 
  coBWWBWW = ((-3*g8)/(2*(g2 - 2*k)) + (3*g10)/(2*(g2 + gp2)*(g2 - 2*k)) - 
      (g4*(9*g4 - 16*g2*k + 64*k2))/(3*(g2 + gp2)) + 
      (183*g6 + 3*g4*gp2 - 328*g4*k + 24*g2*gp2*k + 1648*g2*k2 + 48*gp2*k2)/
       48)*v2;
 
  coBZZBZZ = ((-g6 - 3*g4*gp2 - 3*g2*gp4 - gp6 - 6*g4*k - 12*g2*gp2*k - 
       6*gp4*k + 64*g2*k2 + 64*gp2*k2 + 32*k3)*v2)/16;
 
  coBhhBhh = 54*k3*v2;
 
  coBttBtt = v2*yt2*((-4*g4*(7*k - 4*yt2))/(3*(g2 + gp2)) + 
      (56*g2*k - 23*g2*yt2 + 17*gp2*yt2 + 36*yt4)/6);
 
  coBZZBtt = -((9*g4 - 6*g2*gp2 + 17*gp4 - 18*g2*k - 18*gp2*k + 72*k2)*v2*yt2)/6;
 
  coBWWBtt = 3*v2*yt2*(-g4 + 2*g2*k - 8*k2 + g2*yt2 + 4*k*yt2);
 
  coBhhBtt = -36*k*v2*yt4;
 
  coBWWBZZ = ((15*g6 - 10*g4*gp2 - g2*gp4 - 12*g4*k - 48*g2*gp2*k - 4*gp4*k + 
       80*g2*k2 + 32*gp2*k2 - 64*k3)*v2)/8;
 
  coBWWBhh = (3*(g2 - 4*k)*k*(3*g2 + 8*k)*v2)/2;
 
  coBZZBhh = (3*(g2 + gp2 - 4*k)*k*(3*g2 + 3*gp2 + 8*k)*v2)/4;
 
  coAWBWW = (6*g8)/((g2 + gp2)*(g2 - 2*k)) - (g4*(77*g2 + gp2))/
      (8*(g2 - 2*k)) + (8*(gp2 + 8*k)*k2)/g2 - 
     (2*g2*(45*g4 - 112*g2*k + 304*k2))/(3*(g2 + gp2)) + 
     (282*g4 + 3*g2*gp2 - 694*g2*k - 6*gp2*k + 1888*k2)/6;
 
  coAZBWW = (-3*g8)/((g2 + gp2)*(g2 - 2*k)) + (g4*(9*g2 + gp2))/
      (8*(g2 - 2*k)) + (-19*g4 - g2*gp2 + 46*g2*k + 2*gp2*k - 192*k2)/2 - 
     (8*gp2*k2)/g2 + (4*g4*(3*g4 - 8*g2*k + 16*k2))/(g4 + 2*g2*gp2 + gp4) + 
     (2*g2*(7*g4 - 20*g2*k + 80*k2))/(g2 + gp2);
 
  coAhBWW = (7*g6)/(8*(g2 - 2*k)) - (3*g6)/(8*k) + (3*g4 - 22*g2*k + 48*k2)/2 - 
     (64*k3)/g2;
 
  coAtBWW = (3*g4*(g2 - yt2))/(2*(g2 - 2*k)) + (96*k2*yt2)/g2 + 
     6*(g2 - 2*k)*(g2 + yt2);
 
  coAWBZZ = -((g2 + gp2)*(15*g4 - 10*g2*gp2 - gp4))/(8*(g2 + gp2 - 2*k)) - 
     (8*g2*(3*g2 - 4*k)*k)/(g2 + gp2) + (33*g4 + 22*g2*gp2 + gp4 - 50*g2*k - 
       2*gp2*k)/2 + (96*g4*k2)/(g4 + 2*g2*gp2 + gp4);
 
  coAZBZZ = (-7*(g6 + 3*g4*gp2 + 3*g2*gp4 + gp6))/(8*(g2 + gp2 - 2*k)) + 
     (g4 + 2*g2*gp2 + gp4 + 16*g2*k + 16*gp2*k - 32*k2)/4 + (32*k3)/(g2 + gp2);
 
  coAhBZZ = (7*(g6 + 3*g4*gp2 + 3*g2*gp4 + gp6))/(16*(g2 + gp2 - 2*k)) - 
     (3*(g6 + 3*g4*gp2 + 3*g2*gp4 + gp6))/(16*k) + 
     (3*g4 + 6*g2*gp2 + 3*gp4 - 22*g2*k - 22*gp2*k + 48*k2)/4 - 
     (32*k3)/(g2 + gp2);
 
  coAtBZZ = ((g2 + gp2)*(9*g4 - 6*g2*gp2 + 17*gp4))/(12*(g2 + gp2 - 2*k)) + 
     (64*g2*k*(g2 + 5*k))/(3*(g2 + gp2)) + 
     (9*g4 + 66*g2*gp2 - 7*gp4 - 256*g2*k + 64*gp2*k - 400*k2)/12 - 
     (256*g4*k2)/(3*(g4 + 2*g2*gp2 + gp4)) + 
     ((9*g4 - 30*g2*gp2 + 25*gp4)*(3*g4 + 6*g2*gp2 + 3*gp4 - 8*g2*k - 
        8*gp2*k + 16*k2))/(12*(g2 + gp2)*(g2 + gp2 - 8*yt2));
 
  coAWBhh = (-27*g6)/(2*(g2 - 2*k)) + (3*(9*g4 + 64*k2))/2;
 
  coAZBhh = (-27*(g6 + 3*g4*gp2 + 3*g2*gp4 + gp6))/(4*(g2 + gp2 - 2*k)) + 
     (3*(9*g4 + 18*g2*gp2 + 9*gp4 + 64*k2))/4;
 
  coAhBhh = -108*k2;
 
  coAtBhh = 108*k*yt2;
 
  coAWBtt = 3*(6*g2*k - 3*g2*yt2 + 2*k*yt2 - 5*yt4);
 
  coAZBtt = -((9*g4 - 30*g2*gp2 + 25*gp4)*(g2 + gp2 - 8*k))/
      (32*(g2 + gp2 - 8*yt2)) + (128*g4*(k - yt2)*yt2)/
      (3*(g4 + 2*g2*gp2 + gp4)) + (27*g4 - 90*g2*gp2 + 75*gp4 - 888*g2*k + 
       1032*gp2*k - 216*g2*yt2 - 216*gp2*yt2 + 448*k*yt2 - 448*yt4)/96 + 
     (16*g2*(3*g2*k - 10*k*yt2 + 10*yt4))/(3*(g2 + gp2));
 
  coAhBtt = 3*(6*k - 5*yt2)*yt2 + (6*yt6)/k;
 
  coAtBtt = ((9*g4 - 30*g2*gp2 + 25*gp4)*(g2 + gp2 - 8*k))/
      (16*(g2 + gp2 - 8*yt2)) - (16*g4*(9*k - 8*yt2))/(3*(g2 + gp2)) + 
     (-27*g4 + 90*g2*gp2 - 75*gp4 + 2088*g2*k - 216*gp2*k - 2048*g2*yt2 - 
       192*gp2*yt2 - 1152*k*yt2 + 288*yt4)/48;
  
  coAWAW = (-2*(10*g2 - 2*gp2 - 203*k))/(3*v2) + (4*g2*(7*g2 - 64*k))/
      (3*(g2 + gp2)*v2) - (g2*(21*g2 + gp2))/(2*(g2 - 2*k)*v2) - 
     (48*g4*k)/((g4 + 2*g2*gp2 + gp4)*v2) + (4*k*(gp2 + 8*k))/(g2*v2);
 
  coAZAZ = ((-3*(g4 + 2*g2*gp2 + gp4))/(2*(g2 + gp2 - 2*k)) + 
      (4*g2 + 4*gp2 + 9*k)/3 + (16*k2)/(g2 + gp2))/v2;
 
  coAhAh = 0;
 
  coAtAt = ((128*g4*k)/(3*(g4 + 2*g2*gp2 + gp4)) - (8*g2*(g2 + 60*k))/
       (9*(g2 + gp2)) + (-83*g2 - 187*gp2 + 516*k - 324*yt2)/18 - 
      (2*(9*g4 - 30*g2*gp2 + 25*gp4)*k)/(3*(g2 + gp2)*(g2 + gp2 - 8*yt2)))/v2;
 
  coAZAt = (9*g4 - 6*g2*gp2 + 17*gp4)/(3*(g2 + gp2 - 2*k)*v2) - 
     (128*g4*k)/(3*(g4 + 2*g2*gp2 + gp4)*v2) + (16*g2*(11*g2 + 30*k))/
      (9*(g2 + gp2)*v2) - (253*g2 - 187*gp2 + 300*k - 54*yt2)/(18*v2) + 
     (2*(9*g4 - 30*g2*gp2 + 25*gp4)*k)/(3*(g2 + gp2)*v2*(g2 + gp2 - 8*yt2));
 
  coAWAt = (6*g2*(g2 - yt2))/((g2 - 2*k)*v2) + (48*k*yt2)/(g2*v2) + 
     (11*g2 + 9*yt2)/v2;
 
  coAhAt = (9*yt2)/v2;
 
  coAWAZ = (-12*g6)/((g2 + gp2)*(g2 - 2*k)*v2) + (g2*(9*g2 + gp2))/
      (2*(g2 - 2*k)*v2) - (15*g4 - 10*g2*gp2 - gp4)/(2*(g2 + gp2 - 2*k)*v2) - 
     (4*gp2*k)/(g2*v2) + (60*g4*k)/((g4 + 2*g2*gp2 + gp4)*v2) - 
     (4*(5*g2 + gp2 + 33*k))/(3*v2) + (4*g2*(5*g2 + 54*k))/(3*(g2 + gp2)*v2);
 
  coAWAh = (-4*(2*g2 - 9*k))/(3*v2) + (5*g4)/(2*(g2 - 2*k)*v2) - (32*k2)/(g2*v2);
 
  coAZAh = (-2*(2*g2 + 2*gp2 - 9*k))/(3*v2) + (5*(g4 + 2*g2*gp2 + gp4))/
      (4*(g2 + gp2 - 2*k)*v2) - (16*k2)/((g2 + gp2)*v2);
 
  coBWW = v2*((-3*g10)/(2*(g2 + gp2)*(g2 - 2*k)) + 
      (g4*(63*g4 - 160*g2*k + 400*k2))/(6*(g2 + gp2)) - 
      (g6*(41*g2 + 3*gp2 - 36*yt2))/(48*(g2 - 2*k)) + 
      (-457*g6 - 3*g4*gp2 + 692*g4*k - 12*g2*gp2*k - 1328*g2*k2 + 48*gp2*k2 + 
        384*k3 + 36*g4*yt2 + 144*g2*k*yt2 - 576*k2*yt2)/24);
 
  coBZZ = v2*((8*g4*k2*(9*g2 - 16*yt2))/(3*(g4 + 2*g2*gp2 + gp4)) - 
      ((9*g4 - 30*g2*gp2 + 25*gp4)*(3*g4 + 6*g2*gp2 + 3*gp4 - 8*g2*k - 
         8*gp2*k + 16*k2))/(192*(g2 + gp2 - 8*yt2)) - 
      ((g2 + gp2)*(113*g6 + 113*g4*gp2 + 143*g2*gp4 + 71*gp6 - 36*g4*yt2 + 
         24*g2*gp2*yt2 - 68*gp4*yt2))/(96*(g2 + gp2 - 2*k)) - 
      (2*g2*k*(9*g4 + 50*g2*k - 16*g2*yt2 - 80*k*yt2))/(3*(g2 + gp2)) + 
      (-793*g6 - 883*g4*gp2 - 1579*g2*gp4 - 1201*gp6 + 1288*g4*k + 
        816*g2*gp2*k + 2600*gp4*k + 3024*g2*k2 - 5040*gp2*k2 + 1536*k3 + 
        360*g4*yt2 + 336*g2*gp2*yt2 + 488*gp4*yt2 - 2048*g2*k*yt2 + 
        512*gp2*k*yt2 - 3200*k2*yt2)/192);
 
  coBhh = (-81*g6*v2)/16 - (81*g4*gp2*v2)/16 - (81*g2*gp4*v2)/16 - 
     (27*gp6*v2)/16 + (27*g8*v2)/(16*g2 + 16*gp2 - 32*k) + 
     (108*g6*gp2*v2)/(16*g2 + 16*gp2 - 32*k) + (162*g4*gp4*v2)/
      (16*g2 + 16*gp2 - 32*k) + (108*g2*gp6*v2)/(16*g2 + 16*gp2 - 32*k) + 
     (27*gp8*v2)/(16*g2 + 16*gp2 - 32*k) + (27*g8*v2)/(8*g2 - 16*k) - 
     (153*g4*k*v2)/4 - (51*g2*gp2*k*v2)/2 - (51*gp4*k*v2)/4 + 27*g2*k2*v2 + 
     9*gp2*k2*v2 - 216*k3*v2 - 108*k2*v2*yt2 + 54*k*v2*yt4;
 
  coBtt = v2*(((g2 + gp2)*(9*g4 - 30*g2*gp2 + 25*gp4)*(g2 + gp2 - 8*k))/
       (256*(g2 + gp2 - 8*yt2)) + (8*g4*(13*k - 12*yt2)*yt2)/(3*(g2 + gp2)) + 
      (-27*g6 + 63*g4*gp2 + 15*g2*gp4 - 75*gp6 + 216*g4*k - 720*g2*gp2*k + 
        600*gp4*k - 216*g4*yt2 + 720*g2*gp2*yt2 - 600*gp4*yt2 - 
        35264*g2*k*yt2 + 320*gp2*k*yt2 + 38400*g2*yt4 - 128*gp2*yt4 - 
        9216*k*yt4 + 9216*yt6)/768);
 
  coAW = (-2*g4*(35*g2 - 107*k))/(3*(g2 + gp2)) - 
     (6*g8)/((g2 + gp2)*(g2 - 2*k)) - ((3*g2 + gp2)*(23*g4 + gp4))/
      (16*(g2 + gp2 - 2*k)) + (3*g8)/((g2 + gp2)*k) + 
     (72*g6*k)/(g4 + 2*g2*gp2 + gp4) + (-100*g4 + 26*g2*gp2 + gp4 + 
       112*g2*k - 240*gp2*k + 192*k2 - 21*g2*yt2 - 288*k*yt2)/12 - 
     (2*k*(gp4 + 64*k2 - 12*yt4))/g2 - 
     (g2*(208*g4 - 51*g2*gp2 - 3*gp4 - 54*g2*yt2 + 36*yt4))/(24*(g2 - 2*k)) + 
     (129*g6 - 3*g4*gp2 - 21*g2*gp4 - gp6 - 12*g4*yt2 + 48*g2*yt4 + 48*yt6)/
      (32*k);
 
  coAZ = -(g2*(69*g4 - 16*g2*gp2 - gp4))/(16*(g2 - 2*k)) - (gp4*k)/g2 + 
     (4*g4*k*(9*g2 - 16*yt2))/(g4 + 2*g2*gp2 + gp4) - 
     ((9*g4 - 30*g2*gp2 + 25*gp4)*k)/(24*(g2 + gp2 - 8*yt2)) - 
     (497*g6 + 83*g4*gp2 - 49*g2*gp4 - 67*gp6 - 108*g4*yt2 + 360*g2*gp2*yt2 - 
       300*gp4*yt2)/(96*(g2 + gp2 - 2*k)) - 
     (105*g6 + 324*g4*k + 576*k3 - 104*g4*yt2 - 720*g2*k*yt2)/
      (9*(g2 + gp2)) + (-300*g4 + 48*g2*gp2 + 408*gp4 + 4659*g2*k - 
       2781*gp2*k + 576*k2 - 976*g2*yt2 + 64*gp2*yt2 - 1800*k*yt2 - 108*yt4)/
      72 + (g4*(9*g4 - 32*yt4))/(6*(g2 + gp2)*k) + 
     (387*g6 - 369*g4*gp2 - 519*g2*gp4 - 315*gp6 - 36*g4*yt2 + 
       312*g2*gp2*yt2 - 164*gp4*yt2 + 1168*g2*yt4 - 112*gp2*yt4)/(192*k);
 
  coAh = (-13*g4)/8 - (13*g2*gp2)/12 - (13*gp4)/24 + 
     g6/(16*g2 + 16*gp2 - 32*k) + (3*g4*gp2)/(16*g2 + 16*gp2 - 32*k) + 
     (3*g2*gp4)/(16*g2 + 16*gp2 - 32*k) + gp6/(16*g2 + 16*gp2 - 32*k) + 
     g6/(8*g2 - 16*k) - (9*g6)/(16*k) - (9*g4*gp2)/(16*k) - 
     (9*g2*gp4)/(16*k) - (3*gp6)/(16*k) - 11*g2*k - (11*gp2*k)/3 - 120*k2 - 
     (64*k3)/g2 - (32*k3)/(g2 + gp2) - (9*yt4)/2 + (6*yt6)/k;
 
  coAt = 6*g4 - (11*g2*gp2)/6 + (187*gp4)/36 - (9*g6)/(24*g2 + 24*gp2 - 48*k) - 
     (3*g4*gp2)/(24*g2 + 24*gp2 - 48*k) - (11*g2*gp4)/
      (24*g2 + 24*gp2 - 48*k) - (17*gp6)/(24*g2 + 24*gp2 - 48*k) - 
     (3*g6)/(4*g2 - 8*k) - (9*g6)/(16*k) - (g4*gp2)/(16*k) - 
     (11*g2*gp4)/(48*k) - (17*gp6)/(48*k) - (55*g2*k)/12 + (25*gp2*k)/12 + 
     (48*g4*k)/(9*g2 + 9*gp2) - (9*g4*k)/(12*g2 + 12*gp2 - 96*yt2) + 
     (30*g2*gp2*k)/(12*g2 + 12*gp2 - 96*yt2) - 
     (25*gp4*k)/(12*g2 + 12*gp2 - 96*yt2) - (991*g2*yt2)/36 + 
     (94*gp2*yt2)/9 + (304*g4*yt2)/(9*g2 + 9*gp2) - 
     (18*g4*yt2)/(24*g2 + 24*gp2 - 48*k) - (132*g2*gp2*yt2)/
      (24*g2 + 24*gp2 - 48*k) + (14*gp4*yt2)/(24*g2 + 24*gp2 - 48*k) + 
     (3*g4*yt2)/(4*g2 - 8*k) - (3*g4*yt2)/(4*k) + (13*g2*gp2*yt2)/(4*k) - 
     (41*gp4*yt2)/(24*k) - 86*k*yt2 + (1440*g2*k*yt2)/(9*g2 + 9*gp2) - 
     (128*g4*k*yt2)/(g4 + 2*g2*gp2 + gp4) - 3*yt4 - (6*g2*yt4)/(4*g2 - 8*k) + 
     (41*g2*yt4)/(3*k) - (7*gp2*yt4)/(6*k) + (24*k*yt4)/g2 - 
     (32*g4*yt4)/(3*g2*k + 3*gp2*k) + (27*yt6)/(2*k);
 
  co1 = v2*((3*g10)/(4*(g2 + gp2)*(g2 - 2*k)) - 
      ((9*g4 - 30*g2*gp2 + 25*gp4)*(9*g4 + 18*g2*gp2 + 9*gp4 - 20*g2*k - 
         20*gp2*k + 64*k2))/(768*(g2 + gp2 - 8*yt2)) - 
      (2*g4*k*(9*g4 - 18*g2*k + 32*k*yt2 - 32*yt4))/(g4 + 2*g2*gp2 + gp4) + 
      (757*g8 + 1487*g6*gp2 + 1859*g4*gp4 + 1361*g2*gp6 + 376*gp8 - 
        108*g6*yt2 - 36*g4*gp2*yt2 - 132*g2*gp4*yt2 - 204*gp6*yt2 + 
        72*g4*yt4 + 528*g2*gp2*yt4 - 56*gp4*yt4)/(192*(g2 + gp2 - 2*k)) + 
      (426*g8 - 2145*g6*k - 36*g4*k2 + 1728*k4 + 2200*g4*k*yt2 + 
        2880*g2*k2*yt2 - 1280*g4*yt4 - 2880*g2*k*yt4)/(36*(g2 + gp2)) - 
      (g4*(9*g6 - 64*yt6))/(12*(g2 + gp2)*k) + 
      (k*(gp6 - 2*gp4*k + 384*k3 + 48*k*yt4 - 48*yt6))/(4*g2) + 
      (g2*(1370*g6 + 57*g4*gp2 - 60*g2*gp4 - 3*gp6 - 216*g4*yt2 + 72*g2*yt4 + 
         144*yt6))/(192*(g2 - 2*k)) + (-201*g6 - 25539*g4*gp2 - 
        40587*g2*gp4 - 26817*gp6 - 23964*g4*k - 21624*g2*gp2*k + 
        49764*gp4*k + 70272*g2*k2 - 69888*gp2*k2 + 248832*k3 - 10584*g4*yt2 + 
        26544*g2*gp2*yt2 - 13288*gp4*yt2 - 105376*g2*k*yt2 + 
        36736*gp2*k*yt2 - 375552*k2*yt2 + 111872*g2*yt4 - 3008*gp2*yt4 + 
        168192*k*yt4 + 19008*yt6)/2304 - (1449*g8 - 1014*g4*gp4 - 
        840*g2*gp6 - 315*gp8 - 324*g6*yt2 + 252*g4*gp2*yt2 + 60*g2*gp4*yt2 - 
        300*gp6*yt2 + 144*g4*yt4 + 2304*g2*gp2*yt4 - 768*gp4*yt4 + 
        5536*g2*yt6 - 448*gp2*yt6 + 5184*yt8)/(768*k));

  /* Two-loop non-QCD contribution: */
  M2h2NONQCD = 
    coMhhhhh * Mhhhhh + 
    coUhhhh * Uhhhh +
    coShhh * Shhh +
    coIhhh * Ihhh +
    coMtttt0 * Mtttt0 +
    coTt0t * Tt0t +
    coTbar0tt * Tbar0tt +
    coMttttZ * MttttZ +
    coUtttZ * UtttZ +
    coSttZ * SttZ +
    coTttZ * TttZ +
    coTZtt * TZtt +
    coIttZ * IttZ +
    coMtttth * Mtttth +
    coUttth * Uttth +
    coShtt * Shtt +
    coTtth * Ttth +
    coThtt * Thtt +
    coIhtt * Ihtt +
    coMtZtZt * MtZtZt +
    coUZZtt * UZZtt +
    coMththt * Mththt +
    coUhhtt * Uhhtt +
    coMtWtW0 * MtWtW0 +
    coUWW0t * UWW0t +
    coUtt0W * Utt0W +
    coS0tW * S0tW +
    coTW0t * TW0t +
    coTt0W * Tt0W +
    coI0tW * I0tW +
    coMWWWWZ * MWWWWZ +
    coUWWWZ * UWWWZ +
    coSWWZ * SWWZ +
    coTWWZ * TWWZ +
    coTZWW * TZWW +
    coIWWZ * IWWZ +
    coMWWWW0 * MWWWW0 +
    coTW0W * TW0W +
    coTbar0WW * Tbar0WW +
    coMWWWWh * MWWWWh +
    coUWWWh * UWWWh +
    coShWW * ShWW +
    coTWWh * TWWh +
    coIhWW * IhWW +
    coMZZZZh * MZZZZh +
    coUZZZh * UZZZh +
    coShZZ * ShZZ +
    coThZZ * ThZZ +
    coTZZh * TZZh +
    coIhZZ * IhZZ +
    coMWZWZW * MWZWZW +
    coUZZWW * UZZWW +
    coMhWhWW * MhWhWW +
    coUhhWW * UhhWW +
    coMhZhZZ * MhZhZZ +
    coUhhZZ * UhhZZ +
    coUWW00 * UWW00 +
    coUZZ00 * UZZ00 +
    coS00W * S00W +
    coTW00 * TW00 +
    coS00Z * S00Z +
    coTZ00 * TZ00 +
    coI00W * I00W +
    coI00Z * I00Z +
    coI00h * I00h +
    coI00t * I00t +
    coI0hW * I0hW +
    coI0hZ * I0hZ +
    coI0WZ * I0WZ +
    coBWWBWW * BWW * BWW +
    coBZZBZZ * BZZ * BZZ +
    coBhhBhh * Bhh * Bhh +
    coBttBtt * Btt * Btt +
    coBZZBtt * BZZ * Btt +
    coBWWBtt * BWW * Btt +
    coBhhBtt * Bhh * Btt +
    coBWWBZZ * BWW * BZZ +
    coBWWBhh * BWW * Bhh +
    coBZZBhh * BZZ * Bhh +
    coAWBWW * AW * BWW +
    coAZBWW * AZ * BWW +
    coAhBWW * Ah * BWW +
    coAtBWW * At * BWW +
    coAWBZZ * AW * BZZ +
    coAZBZZ * AZ * BZZ +
    coAhBZZ * Ah * BZZ +
    coAtBZZ * At * BZZ +
    coAWBhh * AW * Bhh +
    coAZBhh * AZ * Bhh +
    coAhBhh * Ah * Bhh +
    coAtBhh * At * Bhh +
    coAWBtt * AW * Btt +
    coAZBtt * AZ * Btt +
    coAhBtt * Ah * Btt +
    coAtBtt * At * Btt +
    coAWAW * AW * AW +
    coAZAZ * AZ * AZ +
    coAhAh * Ah * Ah +
    coAtAt * At * At +
    coAZAt * AZ * At +
    coAWAt * AW * At +
    coAhAt * Ah * At +
    coAWAZ * AW * AZ +
    coAWAh * AW * Ah +
    coAZAh * AZ * Ah +
    coBWW * BWW +
    coBZZ * BZZ +
    coBhh * Bhh +
    coBtt * Btt +
    coAW * AW +
    coAZ * AZ +
    coAh * Ah +
    coAt * At +
    co1;

  return M2h2NONQCD;
}

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
/* Partial 3-loop self energy function for SM Higgs; leading QCD terms. */
/* 1407.4336 equation (3.3).                                            */

SMDR_COMPLEX SMDR_Mh_Pi3QCDh (void)
{
  return T * (g34 * yt2 * 
	      (248.12151804320070541 + 839.19661693776140186*lnbart
	       + 160*lnbart2 - 736*lnbart3) +
	      g32 * yt4 * (2764.3651243340150442 + 
			   1283.7156382854997277*lnbart - 
			   360*lnbart2 + 240*lnbart3));
}

/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
/* Partial 3-loop self energy function for H; leading non-QCD terms. */
/* 1407.4336 equation (3.4).                                         */

SMDR_COMPLEX SMDR_Mh_Pi3nonQCDh (void)
{
  return T * yt6 * (-3199.0165548150885 + 36.*lnbarh
		    - 2653.5107656974665*lnbart + 756.*lnbarh*lnbart
		    + 13.5*lnbart2 + 324.*lnbarh*lnbart2 - 225.*lnbart3);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Performs all needed basis function evaluations for the 2- and 3-
   loop terms. */

int SMDR_Mh_DoTSILh (float loopOrder, SMDR_REAL s)
{
  TSIL_DATA bar;
  int success = 1;

  /* Do all 1-loop functions, of course: */
  AW = TSIL_A (W, Q2);
  AZ = TSIL_A (Z, Q2);
  Bhh = TSIL_B (h, h, s, Q2);
  BWW = TSIL_B (W, W, s, Q2);
  BZZ = TSIL_B (Z, Z, s, Q2);
  Btt = TSIL_B (T, T, s, Q2);
  Bbb = TSIL_B (b, b, s, Q2);
  Btautau = TSIL_B (tau, tau, s, Q2);

  /* Then do 2-loop evals if needed... */
  if (loopOrder > 1.1) {

    /* These are needed for the QCD contribution: */
    At = TSIL_A (T, Q2);
    lnbart = At/T + 1;
    lnbart2 = lnbart*lnbart;
    lnbart3 = lnbart2*lnbart;

    success *= TSIL_Manalytic (T, T, T, T, 0, s, &Mtttt0);
    success *= TSIL_Tanalytic (T, 0, T, s, Q2, &Tt0t);
    success *= TSIL_Tbaranalytic (0, T, T, s, Q2, &Tbar0tt);
    if (b > 0.0001) {
      Ab = TSIL_A (b, Q2);
      success *= TSIL_Manalytic (b, b, b, b, 0, s, &Mbbbb0);
      success *= TSIL_Tanalytic (b, 0, b, s, Q2, &Tb0b);
      success *= TSIL_Tbaranalytic (0, b, b, s, Q2, &Tbar0bb);
    } else {
      /* Shouldn't actually matter. */
      Ab = 0; Mbbbb0 = 0; Tb0b = 0; Tbar0bb = 0; 
    }
  }

  if (loopOrder > 1.6) {

    /* These are needed for nonQCD part: */
    Ah = TSIL_A (h, Q2);
    lnbarh = Ah/h + 1;
    Ihhh = TSIL_I2 (h, h, h, Q2);
    IttZ = TSIL_I2 (T, T, Z, Q2);
    Ihtt = TSIL_I2 (T, T, h, Q2);
    I0tW = TSIL_I2 (0.0, T, W, Q2);
    IWWZ = TSIL_I2 (W, W, Z, Q2);
    IhWW = TSIL_I2 (W, W, h, Q2);
    IhZZ = TSIL_I2 (Z, Z, h, Q2);
    I00W = TSIL_I2 (0.0, 0.0, W, Q2);
    I00Z = TSIL_I2 (0.0, 0.0, Z, Q2);
    I00h = TSIL_I2 (0.0, 0.0, h, Q2);
    I00t = TSIL_I2 (0.0, 0.0, T, Q2);
    I0hW = TSIL_I2 (0.0, h, W, Q2);
    I0hZ = TSIL_I2 (0.0, h, Z, Q2);
    I0WZ = TSIL_I2 (0.0, W, Z, Q2);

    success *= TSIL_Manalytic (W, W, W, W, 0.0, s, &MWWWW0);
    success *= TSIL_Tanalytic (W, W, 0.0, s, Q2, &TW0W);
    success *= TSIL_Tbaranalytic (0.0, W, W, s, Q2, &Tbar0WW);
    success *= TSIL_Uanalytic (W, W, 0.0, 0.0, s, Q2, &UWW00);
    success *= TSIL_Uanalytic (Z, Z, 0.0, 0.0, s, Q2, &UZZ00);
    success *= TSIL_Sanalytic (0.0, 0.0, W, s, Q2, &S00W);
    success *= TSIL_Sanalytic (0.0, 0.0, Z, s, Q2, &S00Z);
    success *= TSIL_Tanalytic (W, 0.0, 0.0, s, Q2, &TW00);
    success *= TSIL_Tanalytic (Z, 0.0, 0.0, s, Q2, &TZ00);

    TSIL_SetParameters (&bar, h, Z, h, Z, Z, Q2);
    TSIL_Evaluate (&bar, s);
    MhZhZZ = TSIL_GetFunction (&bar, "M");
    UhhZZ  = TSIL_GetFunction (&bar, "Uxzuv");

    TSIL_SetParameters (&bar, h, W, h, W, W, Q2);
    TSIL_Evaluate (&bar, s);
    MhWhWW = TSIL_GetFunction (&bar, "M");
    UhhWW  = TSIL_GetFunction (&bar, "Uxzuv");

    TSIL_SetParameters (&bar, W, Z, W, Z, W, Q2);
    TSIL_Evaluate (&bar, s);
    MWZWZW = TSIL_GetFunction (&bar, "M");
    UZZWW  = TSIL_GetFunction (&bar, "Uyuzv");

    TSIL_SetParameters (&bar, Z, Z, Z, Z, h, Q2);
    TSIL_Evaluate (&bar, s);
    MZZZZh = TSIL_GetFunction (&bar, "M");
    UZZZh  = TSIL_GetFunction (&bar, "Uyuzv");
    ShZZ   = TSIL_GetFunction (&bar, "Suxv");
    TZZh   = TSIL_GetFunction (&bar, "Tuxv");
    ThZZ   = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, W, W, W, W, h, Q2);
    TSIL_Evaluate (&bar, s);
    MWWWWh = TSIL_GetFunction (&bar, "M");
    UWWWh  = TSIL_GetFunction (&bar, "Uyuzv");
    ShWW   = TSIL_GetFunction (&bar, "Suxv");
    TWWh   = TSIL_GetFunction (&bar, "Tuxv");
    ThWW   = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, W, W, W, W, Z, Q2);
    TSIL_Evaluate (&bar, s);
    MWWWWZ = TSIL_GetFunction (&bar, "M");
    UWWWZ  = TSIL_GetFunction (&bar, "Uyuzv");
    SWWZ   = TSIL_GetFunction (&bar, "Suxv");
    TWWZ   = TSIL_GetFunction (&bar, "Tuxv");
    TZWW   = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, T, W, T, W, 0.0, Q2);
    TSIL_Evaluate (&bar, s);
    MtWtW0  = TSIL_GetFunction (&bar, "M");
    UWW0t   = TSIL_GetFunction (&bar, "Uyuzv");
    Utt0W   = TSIL_GetFunction (&bar, "Uxzuv");
    S0tW    = TSIL_GetFunction (&bar, "Svyz");
    TW0t    = TSIL_GetFunction (&bar, "Tyzv");
    Tt0W    = TSIL_GetFunction (&bar, "Txuv");
    Tbar0tW = TSIL_GetFunction (&bar, "TBARvxu");

    TSIL_SetParameters (&bar, T, h, T, h, T, Q2);
    TSIL_Evaluate (&bar, s);
    Mththt = TSIL_GetFunction (&bar, "M");
    Uhhtt  = TSIL_GetFunction (&bar, "Uyuzv");

    TSIL_SetParameters (&bar, T, Z, T, Z, T, Q2);
    TSIL_Evaluate (&bar, s);
    MtZtZt = TSIL_GetFunction (&bar, "M");
    UZZtt  = TSIL_GetFunction (&bar, "Uyuzv");

    TSIL_SetParameters (&bar, T, T, T, T, h, Q2);
    TSIL_Evaluate (&bar, s);
    Mtttth = TSIL_GetFunction (&bar, "M");
    Uttth  = TSIL_GetFunction (&bar, "Uyuzv");
    Shtt   = TSIL_GetFunction (&bar, "Suxv");
    Ttth   = TSIL_GetFunction (&bar, "Tuxv");
    Thtt   = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, T, T, T, T, Z, Q2);
    TSIL_Evaluate (&bar, s);
    MttttZ = TSIL_GetFunction (&bar, "M");
    UtttZ  = TSIL_GetFunction (&bar, "Uyuzv");
    SttZ   = TSIL_GetFunction (&bar, "Suxv");
    TttZ   = TSIL_GetFunction (&bar, "Tuxv");
    TZtt   = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, h, h, h, h, h, Q2);
    TSIL_Evaluate (&bar, s);
    Mhhhhh = TSIL_GetFunction (&bar, "M");
    Uhhhh  = TSIL_GetFunction (&bar, "Uyuzv");
    Shhh   = TSIL_GetFunction (&bar, "Suxv");
  }

  /* Implement status codes eventually? */
  return success;
}

