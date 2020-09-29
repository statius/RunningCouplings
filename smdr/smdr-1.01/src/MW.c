/* W boson complex pole mass calculation from 1503.03782. */

#include "smdr_internal.h"

#define Nc 3
#define nQ 3
#define nL 3
#define Nc2 9
#define nL2 9

/* Local functions; not needed or useful elsewhere: */
SMDR_COMPLEX SMDR_MW_W1loop (void);
SMDR_COMPLEX SMDR_MW_W1loopALT (void);
SMDR_COMPLEX SMDR_MW_W2loopQCD (void);
SMDR_COMPLEX SMDR_MW_W2loopnonQCD (void);
int SMDR_MW_DoTSILW (float);

/* #define global variables used in this file, for convenience with safety. */
#define AW SMDR_MW_AW
#define AZ SMDR_MW_AZ
#define Ah SMDR_MW_Ah
#define At SMDR_MW_At
#define Ab SMDR_MW_Ab
#define Atau SMDR_MW_Atau
#define B0t SMDR_MW_B0t
#define Bbt SMDR_MW_Bbt
#define B0tau SMDR_MW_B0tau
#define B00 SMDR_MW_B00
#define BWZ SMDR_MW_BWZ
#define BhW SMDR_MW_BhW
#define B0h SMDR_MW_B0h
#define B0Z SMDR_MW_B0Z
#define Bht SMDR_MW_Bht
#define BtZ SMDR_MW_BtZ
#define I00h SMDR_MW_I00h
#define I00t SMDR_MW_I00t
#define I00W SMDR_MW_I00W
#define I00Z SMDR_MW_I00Z
#define I0hW SMDR_MW_I0hW
#define I0hZ SMDR_MW_I0hZ
#define I0tW SMDR_MW_I0tW
#define I0WZ SMDR_MW_I0WZ
#define Ihhh SMDR_MW_Ihhh
#define Ihtt SMDR_MW_Ihtt
#define IhWW SMDR_MW_IhWW
#define IhZZ SMDR_MW_IhZZ
#define IttZ SMDR_MW_IttZ
#define IWWZ SMDR_MW_IWWZ
#define M00tt0 SMDR_MW_M00tt0
#define M00000W SMDR_MW_M00000W
#define Tt00 SMDR_MW_Tt00
#define M0000Z SMDR_MW_M0000Z
#define M000W0 SMDR_MW_M000W0
#define M00tW0 SMDR_MW_M00tW0
#define M0WW0W SMDR_MW_M0WW0W
#define M0W0Z0 SMDR_MW_M0W0Z0
#define UhW00 SMDR_MW_UhW00
#define UZW00 SMDR_MW_UZW00
#define Tbar0WZ SMDR_MW_Tbar0WZ
#define Tbar0hW SMDR_MW_Tbar0hW
#define TW00 SMDR_MW_TW00
#define TZ00 SMDR_MW_TZ00
#define Th00 SMDR_MW_Th00
#define Th0t SMDR_MW_Th0t
#define Th0W SMDR_MW_Th0W
#define Tt0h SMDR_MW_Tt0h
#define Tt0Z SMDR_MW_Tt0Z
#define TZ0t SMDR_MW_TZ0t
#define TZ0W SMDR_MW_TZ0W
#define MWZZWW SMDR_MW_MWZZWW
#define UZWWZ SMDR_MW_UZWWZ
#define UWZWW SMDR_MW_UWZWW
#define SWZZ SMDR_MW_SWZZ
#define MWWZZh SMDR_MW_MWWZZh
#define UZWhW SMDR_MW_UZWhW
#define UWZhZ SMDR_MW_UWZhZ
#define ShWZ SMDR_MW_ShWZ
#define TZhW SMDR_MW_TZhW
#define ThWZ SMDR_MW_ThWZ
#define MhZWWZ SMDR_MW_MhZWWZ
#define UWhZZ SMDR_MW_UWhZZ
#define UhWWZ SMDR_MW_UhWWZ
#define MhWWZW SMDR_MW_MhWWZW
#define MhWWhW SMDR_MW_MhWWhW
#define UWhWW SMDR_MW_UWhWW
#define UhWhW SMDR_MW_UhWhW
#define ShhW SMDR_MW_ShhW
#define MhhWWh SMDR_MW_MhhWWh
#define UWhhh SMDR_MW_UWhhh
#define M0WWZW SMDR_MW_M0WWZW
#define M0WWhW SMDR_MW_M0WWhW
#define M0ZtW0 SMDR_MW_M0ZtW0
#define UWZ00 SMDR_MW_UWZ00
#define UZW0t SMDR_MW_UZW0t
#define U0t0W SMDR_MW_U0t0W
#define M0WtZt SMDR_MW_M0WtZt
#define UWZtt SMDR_MW_UWZtt
#define U0ttZ SMDR_MW_U0ttZ
#define M0Wtht SMDR_MW_M0Wtht
#define UWhtt SMDR_MW_UWhtt
#define U0tht SMDR_MW_U0tht
#define UhW0t SMDR_MW_UhW0t
#define SttW SMDR_MW_SttW
#define M0tW0t SMDR_MW_M0tW0t
#define UW0tt SMDR_MW_UW0tt
#define M00ttZ SMDR_MW_M00ttZ

SMDR_COMPLEX AW, AZ, Ah, At, Ab, Atau,
  B0t, Bbt, B0tau, B00, BWZ, BhW, B0h, B0Z, Bht, BtZ,
  I00h, I00t, I00W, I00Z, I0hW, I0hZ, I0tW, I0WZ,
  Ihhh, Ihtt, IhWW, IhZZ, IttZ, IWWZ,
  M00tt0, M00000W, Tt00,
  M0000Z, M000W0, M00tW0, M0WW0W, M0W0Z0,
  UhW00, UZW00,
  Tbar0WZ, Tbar0hW,
  TW00, TZ00,
  Th00, Th0t, Th0W, Tt0h, Tt0Z, TZ0t, TZ0W,
  MWZZWW, UZWWZ, UWZWW, SWZZ,
  MWWZZh, UZWhW, UWZhZ, ShWZ, TZhW, ThWZ,   
  MhZWWZ, UWhZZ, UhWWZ,
  MhWWZW,
  MhWWhW, UWhWW, UhWhW, ShhW,
  MhhWWh, UWhhh,
  M0WWZW, M0WWhW,
  M0ZtW0, UWZ00, UZW0t, U0t0W,
  M0WtZt, UWZtt, U0ttZ,
  M0Wtht, UWhtt, U0tht, UhW0t, SttW,
  M0tW0t, UW0tt,
  M00ttZ;

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/*  Computes W mass at up to two loops, from 1503.03782. 
    The argument
    loopOrder may take the folowing values:

    0    tree level
    1    1-loop
    1.5  1-loop plus 2-loop QCD corrections
    2    full 2-loop

    Results are used to set the global variables:

    SMDR_MW_pole,
    SMDR_GammaW_pole,
    SMDR_MW_BreitWigner,
    SMDR_GammaW_BreitWigner
*/

void SMDR_Eval_MW_pole (SMDR_REAL Q_eval,
                        float loopOrder, 
                        SMDR_REAL *MWpoleresult, 
                        SMDR_REAL *GammaWpoleresult,
                        SMDR_REAL *MWBreitWignerresult, 
                        SMDR_REAL *GammaWBreitWignerresult)
{
  SMDR_COMPLEX CM2W;
  char funcname[] = "SMDR_Eval_MW_pole";

  /* If Q_eval is negative, then we just use the current Q.
     Otherwise, we run all the parameters from Q_in to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

  if ( (TSIL_FABS(loopOrder) > 0.0001) &&
       (TSIL_FABS(loopOrder-1) > 0.0001) &&
       (TSIL_FABS(loopOrder-1.5) > 0.0001) &&
       (TSIL_FABS(loopOrder-2) > 0.0001) )
    SMDR_Error (funcname, 
    "Invalid loopOrder specified, should be 0, 1, 1.5, or 2", 3);

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

  /* Tree-level result: */
  CM2W = W;

  /* Perform needed TSIL evaluations: */
  SMDR_MW_DoTSILW (loopOrder);

  if (loopOrder > 0.999)
    CM2W += ONELOOPFACTOR * SMDR_MW_W1loop ();

  if (loopOrder > 1.499)
    CM2W += TWOLOOPFACTOR * SMDR_MW_W2loopQCD ();

  if (loopOrder > 1.999)
    CM2W += TWOLOOPFACTOR * SMDR_MW_W2loopnonQCD ();
  
  *MWpoleresult = SMDR_SQRT(SMDR_CREAL(CM2W));

  *GammaWpoleresult = -SMDR_CIMAG(CM2W)/(*MWpoleresult);

  *MWBreitWignerresult = SMDR_SQRT(
    (*MWpoleresult) * (*MWpoleresult) + 
    (*GammaWpoleresult) * (*GammaWpoleresult));

  *GammaWBreitWignerresult = (*GammaWpoleresult) * (*MWpoleresult)/
    (*MWBreitWignerresult);

  return;
}

/* ------------------------------------------------------------ */
/* From 1503.03782 equations (2.23) - (2.25)                    */

SMDR_COMPLEX SMDR_MW_W1loop ()
{
  SMDR_COMPLEX result;

  result = 
    (2*W*(-72*W2 + 3*h*Z + 128*W*Z + 3*Z2))/(9*v2*Z) + 
    ((3*W - h)*Ah)/(3*v2) + 
    ((48*W2 + h*Z - 36*W*Z + Z2)*AW)/(3*v2*Z) + 
    ((24*W2 - 8*W*Z - Z2)*AZ)/(3*v2*Z) + 
    ((4*h*W - h2 - 12*W2)*BhW)/(3*v2) + 
    ((4*W - Z)*(12*W2 + 20*W*Z + Z2)*BWZ)/(3*v2*Z);

  /* Light quark and lepton loops. */
  result += ((nL - 1) + Nc * (nQ - 1)) * (4*W2*(1 - 3*B00))/(9*v2);

  /* Top/bottom loop: */
  result += Nc * (Ab*(2*b - 2*T - 4*W) + At*(-2*b + 2*T - 4*W) 
    - 4*b*W - 4*T*W + (4*W2)/3 + Bbt*(2*b2 - 4*b*T + 2*T2 + 2*b*W 
    + 2*T*W - 4*W2))/(3*v2);

  /* Tau/nutau loop: */
  result += (Atau*(2*tau - 4*W) - 4*tau*W + B0tau*(2*tau2 + 2*tau*W - 4*W2) 
     + (4*W2)/3)/(3*v2);

  /* Alternate form for top/bottom loop with b = 0:
  result += Nc * ((4*W*(W - 3*T) + 6*(T - 2*W)*At + 
            6*(T2 + T*W - 2*W2)*B0t)/(9*v2));
  */

  return result;
}

/* ------------------------------------------------------------ */
/* From 1503.03782 equations (2.23) - (2.25)                    */

SMDR_COMPLEX SMDR_MW_W1loopALT ()
{
  SMDR_COMPLEX result, TMSbar;

  TMSbar = T;
  T = 173.34 * 173.34;
  T2 = T*T;
  SMDR_MW_DoTSILW(1);

  result = 
    (2*W*(-72*W2 + 3*h*Z + 128*W*Z + 3*Z2))/(9*v2*Z) + 
    ((3*W - h)*Ah)/(3*v2) + 
    ((48*W2 + h*Z - 36*W*Z + Z2)*AW)/(3*v2*Z) + 
    ((24*W2 - 8*W*Z - Z2)*AZ)/(3*v2*Z) + 
    ((4*h*W - h2 - 12*W2)*BhW)/(3*v2) + 
    ((4*W - Z)*(12*W2 + 20*W*Z + Z2)*BWZ)/(3*v2*Z);

  result += (nL + Nc * (nQ - 1)) * (4*W2*(1 - 3*B00))/(9*v2);

  result += Nc * ((4*W*(W - 3*T) + 6*(T - 2*W)*At + 
            6*(T2 + T*W - 2*W2)*B0t)/(9*v2));

  /* Now put in the terms in the expansion of TMSbar around T. */
  result += (6*(-T + TMSbar)*(At - 2*W + B0t*(T + W)))/v2; 
/*
  result += (6*(-T + TMSbar)*(-T + TMSbar)*(At + B0t*T - W))/(v2*(T - W)); 
  result += (2*(B0t + At/T)*(-T + TMSbar)*(-T + TMSbar)*(-T + TMSbar))/
            (v2*(T - W));
  result += -((-T + TMSbar)*(-T + TMSbar)*(-T + TMSbar)*(-T + TMSbar)*W)/
            (2*T*v2*(T - W)*(T - W)); 
*/

  T = TMSbar;
  T2 = T*T;
  SMDR_MW_DoTSILW(1);

  return result;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Two-loop W boson pole mass QCD contribution.                        */
/* From 1503.03782 equations (2.26)                                    */

SMDR_COMPLEX SMDR_MW_W2loopQCD (void)
{
  SMDR_COMPLEX result;

  result = Nc * (g32/v2) * (
    (2*(17*T - 78*W)*W)/9 + (8*(5*T + 14*W)*At)/3 - 
    (32*(5*T + 3*W)*At*At)/(9*T) + (8*(5*T2 + 7*T*W - 2*W2)*B0t)/3 - 
    (16*(9*T2 + 14*T*W + 4*W2)*At*B0t)/(9*T) - 
    (8*W*(5*T + 4*W)*B0t*B0t)/9 - (16*(T - W)*(T - W)*(T + 2*W)*M00tt0)/9 + 
    (32*(T - 2*W)*(T + W)*Tt00)/9);

  result += Nc * (nQ - 1) * (g32/v2) * 
            W2 * (-124./9. - (16*B00)/3 - (32*W*M00000W)/9);  

  return (result);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/*  Two-loop non-QCD contribution to the W boson pole mass.           */ 
/* From 1503.03782 ancillary file "coefficients.txt" and eq. (2.27)   */

SMDR_COMPLEX SMDR_MW_W2loopnonQCD (void) 
{
SMDR_COMPLEX result;

SMDR_REAL co1, coAh, coAhAh, 
 coAhAt, coAhAW, coAhAZ, coAhB00, coAhB0t, coAhBhW, coAhBWZ, coAt, coAtAt, 
 coAtAW, coAtAZ, coAtB00, coAtB0t, coAtBhW, coAtBWZ, coAW, coAWAW, coAWAZ, 
 coAWB00, coAWB0t, coAWBhW, coAWBWZ, coAZ, coAZAZ, coAZB00, coAZB0t, coAZBhW, 
 coAZBWZ, coB00, coB00B00, coB00B0t, coB00BhW, coB00BWZ, coB0h, coB0t, 
 coB0tB0t, coB0tBhW, coB0tBWZ, coB0Z, coBht, coBhW, coBhWBhW, coBhWBWZ, coBtZ, 
 coBWZ, coBWZBWZ, coI00h, coI00t, coI00W, coI00Z, coI0hW, coI0hZ, coI0tW, 
 coI0WZ, coIhhh, coIhtt, coIhWW, coIhZZ, coIttZ, coIWWZ, coM00000, coM0000Z, 
 coM000W0, coM00tt0, coM00ttZ, coM00tW0, coM0tW0t, coM0W0Z0, coM0Wtht, 
 coM0WtZt, coM0WW0W, coM0WWhW, coM0WWZW, coM0ZtW0, coMhhWWh, coMhWWhW, 
 coMhWWZW, coMhZWWZ, coMWWZZh, coMWZZWW, coShhW, coShWZ, coSttW, coSWZZ, 
 coTbar0hW, coTbar0WZ, coTh00, coTh0t, coTh0W, coThWZ, coTt00, coTt0h, coTt0Z, 
 coTW00, coTZ00, coTZ0t, coTZ0W, coTZhW, coU0t0W, coU0tht, coU0ttZ, coUhW00, 
 coUhW0t, coUhWhW, coUhWWZ, coUW0tt, coUWhhh, coUWhtt, coUWhWW, coUWhZZ, 
 coUWZ00, coUWZhZ, coUWZtt, coUWZWW, coUZW00, coUZW0t, coUZWhW, coUZWWZ;

co1 = (64*nL2*W3)/(81*v4) + (128*Nc*nL*(nQ - 1)*W3)/(81*v4) +
     (64*Nc2*(nQ - 1)*(nQ - 1)*W3)/(81*v4) +
     (Nc2*((-8*T2*W)/9 - (56*T*W2)/27 + (64*W3)/81))/v4 +
     (Nc*nL*((-56*T*W2)/27 + (128*W3)/81))/v4 +
     (Nc2*(nQ - 1)*((-56*T*W2)/27 + (128*W3)/81))/v4 +
     (Nc*(nQ - 1)*((-7*h3)/144 + (7*h4)/(144*(h - 4*W)) - (7*h2*W)/36 +
        (43*h*W2)/27 - (17081*W3)/243 - (10880*W5)/(81*Z2) +
        (56512*W4)/(243*Z) + (10300*W2*Z)/243 - (989*W*Z2)/162 +
        (8*Z3)/9))/v4 + (nL*((-7*h3)/144 + (7*h4)/(144*(h - 4*W)) -
        (7*h2*W)/36 + (43*h*W2)/27 - (20515*W3)/81 - (2176*W5)/(9*Z2) +
        (13088*W4)/(27*Z) + (2644*W2*Z)/27 - (61*W*Z2)/2 + 8*Z3))/v4 +
     (Nc*((-7*h3)/144 + (31*h2*T)/32 - 4*h*T2 - (434*T3)/9 +
        (h*(14*h3 + 9*h2*T + 96*h*T2 - 1344*T3))/(288*(h - 4*W)) -
        (2*h2*T2)/(3*W) - (2*h*T3)/(3*W) + (65*h2*W)/36 -
        (437*h*T*W)/72 - (3137*T2*W)/324 + (24*T3*W)/h - (11*h*W2)/27 -
        (3982*T*W2)/81 - (16*T2*W2)/h + (42149*W3)/486 - (64*W4)/(9*T) +
        (4096*T2*W4)/(27*Z3) - (4096*T*W5)/(27*Z3) -
        (15872*T2*W3)/(81*Z2) + (5120*T*W4)/(81*Z2) -
        (8000*W5)/(81*Z2) + (4864*T2*W2)/(81*Z) + (22688*T*W3)/(81*Z) +
        (17128*W4)/(243*Z) - (307*T2*Z)/54 - (2*T3*Z)/(3*W) -
        (445*T*W*Z)/108 - (13861*W2*Z)/486 + (160*W3*Z)/(9*T) +
        (49*T*Z2)/8 - (2*T2*Z2)/(3*W) + (161*W*Z2)/324 -
        (140*W2*Z2)/(9*T) + (737*Z3)/288 + (44*W*Z3)/(9*T) +
        (2*T2*Z*(2*T + Z)*(2*T + Z))/(9*(T - W)*(T - W)) +
        (2*T*(255*T3 - 40*T2*Z - 37*T*Z2 - 7*Z3))/(27*(T - W)) +
        ((-8*W + 5*Z)*(-8*W + 5*Z)*(-112*W3 - 112*W2*Z + 16*W*Z2 + 3*Z3))/
         (96*(4*T - Z)*Z)))/v4 + ((-19577*h3)/4608 + (2*h4)/W +
       (46799*h2*W)/1152 - (97579*h*W2)/864 + (236725*W3)/81 + (4*W4)/h -
       (9*h5)/(512*(h - 4*W)*(h - 4*Z)) -
       (h*(157*h3 + 340*h2*W - 6704*h*W2 + 24576*W3))/(1536*(h - 4*Z)) +
       (3576*W6)/Z3 - (4*h3*W2)/Z2 - (3*h2*W3)/Z2 +
       (7528*W5)/(3*Z2) - (47*h4)/(64*Z) - (397*h3*W)/(48*Z) -
       (47*h2*W2)/(36*Z) + (59*h*W3)/(3*Z) - (21274*W4)/(3*Z) -
       (5977*h2*Z)/1152 - (h3*Z)/(3*W) + (2383*h*W*Z)/288 -
       (309515*W2*Z)/216 + (437*h*Z2)/144 - (2*h2*Z2)/(3*W) +
       (2219*W*Z2)/9 + (20*W2*Z2)/h + (1525*Z3)/36 - (h*Z3)/(3*W) -
       (18*W*Z3)/h + (2*Z4)/W + (h*(3384*h4 + 13601*h3*Z + 2436*h2*Z2 -
          9888*h*Z3 + 12288*Z4))/(4608*(h - 4*W)*Z))/v4;
 
coAh = (nL*(-h2/36 + h3/(36*(h - 4*W)) - (11*h*W)/27 + (10*W2)/9))/
      v4 + (Nc*(nQ - 1)*(-h2/36 + h3/(36*(h - 4*W)) - (11*h*W)/27 +
        (10*W2)/9))/v4 + (Nc*((335*h2)/36 - (106*h*T)/9 + (59*T2)/6 +
        (h*(13*h2 - 36*h*T + 48*T2))/(36*(h - 4*W)) -
        h4/(4*(h - 4*W)*(h + T - W)) + (4*h3)/(3*W) - (4*h2*T)/(3*W) +
        (4*h*T2)/(3*W) - (1001*h*W)/27 + (200*T*W)/9 - (4*T2*W)/h +
        (229*W2)/9 - (16*h4 + 97*h3*W - 552*h2*W2 + 732*h*W3 - 292*W4)/
         (12*(h + T - W)*W)))/v4 + ((1565*h2)/288 - (2*h3)/W -
       (3131*h*W)/432 + (1253*W2)/18 - (206*W3)/(3*h) -
       (h*(2*h2 - 41*h*W + 384*W2))/(48*(h - 4*Z)) + (4*h2*W2)/Z2 +
       (8*h*W3)/Z2 + (3*h3)/(4*Z) + (17*h2*W)/(3*Z) - (176*h*W2)/(9*Z) +
       (320*W3)/(9*Z) + (320*W4)/(3*h*Z) + (19*h*Z)/36 + (h2*Z)/(3*W) +
       (841*W*Z)/36 - (22*Z2)/9 + (h*Z2)/(3*W) + (11*W*Z2)/h -
       (h*(216*h3 + 1489*h2*Z - 408*h*Z2 + 192*Z3))/(288*(h - 4*W)*Z))/v4;
 
coAt = (Nc*nL*((16*T*W)/27 - (80*W2)/27))/v4 +
     (Nc2*(nQ - 1)*((16*T*W)/27 - (80*W2)/27))/v4 +
     (Nc2*((20*T2)/9 + (40*T*W)/27 - (80*W2)/27))/v4 +
     (Nc*((-21*h2)/2 + (50*h*T)/9 + (6650*T2)/27 -
        (h*(h2 - 6*h*T - 8*T2))/(2*(h - 4*W)) +
        h4/(3*(h - 4*W)*(h + T - W)) - (4*h3)/(3*W) + (4*h2*T)/(3*W) +
        (193*h*W)/9 + (18926*T*W)/81 - (24*T2*W)/h - (4156*W2)/81 -
        (16*T*W2)/h + (64*W3)/(9*T) + (h*(4*h3 + 28*h2*W - 110*h*W2 +
           73*W3))/(3*(h + T - W)*W) - (4096*T*W4)/(27*Z3) -
        (4096*W5)/(27*Z3) + (15872*T*W3)/(81*Z2) - (4096*W4)/(81*Z2) -
        (16*T3)/(9*Z) - (16*T2*W)/(9*Z) - (4432*T*W2)/(81*Z) +
        (33152*W3)/(81*Z) - (16*W4)/(9*T*Z) + (127*T*Z)/9 - (83*W*Z)/9 -
        (281*W2*Z)/(9*T) - (871*Z2)/36 + (4*T*Z2)/(3*W) -
        (86*W*Z2)/(3*T) + (155*Z3)/(9*T) - (4*Z3)/(3*W) +
        (4*Z4)/(3*T*W) - (4*T2*(2*T + Z)*(2*T + Z))/(9*(T - W)*(T - W)) +
        (4*T*(12*T3 - 1556*T2*Z + 43*T*Z2 + 7*Z3))/(27*(T - W)*Z) -
        ((-8*W + 5*Z)*(-8*W + 5*Z)*(192*W3 + 464*W2*Z + 136*W*Z2 + 9*Z3))/
         (108*(4*T - Z)*Z2) + (16*W6 - 80*W5*Z + 473*W4*Z2 -
          215*W3*Z3 - 325*W2*Z4 + 143*W*Z5 + 12*Z6)/
         (9*T*W*(-T + W - Z)*Z)))/v4;
 
coAW = (Nc*(nQ - 1)*(h2/18 - h3/(18*(h - 4*W)) + (26*h*W)/27 -
        (8992*W2)/243 + (7040*W4)/(81*Z2) - (13120*W3)/(243*Z) +
        (1808*W*Z)/243 + (22*Z2)/81))/v4 +
     (nL*(h2/18 - h3/(18*(h - 4*W)) + (26*h*W)/27 + (320*W2)/27 +
        (1408*W4)/(9*Z2) - (4928*W3)/(27*Z) + (464*W*Z)/27 + (2*Z2)/3))/
      v4 + (Nc*(h2/18 + (23*h*T)/18 - (443*T2)/9 -
        (h*(h2 - 21*h*T + 72*T2))/(18*(h - 4*W)) + (52*T3)/(T - W) +
        (26*h*W)/27 + (4870*T*W)/81 + (16*T2*W)/h - (24106*W2)/243 +
        (4096*T*W4)/(27*Z3) - (12032*T*W3)/(81*Z2) +
        (3968*W4)/(81*Z2) - (64*T2*W)/(9*Z) - (9536*T*W2)/(81*Z) +
        (12320*W3)/(243*Z) + (2*T*Z)/9 + (5414*W*Z)/243 + (221*Z2)/162 +
        ((-8*W + 5*Z)*(-8*W + 5*Z)*(12*W2 + 20*W*Z + Z2))/(54*(4*T - Z)*Z)))/v4 +
     ((-3451*h2)/96 + h3/W + (24577*h*W)/216 - (40271*W2)/27 +
       (20*W3)/h - (h2*(h - 12*W))/(32*(h - 4*Z)) -
       (3*h4)/(32*(h - 4*W)*(h - 4*Z)) - (7664*W5)/(3*Z3) +
       (4*h2*W2)/Z2 - (2944*W4)/(3*Z2) + (3*h3)/(4*Z) + (43*h2*W)/Z -
       (556*h*W2)/(3*Z) + (41908*W3)/(9*Z) + (17*h*Z)/72 + (h2*Z)/(3*W) -
       (24745*W*Z)/54 - (152*Z2)/9 + (h*Z2)/(3*W) - (20*W*Z2)/h + Z3/W -
       (h*(216*h3 - 149*h2*Z - 12*h*Z2 - 1536*Z3))/(288*(h - 4*W)*Z))/v4;
 
coAZ = (nL*((250*W2)/3 + (320*W4)/(9*Z2) - (2080*W3)/(27*Z) -
        (704*W*Z)/27 - (8*Z2)/3))/v4 +
     (Nc*(nQ - 1)*((6682*W2)/243 + (1600*W4)/(81*Z2) -
        (5984*W3)/(243*Z) - (2048*W*Z)/243 - (88*Z2)/81))/v4 +
     ((83*h2)/16 - (871*h*W)/24 + (37027*W2)/27 +
       (9*h4)/(16*(h - 4*W)*(h - 4*Z)) - (h*(7*h2 - 22*h*W + 912*W2))/
        (48*(h - 4*Z)) + (640*W5)/(3*Z3) + (4*h2*W2)/Z2 -
       (8*h*W3)/Z2 + (928*W4)/Z2 + h3/(2*Z) + (14*h2*W)/(3*Z) -
       (184*h*W2)/(9*Z) + (25304*W3)/(9*Z) - (101*h*Z)/36 + (h2*Z)/(3*W) -
       (13283*W*Z)/54 + (12*W2*Z)/h - 43*Z2 + (h*Z2)/(3*W) +
       (18*W*Z2)/h - (2*Z3)/W - (h*(24*h3 + 61*h2*Z - 12*h*Z2 +
          128*Z3))/(48*(h - 4*W)*Z))/v4 +
     (Nc*((775*T2)/54 + (512*T*W)/81 - (2039*W2)/243 +
        (2048*T*W4)/(27*Z3) - (16*T4)/(9*Z2) - (16*T3*W)/(9*Z2) -
        (16*T2*W2)/(9*Z2) - (8080*T*W3)/(81*Z2) + (304*W4)/(81*Z2) -
        (232*T3)/(27*Z) - (40*T2*W)/(27*Z) + (3488*T*W2)/(81*Z) +
        (352*W3)/(243*Z) - (1187*T*Z)/54 + (4*T2*Z)/(3*W) -
        (9581*W*Z)/243 + (10381*Z2)/648 - (4*T*Z2)/(3*W) + (4*Z3)/(3*W) -
        (4*T3*(2*T + Z)*(2*T + Z))/(9*(T - W)*(T - W)*Z) -
        ((-8*W + 5*Z)*(-8*W + 5*Z)*(48*W2 - 40*W*Z + Z2))/(216*(4*T - Z)*Z) +
        (2*T2*(2*T + Z)*(12*T2 + 64*T*Z + 17*Z2))/(27*(T - W)*Z2) +
        (16*W6 - 80*W5*Z + 473*W4*Z2 - 215*W3*Z3 - 325*W2*Z4 +
          143*W*Z5 + 12*Z6)/(9*W*(-T + W - Z)*Z2)))/v4;
 
coB00 = (-80*nL2*W3)/(27*v4) - (160*Nc*nL*(nQ - 1)*W3)/
      (27*v4) - (80*Nc2*(nQ - 1)*(nQ - 1)*W3)/(27*v4) +
     (Nc*nL*((8*T*W2)/9 - (80*W3)/27))/v4 +
     (Nc2*(nQ - 1)*((8*T*W2)/9 - (80*W3)/27))/v4 +
     (Nc*(nQ - 1)*(-h3/48 + h4/(48*(h - 4*W)) - (h2*W)/12 - (7*h*W2)/9 -
        (1804*W3)/27 + (1216*W4)/(27*Z) - (8*W2*Z)/27))/v4 +
     (nL*(-h3/48 + h4/(48*(h - 4*W)) - (h2*W)/12 - (7*h*W2)/9 -
        (1964*W3)/27 + (448*W4)/(9*Z) + (8*W2*Z)/9))/v4;
 
coB0h = (nL*(-h3/16 + h4/(16*(h - 4*W)) + (7*h2*W)/36 -
        (11*h*W2)/9 + (10*W3)/9))/v4 +
     (Nc*(nQ - 1)*(-h3/16 + h4/(16*(h - 4*W)) + (7*h2*W)/36 -
        (11*h*W2)/9 + (10*W3)/9))/v4;
 
coB0t = (Nc*nL*((16*T2*W)/27 + (16*T*W2)/27 - (80*W3)/27))/v4 +
     (Nc2*(nQ - 1)*((16*T2*W)/27 + (16*T*W2)/27 - (80*W3)/27))/v4 +
     (Nc2*((20*T3)/9 + (76*T2*W)/27 + (40*T*W2)/27 - (80*W3)/27))/v4 +
     (Nc*(-h3/48 + (h2*T)/24 - (17*h*T2)/18 - (5*T3)/3 +
        (h2*(h - 4*T)*(h + 2*T))/(48*(h - 4*W)) - (h2*W)/12 -
        (17*h*T*W)/18 - (5*T2*W)/9 - (7*h*W2)/9 + (826*T*W2)/27 -
        (1960*W3)/27 + (8*T3*W)/(3*Z) - (328*T2*W2)/(27*Z) -
        (368*T*W3)/(9*Z) + (1312*W4)/(27*Z) - (53*T2*Z)/54 - (T*W*Z)/6 +
        (25*W2*Z)/27 - (25*T*Z2)/72 - (5*W*Z2)/72 - (25*Z3)/288 -
        ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(288*(4*T - Z))))/v4;
 
coB0Z = (Nc*(nQ - 1)*((-82*W3)/9 - (208*W4)/(27*Z) + (158*W2*Z)/9 -
        (44*W*Z2)/27 + (8*Z3)/9))/v4 +
     (nL*(-2*W3 - (112*W4)/(9*Z) + (74*W2*Z)/3 - (164*W*Z2)/9 + 8*Z3))/v4;
 
coBht = (Nc*((893*h3)/48 - (455*h2*T)/24 + (83*h*T2)/6 -
       (23*T3)/6 + (h*(35*h3 - 66*h2*T + 96*h*T2 - 32*T3))/
        (48*(h - 4*W)) - h5/(2*(h - 4*W)*(h + T - W)) + (8*h4)/(3*W) -
       (8*h3*T)/(3*W) + (2*h2*T2)/W - (2*h*T3)/(3*W) - (2561*h2*W)/36 +
       (671*h*T*W)/18 - (61*T2*W)/18 + (409*h*W2)/9 + (55*T*W2)/9 +
       (10*W3)/9 - (h*(16*h4 + 97*h3*W - 552*h2*W2 + 732*h*W3 -
          292*W4))/(6*(h + T - W)*W)))/v4;
 
coBhW = (nL*(-h3/72 + h4/(72*(h - 4*W)) - (43*h2*W)/54 +
        (98*h*W2)/27 - (40*W3)/3))/v4 +
     (Nc*(nQ - 1)*(-h3/72 + h4/(72*(h - 4*W)) - (43*h2*W)/54 +
        (98*h*W2)/27 - (40*W3)/3))/v4 +
     (Nc*(-h3/72 - (41*h2*T)/18 + (4*h*T2)/3 + (h3*(h - 12*T))/
         (72*(h - 4*W)) - (43*h2*W)/54 + (70*h*T*W)/9 - (16*T2*W)/3 +
        (98*h*W2)/27 - (16*T*W2)/3 + (16*T2*W2)/h - (40*W3)/3))/v4 +
     ((161*h3)/96 - (4435*h2*W)/216 + (7895*h*W2)/108 - (1082*W3)/3 -
       (40*W4)/h - (h2*(h2 - 4*h*W + 12*W2))/(16*(h - 4*Z)) + h4/(8*Z) -
       (13*h3*W)/(6*Z) + (362*h2*W2)/(9*Z) - (1496*h*W3)/(9*Z) +
       (1216*W4)/(3*Z) + (25*h2*Z)/18 - (44*h*W*Z)/9 + (17*W2*Z)/3 -
       (9*h*Z2)/2 + 14*W*Z2 - (20*W2*Z2)/h -
       (h2*(36*h3 - 271*h2*Z - 24*h*Z2 + 144*Z3))/(288*(h - 4*W)*Z))/v4;
 
coBtZ = (Nc*((-1411*T3)/54 + (313*T2*W)/54 + (1123*T*W2)/27 +
       (190*W3)/9 - (64*W4)/(9*T) - (32*T4)/(9*Z) - (56*T3*W)/(9*Z) -
       (320*T2*W2)/(27*Z) + (440*T*W3)/(27*Z) - (112*W4)/(9*Z) +
       (1331*T2*Z)/54 - (2*T3*Z)/(3*W) + (158*T*W*Z)/9 - (560*W2*Z)/9 +
       (160*W3*Z)/(9*T) - (1061*T*Z2)/27 + (2*T2*Z2)/W - 56*W*Z2 -
       (140*W2*Z2)/(9*T) + (326*Z3)/9 - (8*T*Z3)/(3*W) +
       (44*W*Z3)/(9*T) + (8*Z4)/(3*W) - (2*T2*(4*T - Z)*(2*T + Z)*(2*T + Z))/
        (9*(T - W)*(T - W)) + (2*T*(2*T + Z)*(24*T3 + 122*T2*Z + 5*T*Z2 -
          7*Z3))/(27*(T - W)*Z) + (2*(16*W6 - 80*W5*Z + 473*W4*Z2 -
          215*W3*Z3 - 325*W2*Z4 + 143*W*Z5 + 12*Z6))/
        (9*W*(-T + W - Z)*Z)))/v4;
 
coBWZ = (nL*((-3376*W3)/27 + (640*W5)/(9*Z2) + (2464*W4)/(27*Z) +
        (4280*W2*Z)/27 - (860*W*Z2)/27 - (8*Z3)/3))/v4 +
     (Nc*(nQ - 1)*((-1040*W3)/81 + (3200*W5)/(81*Z2) +
        (22688*W4)/(243*Z) + (4520*W2*Z)/81 - (3260*W*Z2)/243 -
        (88*Z3)/81))/v4 + (Nc*((956*T*W2)/9 + (2908*W3)/81 +
        (4096*T*W5)/(27*Z3) - (10496*T*W4)/(81*Z2) + (896*W5)/(81*Z2) -
        (3536*T*W3)/(27*Z) + (21536*W4)/(243*Z) - (1312*T*W*Z)/81 +
        (2333*W2*Z)/81 - (5*T*Z2)/27 - (2540*W*Z2)/243 - (277*Z3)/324 +
        ((-4*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z)*(12*W2 + 20*W*Z + Z2))/
         (108*(4*T - Z)*Z)))/v4 + ((45*h3)/64 + (5*h2*W)/2 +
       (149*h*W2)/9 + (57992*W3)/27 +
       (h*(h - 16*W)*(h2 + 80*h*W + 192*W2))/(192*(h - 4*Z)) -
       (512*W6)/Z3 - (992*W5)/(3*Z2) + h4/(8*Z) + (h3*W)/(2*Z) +
       (2*h2*W2)/Z + (8*h*W3)/(3*Z) - (1424*W4)/Z - (11*h2*Z)/16 -
       (52*h*W*Z)/9 - (10204*W2*Z)/27 + (11*h*Z2)/36 + (728*W*Z2)/27 +
       (28*Z3)/9 - (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(24*(h - 4*W)*Z))/v4;
 
coI00h = (-3*h*W)/v4;
 
coI00t = (-8*Nc*T*W)/v4;
 
coI00W = (nL*(-h2/12 + h3/(12*(h - 4*W)) - (h*W)/3 + 60*W2 -
        (32*W3)/Z))/v4 + (Nc*(nQ - 1)*(-h2/12 + h3/(12*(h - 4*W)) -
        (h*W)/3 + 60*W2 - (32*W3)/Z))/v4 +
     ((64*h2)/3 - (272*h*W)/3 + (6334*W2)/9 + (6128*W5)/(3*Z3) +
       (1646*W4)/(3*Z2) - (64*h2*W)/(3*Z) + (272*h*W2)/(3*Z) -
       (32572*W3)/(9*Z) + (944*W*Z)/3 + (64*Z2)/3)/v4;
 
coI00Z = (-12*W2 + (18*W3)/Z)/v4 +
     (Nc*(nQ - 1)*((-784*W2)/81 + (320*W4)/(27*Z2) - (64*W3)/(81*Z) +
        (872*W*Z)/81 + (68*Z2)/27))/v4 +
     (nL*((-400*W2)/9 + (64*W4)/(3*Z2) + (64*W3)/(9*Z) + (56*W*Z)/3 +
        12*Z2))/v4 + (Nc*((80*T2)/27 + (8*T*W)/3 - (8*W2)/81 +
        (64*W4)/(27*Z2) + (16*T*W2)/(27*Z) - (128*W3)/(81*Z) +
        (20*T*Z)/9 + (472*W*Z)/81 + (22*Z2)/27 + (2*T2*(2*T + Z)*(2*T + Z))/
         (9*(T - W)*(T - W)) - (2*T*(2*T + Z)*(26*T + 7*Z))/(27*(T - W))))/v4;
 
coI0hW = (6*h*W + 6*W2)/v4;
 
coI0hZ = (3*h*W + 3*W*Z)/v4;
 
coI0tW = (Nc*(-h2/12 - (h*T)/2 - 10*T2 + (h*(h - 4*T)*(h + 2*T))/
        (12*(h - 4*W)) - (6*T3)/(T - W) - (2*h*T2)/(3*W) - (h*W)/3 +
       12*T*W + (148*W2)/3 + (8*W3)/(3*T) - (16*T*W2)/Z - (32*W3)/Z -
       (2*T*Z)/3 - (2*T2*Z)/(3*W)))/v4;
 
coI0WZ = (18*W2 + (6*W4)/Z2 + (18*W3)/Z + 6*W*Z)/v4;
 
coIhhh = ((11*h2)/8 + (9*h3)/(8*(h - 4*W)) - 6*h*W)/v4;
 
coIhtt = (Nc*((5*h*T)/2 - 8*T2 + (h*(h - 4*T)*T)/(2*(h - 4*W)) +
       2*h*W - 10*T*W + (12*T2*W)/h))/v4;
 
coIhWW = ((241*h2)/32 + h3/(3*W) - (247*h*W)/8 + (511*W2)/6 -
      (18*W3)/h - h3/(4*Z) - (7*h2*W)/(3*Z) + (20*h*W2)/(3*Z) -
      (32*W3)/Z - (5*h*Z)/12 + (h2*Z)/(3*W) + W*Z +
      (h2*(24*h2 - 19*h*Z + 8*Z2))/(96*(h - 4*W)*Z))/v4;
 
coIhZZ = ((2*h2)/3 - (439*h*W)/48 + 8*W2 +
      (h*(2*h2 - 41*h*W + 384*W2))/(48*(h - 4*Z)) + (4*h2*W2)/Z2 +
      (4*h2*W)/(3*Z) - (8*h*W2)/Z - (7*h*Z)/3 + (421*W*Z)/12 + (31*Z2)/6 -
      (9*W*Z2)/h + (h*(h2 - 4*h*Z + 12*Z2))/(8*(h - 4*W)))/v4;
 
coIttZ = (Nc*((-1384*T*W)/81 - (158*W2)/81 + (64*W3)/(9*T) +
       (2048*T*W4)/(27*Z3) - (7936*T*W3)/(81*Z2) + (640*W4)/(27*Z2) +
       (2528*T*W2)/(81*Z) - (1952*W3)/(81*Z) + (T*Z)/2 + (109*W*Z)/81 -
       (32*W2*Z)/(3*T) + (65*Z2)/24 + (44*W*Z2)/(9*T) +
       ((-8*W + 5*Z)*(-8*W + 5*Z)*(48*W2 - 40*W*Z + Z2))/(216*(4*T - Z)*Z)))/v4;
 
coIWWZ = ((17*h2)/12 - (23*h*W)/3 - (4552*W2)/9 - (256*W5)/Z3 +
      (304*W4)/(3*Z2) + h3/(4*Z) + (h2*W)/Z + (8*h*W2)/(3*Z) +
      (968*W3)/(9*Z) + 3*h*Z + (770*W*Z)/9 + 11*Z2 + (h*Z2)/(3*W) +
      Z3/(3*W) - (h*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(12*(h - 4*W)*Z))/v4;
 
coM00000 = (Nc*(nQ - 1)*((64*W4)/27 - (64*W5)/(27*Z)))/v4;
 
coM0000Z = (nL*((16*W4)/3 + 8*W3*Z - (8*W*Z3)/3))/v4 +
     (Nc*(nQ - 1)*((16*W4)/3 + (64*W5)/(27*Z) + (88*W3*Z)/27 -
        (8*W*Z3)/27))/v4;
 
coM000W0 = (nL*((64*W4)/3 - (64*W5)/(3*Z)))/v4 +
     (Nc*(nQ - 1)*((64*W4)/3 - (64*W5)/(3*Z)))/v4;
 
coM00tt0 = (Nc*((32*T3*W)/27 - (32*T*W3)/9 + (64*W4)/27 -
       (32*T3*W2)/(27*Z) + (32*T*W4)/(9*Z) - (64*W5)/(27*Z)))/v4;
 
coM00ttZ = (Nc*((-16*T3*W)/27 + (80*T2*W2)/27 -
       (112*T*W3)/27 + (16*W4)/3 + (32*T3*W2)/(27*Z) - (32*T*W4)/(9*Z) +
       (64*W5)/(27*Z) - (16*T3*Z)/27 + (32*T2*W*Z)/27 - (64*T*W2*Z)/27 +
       (88*W3*Z)/27 - (4*T2*Z2)/27 - (16*T*W*Z2)/27 - (8*W*Z3)/27))/v4;
 
coM00tW0 = (Nc*((32*T3*W)/9 - (32*T*W3)/3 + (64*W4)/9 -
       (32*T3*W2)/(9*Z) + (32*T*W4)/(3*Z) - (64*W5)/(9*Z)))/v4;
 
coM0tW0t = (Nc*((-64*T3*W)/9 - (128*T2*W2)/9 + (64*T*W3)/9 +
       (128*W4)/9 + (64*T3*W2)/(9*Z) + (128*T2*W3)/(9*Z) -
       (64*T*W4)/(9*Z) - (128*W5)/(9*Z)))/v4;
 
coM0W0Z0 = (nL*((128*W4)/3 + (64*W5)/(3*Z)))/v4 +
     (Nc*(nQ - 1)*((128*W4)/3 + (64*W5)/(3*Z)))/v4;
 
coM0Wtht = (Nc*((8*h*T3)/3 - 8*h*T2*W + (16*T3*W)/3 +
       (16*h*T*W2)/3 + (16*T2*W2)/3 - (32*T*W3)/3))/v4;
 
coM0WtZt = (Nc*((16*T3*W)/3 + (16*T2*W2)/3 - (320*T*W3)/9 +
       (224*W4)/9 - (64*T3*W2)/(9*Z) - (128*T2*W3)/(9*Z) +
       (64*T*W4)/(9*Z) + (128*W5)/(9*Z) - (8*T3*Z)/9 + (152*T2*W*Z)/9 -
       (80*T*W2*Z)/9 - (64*W3*Z)/9))/v4;
 
coM0WW0W = (-256*W4 - (256*W6)/Z2 + (512*W5)/Z)/v4;
 
coM0WWhW = ((-16*h3*W)/3 + 32*h2*W2 - (320*h*W3)/3 +
      128*W4 + (16*h3*W2)/(3*Z) - (32*h2*W3)/Z + (320*h*W4)/(3*Z) -
      (128*W5)/Z)/v4;
 
coM0WWZW = (16*W*(W - Z)*(2*W - Z)*(4*W - Z)*
      (12*W2 + 20*W*Z + Z2))/(3*v4*Z2);
 
coM0ZtW0 = (Nc*((64*T2*W2)/9 - (128*T*W3)/9 + (160*W4)/9 +
       (32*T3*W2)/(9*Z) - (32*T*W4)/(3*Z) + (64*W5)/(9*Z) - (8*T3*Z)/9 +
       (16*T2*W*Z)/9 + (8*T*W2*Z)/3 + (64*W3*Z)/9 - (8*T2*Z2)/9 +
       (32*T*W*Z2)/9))/v4;
 
coMhhWWh = (-h4 + 12*h2*W2 - 24*h*W3)/v4;
 
coMhWWhW = (-h4/3 - (8*h3*W)/3 + (8*h2*W2)/3 + (32*h*W3)/3 -
      16*W4)/v4;
 
coMhWWZW = ((16*h3*W)/3 - (104*h2*W2)/3 + 32*h*W3 +
      (544*W4)/3 - (16*h3*W2)/(3*Z) + (32*h2*W3)/Z - (320*h*W4)/(3*Z) +
      (128*W5)/Z - (4*h3*Z)/3 + (32*h2*W*Z)/3 - (80*h*W2*Z)/3 -
      (128*W3*Z)/3 - (2*h2*Z2)/3 - (8*h*W*Z2)/3 - (8*W2*Z2)/3)/v4;
 
coMhZWWZ = ((-4*h3*W)/3 + (8*h2*W2)/3 - (368*h*W3)/3 +
      128*W4 + (16*h2*W*Z)/3 + (64*h*W2*Z)/3 + (544*W3*Z)/3 -
      (2*h2*Z2)/3 - (4*h*W*Z2)/3 - (128*W2*Z2)/3 - (4*h*Z3)/3 -
      (8*W*Z3)/3)/v4;
 
coMWWZZh = ((-2*h3*W)/3 + (4*h2*W2)/3 - (184*h*W3)/3 +
      64*W4 + (8*h2*W*Z)/3 + (32*h*W2*Z)/3 + (272*W3*Z)/3 - (h2*Z2)/3 -
      (2*h*W*Z2)/3 - (64*W2*Z2)/3 - (2*h*Z3)/3 - (4*W*Z3)/3)/v4;
 
coMWZZWW = ((2000*W4)/3 - (256*W6)/Z2 - (1408*W5)/(3*Z) -
      416*W3*Z - (8*W2*Z2)/3 + (64*W*Z3)/3 - 3*Z4)/v4;
 
coShhW = ((-7*h2)/4 + (5*h3)/(4*(h - 4*W)) + (2*h3)/(3*W) - h*W -
      (38*W2)/3)/v4;
 
coShWZ = ((-23*h2)/2 + (490*h*W)/9 - (836*W2)/9 - (8*h2*W2)/Z2 -
      h3/(2*Z) - (22*h2*W)/(3*Z) + (376*h*W2)/(9*Z) - (320*W3)/(3*Z) +
      (4*h*Z)/9 - (2*h2*Z)/(3*W) - (536*W*Z)/9 - (35*Z2)/9 -
      (2*h*Z2)/(3*W) + (h*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(6*(h - 4*W)*Z))/v4;
 
coSttW = (Nc*(-(h*T)/3 - (244*T2)/3 - (h2*T)/(h - 4*W) +
       (80*T3)/(T - W) - (604*T*W)/9 - (2596*W2)/81 + (512*W4)/(27*Z2) -
       (1216*W3)/(81*Z) + (16*T*Z)/9 + (116*W*Z)/81 - Z2/3 +
       ((-8*W + 5*Z)*(-8*W + 5*Z)*(12*W2 + 20*W*Z + Z2))/(27*(4*T - Z)*Z)))/v4;
 
coSWZZ = ((-3*h2)/16 + (19*h*W)/12 - (1294*W2)/3 -
      (h2*(h - 12*Z))/(16*(h - 4*W)) - (h2*(h - 12*W))/(16*(h - 4*Z)) -
      (3*h4)/(16*(h - 4*W)*(h - 4*Z)) - (1088*W4)/(3*Z2) -
      (1696*W3)/(3*Z) + (19*h*Z)/12 + (239*W*Z)/3 + 12*Z2 + (2*Z3)/(3*W))/v4;
 
coTh00 = (nL*(-h3/24 + h4/(24*(h - 4*W)) - (h2*W)/6 +
        (2*h*W2)/3 - (8*W3)/3))/v4 +
     (Nc*(nQ - 1)*(-h3/24 + h4/(24*(h - 4*W)) - (h2*W)/6 + (2*h*W2)/3 -
        (8*W3)/3))/v4;
 
coTh0t = (Nc*((223*h3)/24 - (37*h2*T)/4 + (8*h*T2)/3 + (4*T3)/3 +
       (h2*(3*h2 - 6*h*T + 8*T2))/(8*(h - 4*W)) -
       h5/(4*(h - 4*W)*(h + T - W)) + (4*h4)/(3*W) - (4*h3*T)/(3*W) +
       (4*h2*T2)/(3*W) - (209*h2*W)/6 + (41*h*T*W)/3 + (8*T2*W)/3 +
       27*h*W2 - (4*T*W2)/3 - (8*W3)/3 -
       (h*(16*h4 + 97*h3*W - 552*h2*W2 + 732*h*W3 - 292*W4))/
        (12*(h + T - W)*W)))/v4;
 
coTh0W = ((8*h3)/3 - 16*h2*W + (224*h*W2)/3 - (320*W3)/3 -
      (8*h3*W)/(3*Z) + (16*h2*W2)/Z - (224*h*W3)/(3*Z) + (320*W4)/(3*Z))/v4;
 
coThWZ = ((-229*h3)/32 + (460*h2*W)/9 - (1198*h*W2)/9 +
      (368*W3)/3 - (h*(h - 16*W)*(h2 + 80*h*W + 192*W2))/(96*(h - 4*Z)) -
      (4*h3*W2)/Z2 + (16*h2*W3)/Z2 - h4/(2*Z) - (14*h3*W)/(3*Z) +
      (316*h2*W2)/(9*Z) - (1072*h*W3)/(9*Z) + (64*W4)/Z -
      (277*h2*Z)/72 - (h3*Z)/(3*W) - (100*h*W*Z)/9 + 24*W2*Z -
      (7*h*Z2)/6 - (h2*Z2)/W - 12*W*Z2 - (2*Z3)/3 +
      (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(6*(h - 4*W)*Z))/v4;
 
coTt00 = (Nc*((-8*T3)/3 - (496*T2*W)/27 - (320*T*W2)/27 -
       (256*W3)/27 + (8*T3*W)/(3*Z) + (496*T2*W2)/(27*Z) +
       (320*T*W3)/(27*Z) + (256*W4)/(27*Z)))/v4;
 
coTt0h = (Nc*((28*h3)/3 - (125*h2*T)/12 + (49*h*T2)/6 -
       (10*T3)/3 + (h*(4*h3 - 7*h2*T + 14*h*T2 - 8*T3))/(12*(h - 4*W)) -
       h5/(4*(h - 4*W)*(h + T - W)) + (4*h4)/(3*W) - (4*h3*T)/(3*W) +
       (2*h2*T2)/W - (2*h*T3)/(3*W) - (110*h2*W)/3 + (70*h*T*W)/3 -
       7*T2*W + (73*h*W2)/3 + 2*T*W2 -
       (h*(16*h4 + 97*h3*W - 552*h2*W2 + 732*h*W3 - 292*W4))/
        (12*(h + T - W)*W)))/v4;
 
coTt0Z = (Nc*((-502*T3)/27 - (833*T2*W)/27 + (350*T*W2)/27 +
       (136*W3)/9 - (16*T4)/(9*Z) - (40*T3*W)/(9*Z) -
       (832*T2*W2)/(27*Z) - (272*T*W3)/(9*Z) - (16*W4)/(9*Z) +
       (1193*T2*Z)/54 - (2*T3*Z)/(3*W) + (391*T*W*Z)/54 - (388*W2*Z)/9 -
       (4595*T*Z2)/216 + (2*T2*Z2)/W - (1717*W*Z2)/72 + (1645*Z3)/96 -
       (4*T*Z3)/(3*W) + (4*Z4)/(3*W) - (4*T3*(2*T + Z)*(2*T + Z))/(9*(T - W)*(T - W)) -
       ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(288*(4*T - Z)) +
       (2*T2*(2*T + Z)*(12*T2 + 64*T*Z + 17*Z2))/(27*(T - W)*Z) +
       (16*W6 - 80*W5*Z + 473*W4*Z2 - 215*W3*Z3 - 325*W2*Z4 +
         143*W*Z5 + 12*Z6)/(9*W*(-T + W - Z)*Z)))/v4;
 
coTW00 = ((-8*h*W2)/3 + 216*W3 - (8*W6)/(3*Z3) +
       (640*W5)/(3*Z2) + (8*h*W3)/(3*Z) - (424*W4)/Z - (8*W2*Z)/3)/
      v4 + (Nc*((-50*T2*W)/9 - 4*T*W2 + (592*W3)/81 + (16*W4)/(3*T) -
        (128*W5)/(27*Z2) + (32*T2*W2)/(9*Z) - (16*T*W3)/(9*Z) +
        (544*W4)/(81*Z) + (4*T*W*Z)/9 + (352*W2*Z)/81 + (10*W*Z2)/27))/
      v4 + (Nc*(nQ - 1)*((272*W3)/81 - (640*W5)/(27*Z2) +
        (6176*W4)/(81*Z) + (1496*W2*Z)/81 + (44*W*Z2)/27))/v4 +
     (nL*((-688*W3)/9 - (128*W5)/(3*Z2) + (1312*W4)/(9*Z) +
        (136*W2*Z)/3 + 4*W*Z2))/v4;
 
coTZ00 = (Nc*(nQ - 1)*((560*W3)/27 + (176*W4)/(9*Z) +
        (236*W2*Z)/27 - (16*W*Z2)/27 + (8*Z3)/9))/v4 +
     (nL*(16*W3 + (80*W4)/(3*Z) + 4*W2*Z - (16*W*Z2)/3 + 8*Z3))/v4;
 
coTZ0t = (Nc*((-232*T3)/27 - (416*T2*W)/27 + (136*T*W2)/9 +
       (64*W3)/3 + (64*W4)/(9*T) - (16*T4)/(9*Z) - (16*T3*W)/(9*Z) -
       (280*T2*W2)/(27*Z) - (280*T*W3)/(27*Z) + (416*W4)/(27*Z) +
       (44*T2*Z)/3 + (197*T*W*Z)/27 - (233*W2*Z)/9 - (32*W3*Z)/(9*T) -
       (719*T*Z2)/36 + (4*T2*Z2)/(3*W) - (3425*W*Z2)/108 -
       (52*W2*Z2)/(9*T) + (2761*Z3)/144 - (4*T*Z3)/(3*W) +
       (44*W*Z3)/(9*T) + (4*Z4)/(3*W) - (2*T2*(2*T - Z)*(2*T + Z)*(2*T + Z))/
        (9*(T - W)*(T - W)) + ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/
        (144*(4*T - Z)) + (2*T*(2*T + Z)*(12*T3 + 64*T2*Z - 9*T*Z2 -
          7*Z3))/(27*(T - W)*Z) + (16*W6 - 80*W5*Z + 473*W4*Z2 -
         215*W3*Z3 - 325*W2*Z4 + 143*W*Z5 + 12*Z6)/
        (9*W*(-T + W - Z)*Z)))/v4;
 
coTZ0W = (832*W3 + (640*W5)/(3*Z2) - (2240*W4)/(3*Z) -
      (1040*W2*Z)/3 + (136*W*Z2)/3 + (8*Z3)/3)/v4;
 
coTZhW = ((-49*h3)/192 + (28*h2*W)/3 - 13*h*W2 + (368*W3)/9 +
      (h*(h - 16*W)*(h2 + 80*h*W + 192*W2))/(192*(h - 4*Z)) +
      (4*h2*W2)/(3*Z) - (112*h*W3)/(3*Z) + (64*W4)/Z - (199*h2*Z)/16 +
      45*h*W*Z - 72*W2*Z + (19*h*Z2)/12 - (h2*Z2)/W - (92*W*Z2)/3 -
      (20*Z3)/9 - (h*Z3)/(3*W) + (h*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/
       (12*(h - 4*W)))/v4;
 
coTbar0hW = ((16*h2*W)/3 - (64*h*W2)/3 + 64*W3 -
      (16*h2*W2)/(3*Z) + (64*h*W3)/(3*Z) - (64*W4)/Z)/v4;
 
coTbar0WZ = (-448*W3 + (256*W5)/Z2 + (320*W4)/(3*Z) + 80*W2*Z +
      (16*W*Z2)/3)/v4;
 
coU0t0W = (Nc*(-2*T3 - (20*T2*W)/9 + (62*T*W2)/9 - (32*W3)/3 +
       (8*W4)/(3*T) - (16*T2*W2)/(3*Z) - (32*T*W3)/(3*Z) + (8*T2*Z)/9 -
       (32*T*W*Z)/9))/v4;
 
coU0tht = (Nc*((5*h*T2)/3 - (10*T3)/3 - h*T*W - 6*T2*W +
       2*h*W2 + (4*T*W2)/3))/v4;
 
coU0ttZ = (Nc*((-176*T2*W)/27 - (496*T*W2)/27 - (88*W3)/27 +
       (64*W4)/(9*T) + (160*T2*W2)/(27*Z) + (32*T*W3)/(3*Z) -
       (64*W4)/(27*Z) + (77*T2*Z)/54 + (241*T*W*Z)/54 - (61*W2*Z)/27 -
       (32*W3*Z)/(3*T) + (331*T*Z2)/216 + (23*W*Z2)/24 +
       (44*W2*Z2)/(9*T) + (25*Z3)/288 +
       ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(288*(4*T - Z))))/v4;
 
coUhW00 = (nL*(h3/48 - h4/(48*(h - 4*W)) + (h2*W)/12 - h*W2 +
        (20*W3)/3))/v4 + (Nc*(nQ - 1)*(h3/48 - h4/(48*(h - 4*W)) +
        (h2*W)/12 - h*W2 + (20*W3)/3))/v4;
 
coUhW0t = (Nc*(h3/48 + (5*h2*T)/8 - (13*h*T2)/6 -
       (h2*(h - 4*T)*(h + 2*T))/(48*(h - 4*W)) + (2*h2*T2)/(3*W) +
       (h2*W)/12 - (5*h*T*W)/6 - (2*T2*W)/3 - h*W2 + 2*T*W2 +
       (20*W3)/3))/v4;
 
coUhWhW = ((17*h3)/8 - h4/(8*(h - 4*W)) - h4/(3*W) -
      (29*h2*W)/6 - (2*h*W2)/3 + 8*W3)/v4;
 
coUhWWZ = ((-17*h3)/48 + (53*h2*W)/4 - (109*h*W2)/3 -
      (148*W3)/3 - h4/(16*Z) - (h3*W)/(4*Z) - (7*h2*W2)/(3*Z) +
      (20*h*W3)/(3*Z) - (32*W4)/Z - 4*h2*Z + (32*h*W*Z)/3 + 4*W2*Z +
      (13*h*Z2)/12 - (h2*Z2)/(3*W) + (W*Z2)/3 +
      (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(48*(h - 4*W)*Z))/v4;
 
coUW0tt = (Nc*((-64*T2*W)/9 - (9280*T*W2)/81 - (128*W3)/9 -
       (4096*T*W5)/(27*Z3) + (12032*T*W4)/(81*Z2) - (256*W5)/(9*Z2) +
       (64*T2*W2)/(9*Z) + (9536*T*W3)/(81*Z) + (128*W4)/(3*Z)))/v4;
 
coUWhhh = ((15*h3)/16 + (9*h4)/(16*(h - 4*W)) - (37*h2*W)/4 +
      9*h*W2)/v4;
 
coUWhtt = (Nc*((13*h2*T)/12 - 3*h*T2 + (h2*(h - 4*T)*T)/
        (4*(h - 4*W)) - 5*h*T*W + 20*T2*W + (4*T*W2)/3))/v4;
 
coUWhWW = ((71*h3)/48 + (3*h4)/(16*(h - 4*W)) - (25*h2*W)/12 +
      21*h*W2 - 12*W3 + (2*h2*Z)/3 + (8*h*W*Z)/3 + (8*W2*Z)/3)/v4;
 
coUWhZZ = ((3*h3)/8 - (h2*W)/3 + (179*h*W2)/12 + 64*W3 +
      (h2*(h2 - 4*h*W + 12*W2))/(16*(h - 4*Z)) - h2*Z - 6*h*W*Z -
      (97*W2*Z)/3 + (35*h*Z2)/12 - (41*W*Z2)/3 +
      (h2*(h2 - 4*h*Z + 12*Z2))/(16*(h - 4*W)))/v4;
 
coUWZ00 = (Nc*((-16*T*W2)/3 - (80*W3)/27 + (128*W5)/(27*Z2) -
        (32*T2*W2)/(9*Z) - (32*T*W3)/(9*Z) - (1216*W4)/(81*Z) +
        (8*T2*Z)/9 - (8*T*W*Z)/9 - (344*W2*Z)/27 + (4*T*Z2)/9 +
        (292*W*Z2)/81 + (10*Z3)/27))/v4 +
     (Nc*(nQ - 1)*((1120*W3)/27 + (640*W5)/(27*Z2) - (6080*W4)/(81*Z) -
        (1376*W2*Z)/27 + (1232*W*Z2)/81 + (44*Z3)/27))/v4 +
     (nL*((1184*W3)/9 + (128*W5)/(3*Z2) - (1216*W4)/(9*Z) -
        (1120*W2*Z)/9 + (112*W*Z2)/3 + 4*Z3))/v4;
 
coUWZhZ = (h3/192 - (5*h2*W)/3 - (13*h*W2)/3 - 48*W3 -
      (h*(h - 16*W)*(h2 + 80*h*W + 192*W2))/(192*(h - 4*Z)) +
      (8*h2*W3)/Z2 + (8*h2*W2)/(3*Z) - (16*h*W3)/Z + (h2*Z)/48 +
      (28*h*W*Z)/3 - 60*W2*Z - (35*h*Z2)/12 + (88*W*Z2)/3 + (5*Z3)/3)/v4;
 
coUWZtt = (Nc*((16*T2*W)/3 + (1180*T*W2)/9 - (116*W3)/27 +
       (4096*T*W5)/(27*Z3) - (10496*T*W4)/(81*Z2) + (1280*W5)/(27*Z2) -
       (64*T2*W2)/(9*Z) - (3584*T*W3)/(27*Z) - (4480*W4)/(81*Z) -
       (8*T2*Z)/9 - (700*T*W*Z)/81 - (101*W2*Z)/9 - (41*T*Z2)/27 +
       (700*W*Z2)/81 + (37*Z3)/36 - ((-4*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z)*
         (12*W2 + 20*W*Z + Z2))/(108*(4*T - Z)*Z)))/v4;
 
coUWZWW = ((-16*h2*W)/3 + (64*h*W2)/3 + 320*W3 - (512*W6)/Z3 -
      (832*W5)/(3*Z2) + (16*h2*W2)/(3*Z) - (64*h*W3)/(3*Z) +
      (7696*W4)/(9*Z) + (4*h2*Z)/3 - (16*h*W*Z)/3 + (160*W2*Z)/3 -
      (76*W*Z2)/9 + (13*Z3)/3)/v4;
 
coUZW00 = (nL*((-64*W3)/3 - (80*W4)/(3*Z) + 4*W2*Z))/v4 +
     (Nc*(nQ - 1)*((-64*W3)/3 - (80*W4)/(3*Z) + 4*W2*Z))/v4;
 
coUZW0t = (Nc*((-88*T2*W)/3 - (224*T*W2)/9 - (64*W3)/3 +
       (8*T2*W2)/(3*Z) - (8*T*W3)/Z - (80*W4)/(3*Z) + (26*T2*Z)/3 +
       (2*T*W*Z)/9 + 4*W2*Z + (2*T*Z2)/3 + (2*T2*Z2)/(3*W)))/v4;
 
coUZWhW = ((-17*h3)/48 + (53*h2*W)/4 - (109*h*W2)/3 -
      (148*W3)/3 - h4/(16*Z) - (h3*W)/(4*Z) - (7*h2*W2)/(3*Z) +
      (20*h*W3)/(3*Z) - (32*W4)/Z - 4*h2*Z + (32*h*W*Z)/3 + 4*W2*Z +
      (13*h*Z2)/12 - (h2*Z2)/(3*W) + (W*Z2)/3 +
      (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(48*(h - 4*W)*Z))/v4;
 
coUZWWZ = ((4*h2*W)/3 - (8*h*W2)/3 - (1024*W3)/3 +
      (128*W5)/Z2 + (1184*W4)/(3*Z) - (8*h*W*Z)/3 + (832*W2*Z)/3 -
      (44*W*Z2)/3 - (26*Z3)/3 - Z4/(3*W))/v4;
 
coAhAh = ((-5*h)/4 + (23*h2)/(12*(h - 4*W)) + (4*h2)/(9*W) - 5*W -
      (8*W2)/h)/v4;
 
coAhAt = (Nc*((92*h)/9 - (16*T)/3 + (h*(h - 2*T))/(3*(h - 4*W)) -
       h3/(3*(h - 4*W)*(h + T - W)) + (4*h2)/(3*W) - (4*h*T)/(9*W) -
       (59*W)/3 - (4*T*W)/h - (4*h3 + 28*h2*W - 110*h*W2 + 73*W3)/
        (3*(h + T - W)*W)))/v4;
 
coAhAW = ((325*h)/24 + (19*h2)/(9*W) - (625*W)/18 +
      (290*W2)/(3*h) - (8*h*W2)/Z2 - (2*h2)/Z - (8*h*W)/Z +
      (160*W2)/(9*Z) - (320*W3)/(3*h*Z) + (13*Z)/9 + (7*h*Z)/(9*W) +
      (h*(48*h2 - 13*h*Z + 16*Z2))/(24*(h - 4*W)*Z))/v4;
 
coAhAZ = ((-119*h)/36 - (41*W)/36 + (2*h2 - 41*h*W + 384*W2)/
       (12*(h - 4*Z)) + (4*h*W2)/Z2 - h2/(2*Z) + (2*h*W)/(3*Z) +
      (40*W2)/(9*Z) - (14*Z)/9 - (4*h*Z)/(9*W) + (3*W*Z)/h +
      (h*(3*h - 4*Z)*(2*h + Z))/(12*(h - 4*W)*Z))/v4;
 
coAhB00 = (nL*(h2/24 - h3/(24*(h - 4*W)) - (5*h*W)/18 +
        (2*W2)/3 - (8*W3)/(3*h)))/v4 +
     (Nc*(nQ - 1)*(h2/24 - h3/(24*(h - 4*W)) - (5*h*W)/18 + (2*W2)/3 -
        (8*W3)/(3*h)))/v4;
 
coAhB0t = (Nc*(h2/24 + (29*h*T)/36 - (2*T2)/3 + (4*T3)/(3*h) -
       (h*(h - 4*T)*(h + 2*T))/(24*(h - 4*W)) + (8*h*T2)/(9*W) -
       (5*h*W)/18 - 2*T*W + (8*T2*W)/(3*h) + (8*W2)/3 - (4*T*W2)/(3*h) -
       (8*W3)/(3*h)))/v4;
 
coAhBhW = ((149*h2)/72 - (7*h3)/(24*(h - 4*W)) - (7*h3)/(9*W) -
      (29*h*W)/6 + (14*W2)/3 - (16*W3)/h)/v4;
 
coAhBWZ = ((-45*h2)/32 + (245*h*W)/9 - (214*W2)/3 +
      (368*W3)/(3*h) - ((h - 16*W)*(h2 + 80*h*W + 192*W2))/
       (96*(h - 4*Z)) + (8*h*W3)/Z2 - h3/(4*Z) - (h2*W)/Z +
      (20*h*W2)/(3*Z) - (112*W3)/(3*Z) + (64*W4)/(h*Z) - (557*h*Z)/72 +
      (56*W*Z)/3 + (24*W2*Z)/h + Z2/2 - (7*h*Z2)/(9*W) - (12*W*Z2)/h -
      (2*Z3)/(3*h) + (h*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(12*(h - 4*W)*Z))/v4;
 
coAtAt = (Nc2*((8*T)/9 - (8*T2)/(9*W) + (16*W)/9))/v4 +
     (Nc*((-158*T)/3 - (2*h*T)/(h - 4*W) + (46*T2)/(T - W) - (4366*W)/81 +
        (12*T*W)/h + (32*W2)/(3*T) + (2048*W4)/(27*Z3) -
        (7936*W3)/(81*Z2) + (128*W4)/(9*T*Z2) + (3008*W2)/(81*Z) -
        (64*W3)/(3*T*Z) - (47*Z)/18 - (44*W*Z)/(9*T) +
        ((40*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(18*(4*T - Z)*Z)))/v4;
 
coAtAW = (Nc*((-2*h)/9 - (254*T)/3 + (2*h*(h - 8*T))/(3*(h - 4*W)) +
       (68*T2)/(T - W) - (14*h*T)/(9*W) + (8056*W)/81 + (16*T*W)/h +
       (8*W2)/(3*T) + (4096*W4)/(27*Z3) - (12032*W3)/(81*Z2) -
       (64*T*W)/(3*Z) - (15440*W2)/(81*Z) - (8*Z)/9 - (14*T*Z)/(9*W) +
       (4*(-8*W + 5*Z)*(-8*W + 5*Z)*(12*W2 + 20*W*Z + Z2))/(27*(4*T - Z)*Z2)))/v4;
 
coAtAZ = (Nc*((-74*T)/27 - (1369*W)/81 + (185*W2)/(9*T) +
       (2048*W4)/(27*Z3) + (16*T3)/(9*Z2) + (16*T2*W)/(9*Z2) +
       (16*T*W2)/(9*Z2) - (7792*W3)/(81*Z2) + (16*W4)/(9*T*Z2) +
       (232*T2)/(27*Z) - (32*T*W)/(27*Z) - (760*W2)/(81*Z) + (1411*Z)/54 -
       (4*T*Z)/(9*W) + (302*W*Z)/(9*T) - (155*Z2)/(9*T) + (4*Z2)/(3*W) -
       (4*Z3)/(3*T*W) + (4*T2*(2*T + Z)*(2*T + Z))/(9*(T - W)*(T - W)*Z) +
       ((-8*W + 5*Z)*(-8*W + 5*Z)*(48*W2 - 40*W*Z + Z2))/(54*(4*T - Z)*Z2) -
       (2*T*(2*T + Z)*(12*T2 + 64*T*Z + 17*Z2))/(27*(T - W)*Z2) -
       (16*W6 - 80*W5*Z + 473*W4*Z2 - 215*W3*Z3 - 325*W2*Z4 +
         143*W*Z5 + 12*Z6)/(9*T*W*(-T + W - Z)*Z2)))/v4;
 
coAtB00 = (Nc*nL*((8*T*W)/9 + (32*W2)/9))/v4 +
     (Nc2*(nQ - 1)*((8*T*W)/9 + (32*W2)/9))/v4;
 
coAtB0t = (Nc2*((-8*T2)/9 - (16*T3)/(9*W) + (16*T*W)/9 +
        (32*W2)/9))/v4 + (Nc*((-20*T2)/3 + (6*T3)/(T - W) -
        (508*T*W)/27 - (1240*W2)/27 + (32*W3)/(27*T) + (320*T*W2)/(27*Z) +
        (256*W3)/(9*Z) + (64*W4)/(27*T*Z) - (155*T*Z)/54 - (79*W*Z)/54 -
        (44*W2*Z)/(9*T) - (25*Z2)/72 - ((-4*W + Z)*(8*W + Z)*
          (-8*W + 5*Z)*(-8*W + 5*Z))/(72*(4*T - Z)*Z)))/v4;
 
coAtBhW = (Nc*((-5*h2)/18 - (11*h*T)/9 - (h2*(h - 2*T))/
        (6*(h - 4*W)) + (14*h2*T)/(9*W) - (2*h*W)/9 + (4*T*W)/3 + 8*W2 +
       (16*T*W2)/h))/v4;
 
coAtBWZ = (Nc*((-616*T*W)/9 + (140*W2)/9 + (4096*W5)/(27*Z3) -
       (10496*W4)/(81*Z2) - (16*T*W2)/Z - (5264*W3)/(27*Z) +
       (188*T*Z)/9 + (416*W*Z)/81 + (31*Z2)/27 + (14*T*Z2)/(9*W) -
       ((-4*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z)*(12*W2 + 20*W*Z + Z2))/
        (27*(4*T - Z)*Z2)))/v4;
 
coAWAW = ((-1915*h)/24 + (133*h2)/(9*W) + (9149*W)/18 - (42*W2)/h +
      (1532*W4)/(3*Z3) + (1616*W3)/(3*Z2) - (31*h2)/(3*Z) +
      (172*h*W)/(3*Z) - (21868*W2)/(9*Z) + 201*Z - (10*h*Z)/(9*W) +
      (133*Z2)/(9*W) - (h*(24*h2 - 73*h*Z + 8*Z2))/(8*(h - 4*W)*Z))/v4;
 
coAWAZ = ((431*h)/36 + (1961*W)/9 - (h*(h - 12*W))/(4*(h - 4*Z)) -
      (3*h3)/(4*(h - 4*W)*(h - 4*Z)) - (2176*W4)/(3*Z3) + (8*h*W2)/Z2 -
      (302*W3)/(3*Z2) + h2/Z + (20*h*W)/(3*Z) - (1564*W2)/Z + (139*Z)/3 +
      (7*h*Z)/(9*W) - (12*W*Z)/h + (19*Z2)/(9*W) -
      (h*(12*h2 - 7*h*Z - 44*Z2))/(12*(h - 4*W)*Z))/v4;
 
coAWB00 = (nL*(-h2/12 + h3/(12*(h - 4*W)) + (h*W)/9 + 44*W2 -
        (32*W3)/(3*Z) + (4*W*Z)/9))/v4 +
     (Nc*(nQ - 1)*(-h2/12 + h3/(12*(h - 4*W)) + (h*W)/9 + 44*W2 -
        (32*W3)/(3*Z) + (4*W*Z)/9))/v4;
 
coAWB0t = (Nc*(-h2/12 - (13*h*T)/18 - 2*T2 +
       (h*(h - 4*T)*(h + 2*T))/(12*(h - 4*W)) - (6*T3)/(T - W) -
       (8*h*T2)/(9*W) + (h*W)/9 + 20*T*W + (100*W2)/3 + (8*W3)/(3*T) -
       (32*T2*W)/(3*Z) - (80*T*W2)/(3*Z) - (32*W3)/(3*Z) - (8*T*Z)/9 -
       (8*T2*Z)/(9*W) + (4*W*Z)/9))/v4;
 
coAWBhW = ((751*h2)/72 + (7*h3)/(9*W) - (235*h*W)/6 +
      (322*W2)/3 - (24*W3)/h - h3/(2*Z) + (2*h2*W)/(3*Z) - (8*h*W2)/Z -
      (23*h*Z)/18 + (7*h2*Z)/(9*W) + (10*W*Z)/3 +
      (h2*(12*h2 - 5*h*Z + 4*Z2))/(24*(h - 4*W)*Z))/v4;
 
coAWBWZ = ((17*h2)/6 - (206*h*W)/9 - (5936*W2)/9 -
      (512*W5)/Z3 - (160*W4)/(3*Z2) + h3/(2*Z) + (2*h2*W)/Z +
      (400*W3)/(9*Z) + (70*h*Z)/9 + (944*W*Z)/9 + (178*Z2)/9 +
      (7*h*Z2)/(9*W) + (7*Z3)/(9*W) - (h*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/
       (6*(h - 4*W)*Z))/v4;
 
coAZAZ = ((11*h)/24 + (1933*W)/36 - (h*(h - 12*Z))/(8*(h - 4*W)) -
      (3*h3)/(8*(h - 4*W)*(h - 4*Z)) - (7*h2 - 118*h*W + 768*W2)/
       (24*(h - 4*Z)) - (128*W4)/Z3 - (4*h*W2)/Z2 - (448*W3)/Z2 -
      (4*h*W)/(3*Z) - (1088*W2)/(3*Z) + (185*Z)/18 - (9*W*Z)/h +
      (4*Z2)/(9*W))/v4;
 
coAZB00 = (Nc*(nQ - 1)*((-212*W2)/27 + (176*W4)/(9*Z2) +
        (256*W3)/(27*Z) + (4*W*Z)/27))/v4 +
     (nL*((-92*W2)/9 + (80*W4)/(3*Z2) + (44*W*Z)/9))/v4;
 
coAZB0t = (Nc*((242*T2)/27 - (14*T*W)/27 - (80*W2)/9 -
       (32*W3)/(3*T) - (232*T2*W2)/(27*Z2) - (232*T*W3)/(27*Z2) +
       (464*W4)/(27*Z2) - (448*T2*W)/(27*Z) - (272*T*W2)/(27*Z) +
       (16*W3)/(3*Z) + (64*W4)/(9*T*Z) + (299*T*Z)/108 + (8*T2*Z)/(9*W) +
       (127*W*Z)/108 + (44*W2*Z)/(9*T) + (25*Z2)/144 +
       ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(144*(4*T - Z)*Z)))/v4;
 
coAZBhW = ((-475*h2)/72 + (277*h*W)/18 - (139*W2)/3 +
      (h*(h2 - 4*h*W + 12*W2))/(4*(h - 4*Z)) + (16*h2*W2)/(3*Z2) -
      (64*h*W3)/(3*Z2) + (64*W4)/Z2 - h3/(4*Z) + (15*h2*W)/Z -
      (124*h*W2)/(3*Z) + (272*W3)/(3*Z) + (5*h*Z)/18 - (7*h2*Z)/(9*W) +
      (2*W*Z)/3 - (12*W2*Z)/h + (h2*(3*h - 4*Z)*(2*h + Z))/
       (24*(h - 4*W)*Z))/v4;
 
coAZBWZ = (-h2/48 + (2*h*W)/3 + (5800*W2)/9 +
      ((h - 16*W)*(h2 + 80*h*W + 192*W2))/(48*(h - 4*Z)) - (256*W5)/Z3 -
      (8*h*W3)/Z2 - (1792*W4)/(3*Z2) - (8*h*W2)/(3*Z) -
      (1024*W3)/(3*Z) - (17*h*Z)/12 - (208*W*Z)/3 - (50*Z2)/3 -
      (7*Z3)/(9*W))/v4;
 
coB00B00 = (16*nL2*W3)/(9*v4) + (32*Nc*nL*(nQ - 1)*W3)/
      (9*v4) + (16*Nc2*(nQ - 1)*(nQ - 1)*W3)/(9*v4) +
     (Nc*(nQ - 1)*((-376*W3)/27 + (64*W4)/(9*Z) - (4*W2*Z)/27 +
        (8*W*Z2)/27))/v4 + (nL*((-56*W3)/3 + (32*W4)/(3*Z) -
        (4*W2*Z)/3 + (8*W*Z2)/3))/v4;
 
coB00B0t = (Nc*nL*((8*T2*W)/9 + (8*T*W2)/9 + (32*W3)/9))/v4 + 
      (Nc2*(nQ - 1)*((8*T2*W)/9 + (8*T*W2)/9 + (32*W3)/9))/v4;
 
coB00BhW = (nL*(h3/48 - h4/(48*(h - 4*W)) - (13*h2*W)/36 +
        (7*h*W2)/9 + (4*W3)/3))/v4 +
     (Nc*(nQ - 1)*(h3/48 - h4/(48*(h - 4*W)) - (13*h2*W)/36 +
        (7*h*W2)/9 + (4*W3)/3))/v4;
 
coB00BWZ = (nL*((80*W3)/9 - (16*W4)/(3*Z) - (28*W2*Z)/9 -
        (4*W*Z2)/9))/v4 + (Nc*(nQ - 1)*((80*W3)/9 - (16*W4)/(3*Z) -
        (28*W2*Z)/9 - (4*W*Z2)/9))/v4;
 
coB0tB0t = (Nc2*((-16*T3)/9 - (8*T4)/(9*W) + (8*T*W2)/9 +
        (16*W3)/9))/v4 + (Nc*((-112*T*W2)/27 - (56*W3)/3 +
        (160*T*W3)/(27*Z) + (320*W4)/(27*Z) + (8*T*W*Z)/9 - (4*W2*Z)/27 +
        (8*W*Z2)/27))/v4;
 
coB0tBhW = (Nc*(h3/48 + (61*h2*T)/72 - (55*h*T2)/18 -
       (h2*(h - 4*T)*(h + 2*T))/(48*(h - 4*W)) + (8*h2*T2)/(9*W) -
       (13*h2*W)/36 - (31*h*T*W)/18 + 2*T2*W + (7*h*W2)/9 + (14*T*W2)/3 +
       (4*W3)/3))/v4;
 
coB0tBWZ = (Nc*((-400*T2*W)/9 - 40*T*W2 + (80*W3)/9 -
       (8*T2*W2)/Z - (56*T*W3)/(3*Z) - (16*W4)/(3*Z) + (110*T2*Z)/9 +
       (34*T*W*Z)/9 - (28*W2*Z)/9 + (8*T*Z2)/9 + (8*T2*Z2)/(9*W) -
       (4*W*Z2)/9))/v4;
 
coBhWBhW = ((217*h3)/144 - h4/(16*(h - 4*W)) - (2*h4)/(9*W) -
      (167*h2*W)/36 + 5*h*W2 - 4*W3)/v4;

coBhWBWZ = ((-17*h3)/48 + (749*h2*W)/36 - (599*h*W2)/9 +
      (124*W3)/3 - h4/(16*Z) - (h3*W)/(4*Z) + (3*h2*W2)/Z -
      (44*h*W3)/(3*Z) + (32*W4)/Z - (52*h2*Z)/9 + (160*h*W*Z)/9 -
      (52*W2*Z)/3 + (55*h*Z2)/36 - (4*h2*Z2)/(9*W) - W*Z2 +
      (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(48*(h - 4*W)*Z))/v4;
 
coBWZBWZ = ((2*h2*W)/3 - (4*h*W2)/3 - (3080*W3)/9 -
      (64*W5)/Z2 - (496*W4)/(3*Z) - (4*h*W*Z)/3 + (2384*W2*Z)/9 -
      14*W*Z2 - (55*Z3)/9 - (2*Z4)/(9*W))/v4;
 
result = 
co1 + 
Ah*coAh + 
At*coAt + 
AW*coAW + 
AZ*coAZ + 
B00*coB00 + 
B0h*coB0h + 
B0t*coB0t + 
B0Z*coB0Z + 
Bht*coBht + 
BhW*coBhW + 
BtZ*coBtZ + 
BWZ*coBWZ + 
Ah*Ah*coAhAh + 
Ah*At*coAhAt + 
At*At*coAtAt + 
Ah*AW*coAhAW + 
At*AW*coAtAW + 
AW*AW*coAWAW + 
Ah*AZ*coAhAZ + 
At*AZ*coAtAZ + 
AW*AZ*coAWAZ + 
AZ*AZ*coAZAZ + 
Ah*B00*coAhB00 + 
At*B00*coAtB00 + 
AW*B00*coAWB00 + 
AZ*B00*coAZB00 + 
B00*B00*coB00B00 + 
Ah*B0t*coAhB0t + 
At*B0t*coAtB0t + 
AW*B0t*coAWB0t + 
AZ*B0t*coAZB0t + 
B00*B0t*coB00B0t + 
B0t*B0t*coB0tB0t + 
Ah*BhW*coAhBhW + 
At*BhW*coAtBhW + 
AW*BhW*coAWBhW + 
AZ*BhW*coAZBhW + 
B00*BhW*coB00BhW + 
B0t*BhW*coB0tBhW + 
BhW*BhW*coBhWBhW + 
Ah*BWZ*coAhBWZ + 
At*BWZ*coAtBWZ + 
AW*BWZ*coAWBWZ + 
AZ*BWZ*coAZBWZ + 
B00*BWZ*coB00BWZ + 
B0t*BWZ*coB0tBWZ + 
BhW*BWZ*coBhWBWZ + 
BWZ*BWZ*coBWZBWZ + 
coI00h*I00h + 
coI00t*I00t + 
coI00W*I00W + 
coI00Z*I00Z + 
coI0hW*I0hW + 
coI0hZ*I0hZ + 
coI0tW*I0tW + 
coI0WZ*I0WZ + 
coIhhh*Ihhh + 
coIhtt*Ihtt + 
coIhWW*IhWW + 
coIhZZ*IhZZ + 
coIttZ*IttZ + 
coIWWZ*IWWZ + 
coM00000*M00000W + 
coM0000Z*M0000Z + 
coM000W0*M000W0 + 
coM00tt0*M00tt0 + 
coM00ttZ*M00ttZ + 
coM00tW0*M00tW0 + 
coM0tW0t*M0tW0t + 
coM0W0Z0*M0W0Z0 + 
coM0Wtht*M0Wtht + 
coM0WtZt*M0WtZt + 
coM0WW0W*M0WW0W + 
coM0WWhW*M0WWhW + 
coM0WWZW*M0WWZW + 
coM0ZtW0*M0ZtW0 + 
coMhhWWh*MhhWWh + 
coMhWWhW*MhWWhW + 
coMhWWZW*MhWWZW + 
coMhZWWZ*MhZWWZ + 
coMWWZZh*MWWZZh + 
coMWZZWW*MWZZWW + 
coShhW*ShhW + 
coShWZ*ShWZ + 
coSttW*SttW + 
coSWZZ*SWZZ + 
coTh00*Th00 + 
coTh0t*Th0t + 
coTh0W*Th0W + 
coThWZ*ThWZ + 
coTt00*Tt00 + 
coTt0h*Tt0h + 
coTt0Z*Tt0Z + 
coTW00*TW00 + 
coTZ00*TZ00 + 
coTZ0t*TZ0t + 
coTZ0W*TZ0W + 
coTZhW*TZhW + 
coTbar0hW*Tbar0hW + 
coTbar0WZ*Tbar0WZ + 
coU0t0W*U0t0W + 
coU0tht*U0tht + 
coU0ttZ*U0ttZ + 
coUhW00*UhW00 + 
coUhW0t*UhW0t + 
coUhWhW*UhWhW + 
coUhWWZ*UhWWZ + 
coUW0tt*UW0tt + 
coUWhhh*UWhhh + 
coUWhtt*UWhtt + 
coUWhWW*UWhWW + 
coUWhZZ*UWhZZ + 
coUWZ00*UWZ00 + 
coUWZhZ*UWZhZ + 
coUWZtt*UWZtt + 
coUWZWW*UWZWW + 
coUZW00*UZW00 + 
coUZW0t*UZW0t + 
coUZWhW*UZWhW + 
coUZWWZ*UZWWZ;

  return result;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Performs all needed basis integral evaluations.                    */

int SMDR_MW_DoTSILW (float loopOrder)
{
  TSIL_DATA bar;
  int success = 1;

  /* Do all 1-loop functions: */
  Ah = TSIL_A (h, Q2);
  AW = TSIL_A (W, Q2);
  AZ = TSIL_A (Z, Q2);
  At = TSIL_A (T, Q2);
  Ab = TSIL_A (b, Q2);
  Atau = TSIL_A (tau, Q2);
  B0t = TSIL_B (0, T, W, Q2);
  Bbt = TSIL_B (b, T, W, Q2);
  B0tau = TSIL_B (0, tau, W, Q2);
  B00 = TSIL_B (0, 0, W, Q2);
  BWZ = TSIL_B (W, Z, W, Q2);
  BhW = TSIL_B (h, W, W, Q2);

  /* Then do 2-loop evals if needed... */

  /* These are needed for the QCD contribution: */
  if (loopOrder > 1.1) {
    success *= TSIL_Manalytic (0, 0, T, T, 0, W, &M00tt0);
    success *= TSIL_Manalytic (0, 0, 0, 0, 0, W, &M00000W);
    success *= TSIL_Tanalytic (T, 0, 0, W, Q2, &Tt00);
  }

  /* These are needed for nonQCD part: */
  if (loopOrder > 1.6) {

    B0h = TSIL_B (0, h, W, Q2);
    B0Z = TSIL_B (0, Z, W, Q2);
    Bht = TSIL_B (h, T, W, Q2);
    BtZ = TSIL_B (T, Z, W, Q2);

    I00W = TSIL_I2 (0, 0, W, Q2);
    I00Z = TSIL_I2 (0, 0, Z, Q2);
    I00h = TSIL_I2 (0, 0, h, Q2);
    I00t = TSIL_I2 (0, 0, T, Q2);
    I0hW = TSIL_I2 (0, W, h, Q2);
    I0hZ = TSIL_I2 (0, Z, h, Q2);
    I0WZ = TSIL_I2 (0, W, Z, Q2);
    Ihtt = TSIL_I2 (T, T, h, Q2);
    IhWW = TSIL_I2 (W, W, h, Q2);
    IhZZ = TSIL_I2 (Z, Z, h, Q2);
    IttZ = TSIL_I2 (T, T, Z, Q2);
    I0tW = TSIL_I2 (0, T, W, Q2);
    IWWZ = TSIL_I2 (W, W, Z, Q2);
    Ihhh = TSIL_I2 (h, h, h, Q2);

    success *= TSIL_Manalytic (0, 0, 0, 0, Z, W, &M0000Z);
    success *= TSIL_Manalytic (0, 0, 0, W, 0, W, &M000W0);
    success *= TSIL_Manalytic (0, 0, T, W, 0, W, &M00tW0);
    success *= TSIL_Manalytic (0, W, W, 0, W, W, &M0WW0W);
    success *= TSIL_Manalytic (0, W, 0, Z, 0, W, &M0W0Z0);

    success *= TSIL_Uanalytic (h, W, 0, 0, W, Q2, &UhW00);
    success *= TSIL_Uanalytic (Z, W, 0, 0, W, Q2, &UZW00);

    success *= TSIL_Tanalytic (h, 0, 0, W, Q2, &Th00);
    success *= TSIL_Tanalytic (h, 0, T, W, Q2, &Th0t);
    success *= TSIL_Tanalytic (h, 0, W, W, Q2, &Th0W);
    success *= TSIL_Tanalytic (T, 0, h, W, Q2, &Tt0h);
    success *= TSIL_Tanalytic (T, 0, Z, W, Q2, &Tt0Z);
    success *= TSIL_Tanalytic (W, 0, 0, W, Q2, &TW00);
    success *= TSIL_Tanalytic (Z, 0, 0, W, Q2, &TZ00);
    success *= TSIL_Tanalytic (Z, 0, T, W, Q2, &TZ0t);
    success *= TSIL_Tanalytic (Z, 0, W, W, Q2, &TZ0W);

    success *= TSIL_Tbaranalytic (0, h, W, W, Q2, &Tbar0hW);
    success *= TSIL_Tbaranalytic (0, W, Z, W, Q2, &Tbar0WZ);

    if (1 != success) SMDR_Error ("SMDR_MW_DoTSILW", 
      "Analytic TSIL computation failed. This can't happen! Please contact the authors with details.", 1);

    TSIL_SetParameters (&bar, W, Z, Z, W, W, Q2);
    TSIL_Evaluate (&bar, W);
    MWZZWW = TSIL_GetFunction (&bar, "M");
    UZWWZ  = TSIL_GetFunction (&bar, "Uzxyv");
    UWZWW  = TSIL_GetFunction (&bar, "Uxzuv");
    SWZZ = TSIL_GetFunction (&bar, "Svyz");

    TSIL_SetParameters (&bar, W, W, Z, Z, h, Q2);
    TSIL_Evaluate (&bar, W);
    MWWZZh = TSIL_GetFunction (&bar, "M");
    UZWhW  = TSIL_GetFunction (&bar, "Uzxyv");
    UWZhZ  = TSIL_GetFunction (&bar, "Uxzuv");
    ShWZ = TSIL_GetFunction (&bar, "Svyz");
    TZhW   = TSIL_GetFunction (&bar, "Tuxv");
    ThWZ   = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, h, Z, W, W, Z, Q2);
    TSIL_Evaluate (&bar, W);
    MhZWWZ = TSIL_GetFunction (&bar, "M");
    UWhZZ  = TSIL_GetFunction (&bar, "Uzxyv");
    UhWWZ  = TSIL_GetFunction (&bar, "Uxzuv");

    TSIL_SetParameters (&bar, h, W, W, Z, W, Q2);
    TSIL_Evaluate (&bar, W);
    MhWWZW = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, h, W, W, h, W, Q2);
    TSIL_Evaluate (&bar, W);
    MhWWhW = TSIL_GetFunction (&bar, "M");
    UWhWW  = TSIL_GetFunction (&bar, "Uzxyv");
    UhWhW  = TSIL_GetFunction (&bar, "Uxzuv");
    ShhW = TSIL_GetFunction (&bar, "Suxv");

    TSIL_SetParameters (&bar, h, h, W, W, h, Q2);
    TSIL_Evaluate (&bar, W);
    MhhWWh = TSIL_GetFunction (&bar, "M");
    UWhhh  = TSIL_GetFunction (&bar, "Uzxyv");

    TSIL_SetParameters (&bar, 0, W, W, Z, W, Q2);
    TSIL_Evaluate (&bar, W);
    M0WWZW = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, 0, W, W, h, W, Q2);
    TSIL_Evaluate (&bar, W);
    M0WWhW = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, 0, Z, T, W, 0, Q2);
    TSIL_Evaluate (&bar, W);
    M0ZtW0 = TSIL_GetFunction (&bar, "M");
    UWZ00  = TSIL_GetFunction (&bar, "Uuyxv");
    UZW0t  = TSIL_GetFunction (&bar, "Uyuzv");
    U0t0W  = TSIL_GetFunction (&bar, "Uxzuv");

    TSIL_SetParameters (&bar, 0, W, T, Z, T, Q2);
    TSIL_Evaluate (&bar, W);
    M0WtZt = TSIL_GetFunction (&bar, "M");
    UWZtt  = TSIL_GetFunction (&bar, "Uyuzv");
    U0ttZ  = TSIL_GetFunction (&bar, "Uxzuv");

    TSIL_SetParameters (&bar, 0, W, T, h, T, Q2);
    TSIL_Evaluate (&bar, W);
    M0Wtht = TSIL_GetFunction (&bar, "M");
    UWhtt  = TSIL_GetFunction (&bar, "Uyuzv");
    U0tht  = TSIL_GetFunction (&bar, "Uxzuv");
    UhW0t  = TSIL_GetFunction (&bar, "Uuyxv");
    SttW = TSIL_GetFunction (&bar, "Svyz");

    TSIL_SetParameters (&bar, 0, T, W, 0, T, Q2);
    TSIL_Evaluate (&bar, W);
    M0tW0t = TSIL_GetFunction (&bar, "M");
    UW0tt  = TSIL_GetFunction (&bar, "Uzxyv");

    TSIL_SetParameters (&bar, 0, 0, T, T, Z, Q2);
    TSIL_Evaluate (&bar, W);
    M00ttZ = TSIL_GetFunction (&bar, "M");
  }

  return 0;
}
