/* Z boson complex pole mass calculation from 1505.04833 */

#include "smdr_internal.h"

#define Nc 3
#define nQ 3
#define nu 3
#define nd 3
#define nL 3
#define ne 3

/* Local functions; not needed or useful elsewhere: */
SMDR_COMPLEX SMDR_MZ_Z1loop (void);
SMDR_COMPLEX SMDR_MZ_Z1loopALT (void);
SMDR_COMPLEX SMDR_MZ_Z2loopQCD (void);
SMDR_COMPLEX SMDR_MZ_Z2loopnonQCD (void);
int SMDR_MZ_DoTSILZ (float);

/* #define global variables used in this file, for convenience with safety. */
#define AW SMDR_MZ_AW
#define AZ SMDR_MZ_AZ
#define Ah SMDR_MZ_Ah
#define At SMDR_MZ_At
#define Ab SMDR_MZ_Ab
#define Atau SMDR_MZ_Atau
#define B00 SMDR_MZ_B00
#define BhZ SMDR_MZ_BhZ
#define Btt SMDR_MZ_Btt
#define BWW SMDR_MZ_BWW
#define Bbb SMDR_MZ_Bbb
#define Btautau SMDR_MZ_Btautau
#define I00h SMDR_MZ_I00h
#define I00t SMDR_MZ_I00t
#define I00W SMDR_MZ_I00W
#define I00Z SMDR_MZ_I00Z
#define I0hW SMDR_MZ_I0hW
#define I0hZ SMDR_MZ_I0hZ
#define I0tW SMDR_MZ_I0tW
#define I0WZ SMDR_MZ_I0WZ
#define Ihhh SMDR_MZ_Ihhh
#define Ihtt SMDR_MZ_Ihtt
#define IhWW SMDR_MZ_IhWW
#define IhZZ SMDR_MZ_IhZZ
#define IttZ SMDR_MZ_IttZ
#define IWWZ SMDR_MZ_IWWZ
#define M00000Z SMDR_MZ_M00000Z
#define M0000W SMDR_MZ_M0000W
#define M0000Z SMDR_MZ_M0000Z
#define M0t0tW SMDR_MZ_M0t0tW
#define M0W0W0 SMDR_MZ_M0W0W0
#define M0W0Wt SMDR_MZ_M0W0Wt
#define MhtZtt SMDR_MZ_MhtZtt
#define MhWZWW SMDR_MZ_MhWZWW
#define MhhZZh SMDR_MZ_MhhZZh
#define MhZZhZ SMDR_MZ_MhZZhZ
#define Mtttt0 SMDR_MZ_Mtttt0
#define Mtttth SMDR_MZ_Mtttth
#define MttttZ SMDR_MZ_MttttZ
#define MtWtW0 SMDR_MZ_MtWtW0
#define MWWWW0 SMDR_MZ_MWWWW0
#define MWWWWh SMDR_MZ_MWWWWh
#define MWWWWZ SMDR_MZ_MWWWWZ
#define S00h SMDR_MZ_S00h
#define S00W SMDR_MZ_S00W
#define S0tW SMDR_MZ_S0tW
#define ShhZ SMDR_MZ_ShhZ
#define Shtt SMDR_MZ_Shtt
#define ShWW SMDR_MZ_ShWW
#define SttZ SMDR_MZ_SttZ
#define SWWZ SMDR_MZ_SWWZ
#define SZZZ SMDR_MZ_SZZZ
#define Th00 SMDR_MZ_Th00
#define ThhZ SMDR_MZ_ThhZ
#define Thtt SMDR_MZ_Thtt
#define ThWW SMDR_MZ_ThWW
#define Tt0W SMDR_MZ_Tt0W
#define Ttht SMDR_MZ_Ttht
#define TttZ SMDR_MZ_TttZ
#define TW00 SMDR_MZ_TW00
#define TW0t SMDR_MZ_TW0t
#define TWhW SMDR_MZ_TWhW
#define TWWZ SMDR_MZ_TWWZ
#define TZ00 SMDR_MZ_TZ00
#define Tbar0tt SMDR_MZ_Tbar0tt
#define Tbar0WW SMDR_MZ_Tbar0WW
#define UhZ00 SMDR_MZ_UhZ00
#define UhZhZ SMDR_MZ_UhZhZ
#define UhZtt SMDR_MZ_UhZtt
#define UhZWW SMDR_MZ_UhZWW
#define Utt0W SMDR_MZ_Utt0W
#define Uttht SMDR_MZ_Uttht
#define UtttZ SMDR_MZ_UtttZ
#define UWW00 SMDR_MZ_UWW00
#define UWW0t SMDR_MZ_UWW0t
#define UWWhW SMDR_MZ_UWWhW
#define UWWWZ SMDR_MZ_UWWWZ
#define UZhhh SMDR_MZ_UZhhh
#define UZhtt SMDR_MZ_UZhtt
#define UZhWW SMDR_MZ_UZhWW
#define UZhZZ SMDR_MZ_UZhZZ

SMDR_COMPLEX AW, AZ, Ah, At, Ab, Atau,
  B00, BhZ, Btt, BWW, Bbb, Btautau,
  I00h, I00t, I00W, I00Z, I0hW, I0hZ, I0tW, I0WZ,
  Ihhh, Ihtt, IhWW, IhZZ, IttZ, IWWZ,
  M00000Z, M0000W, M0000Z, M0t0tW, M0W0W0, M0W0Wt,
  MhtZtt, MhWZWW, MhhZZh, MhZZhZ, Mtttt0, Mtttth, 
  MttttZ, MtWtW0, MWWWW0, MWWWWh, MWWWWZ,
  S00h, S00W, S0tW, ShhZ, Shtt, ShWW, SttZ, SWWZ, SZZZ,
  Th00, ThhZ, Thtt, ThWW, Tt0W, Ttht, TttZ, TW00, TW0t, 
  TWhW, TWWZ, TZ00, Tbar0tt, Tbar0WW, 
  UhZ00, UhZhZ, UhZtt, UhZWW, Utt0W, 
  Uttht, UtttZ, UWW00, UWW0t, UWWhW, 
  UWWWZ, UZhhh, UZhtt, UZhWW, UZhZZ;

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/*  Computes Z mass at up to two loops, from 1505.04833. 
    loopOrder may take the following values:

    0    tree level
    1    1-loop
    1.5  1-loop plus 2-loop QCD corrections
    2    full 2-loop
*/

void SMDR_Eval_MZ_pole (SMDR_REAL Q_eval,
                        float loopOrder,
                        SMDR_REAL *MZpoleresult,
                        SMDR_REAL *GammaZpoleresult,
                        SMDR_REAL *MZBreitWignerresult,
                        SMDR_REAL *GammaZBreitWignerresult)
{
  SMDR_COMPLEX CM2Z;
  char funcname[] = "SMDR_Eval_MZ_pole";

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
  CM2Z = Z;

  /* Perform needed TSIL evaluations: */
  SMDR_MZ_DoTSILZ (loopOrder);

  if (loopOrder > 0.999) {
    CM2Z += ONELOOPFACTOR * SMDR_MZ_Z1loop ();
  }

  if (loopOrder > 1.499) {
    CM2Z += TWOLOOPFACTOR * SMDR_MZ_Z2loopQCD ();
  }

  if (loopOrder > 1.999) {
    CM2Z += TWOLOOPFACTOR * SMDR_MZ_Z2loopnonQCD ();
  }

  *MZpoleresult = SMDR_SQRT(SMDR_CREAL(CM2Z));

  *GammaZpoleresult = -SMDR_CIMAG(CM2Z)/(*MZpoleresult);

  SMDR_MZ_BreitWigner = SMDR_SQRT(
    (*MZpoleresult) * (*MZpoleresult) +
    (*GammaZpoleresult) * (*GammaZpoleresult));

  *GammaZBreitWignerresult = (*GammaZpoleresult) * (*MZpoleresult)/
    (*MZBreitWignerresult);

   return;
}

/* ------------------------------------------------------------ */
/* From 1505.04833, eqs. (2.22)-(2.25)                          */

SMDR_COMPLEX SMDR_MZ_Z1loop ()
{
  SMDR_COMPLEX result;

/* Version with yb = ytau = 0:
  result = g2 * (
      (4*W-Z)*(W/Z + 5./3. + Z/(12*W))*BWW +
      (4*h*Z - 12*Z2 - h2)/(12*W) * BhZ +
      ((h-2*Z)/(12*W)) * AZ +
      (4*W/Z - 4./3. - Z/(6*W)) * AW + 
      (3*Z - h)/(12*W) * Ah +
      4*W2/Z - 4*W/3 + 5*Z/9 + h*Z/(6*W) + Z2/(18*W)) +
    Nc*(auL2 + auR2)*((2./3.)*(T-Z)*Btt - (4./3.)*At + (2./9.)*(Z-6*T)) +
    -Nc*4*auL*auR*T*Btt +
    (Nc*(nQ-1)*auL2 + Nc*nQ*adL2 + Nc*(nu-1)*auR2 + Nc*nd*adR2 +
       nL*(aeL2 + anL2) + ne*aeR2)*Z*(2./9. - (2./3.)*B00);
*/

  result = g2 * (
      (4*W-Z)*(W/Z + 5./3. + Z/(12*W))*BWW +
      (4*h*Z - 12*Z2 - h2)/(12*W) * BhZ +
      ((h-2*Z)/(12*W)) * AZ +
      (4*W/Z - 4./3. - Z/(6*W)) * AW + 
      (3*Z - h)/(12*W) * Ah +
      4*W2/Z - 4*W/3 + 5*Z/9 + h*Z/(6*W) + Z2/(18*W)) +
     Nc*(auL2 + auR2)*((2./3.)*(T-Z)*Btt - (4./3.)*At + (2./9.)*(Z-6*T)) +
    -Nc*4*auL*auR*T*Btt +
     Nc*(adL2 + adR2)*((2./3.)*(b-Z)*Bbb - (4./3.)*Ab + (2./9.)*(Z-6*b)) +
    -Nc*4*adL*adR*b*Bbb +
     (aeL2 + aeR2)*((2./3.)*(tau-Z)*Btautau -(4./3.)*Atau + (2./9.)*(Z-6*tau)) +
    -4*aeL*aeR*tau*Btautau +
    (Nc*(nQ-1)*auL2 + Nc*(nQ-1)*adL2 + Nc*(nu-1)*auR2 + Nc*(nd-1)*adR2 +
       (nL-1)*aeL2 + nL*anL2 + (ne-1)*aeR2)*Z*(2./9. - (2./3.)*B00);

  return result;
}

/* ------------------------------------------------------------ */

SMDR_COMPLEX SMDR_MZ_Z1loopALT ()
{
  SMDR_COMPLEX result, TMSbar;
  SMDR_COMPLEX deltaT, deltaT2, deltaT3;
  /* SMDR_COMPLEX deltaT4, deltaT5; */
  SMDR_COMPLEX f1T, f2T;

  TMSbar = T;
  T = 173.34 * 173.34;
  T2 = T*T;
  T3 = T*T2;
  T4 = T*T3;
  SMDR_MZ_DoTSILZ(1);

  deltaT = (TMSbar - T);
  deltaT2 = deltaT * deltaT;
  deltaT3 = deltaT * deltaT2;
  /* deltaT4 = deltaT * deltaT3; */
  /* deltaT5 = deltaT * deltaT4; */
  
  f1T = (-4*At)/3 - (4*T)/3 + (2*Btt*(T - Z))/3 + (2*Z)/9;
  f2T = -2*Btt*T;

  f1T += deltaT * (-4*At + Btt*(4*T - 2*Z) + 4*(-3*T + Z))/(4*T - Z);
  f2T += deltaT * (-4*At + 4*T + 2*Btt*(-6*T + Z))/(4*T - Z);

  f1T += (1/2.) * deltaT2 *
         (4*(2*At - 6*T + 2*Btt*T + Z))/((-4*T + Z)*(-4*T + Z));
  f2T += (1/2.) * deltaT2 *
         (-4*(At*(6*T - 2*Z) + T*(-2*T + Btt*(6*T - 2*Z) + Z)))/
         (T*(-4*T + Z)*(-4*T + Z));

  f1T += (1/6.) * deltaT3 *
         (-8*(At*(2*T + Z) + T*(2*(-7*T + Z) + Btt*(2*T + Z))))/
         (T*(4*T - Z)*(4*T - Z)*(4*T - Z));
  f2T += (1/6.) * deltaT3 *
         (8*(-10*T2 + At*(6*T - 3*Z) + 3*Btt*T*(2*T - Z) + 8*T*Z - Z2))/
         (T*(4*T - Z)*(4*T - Z)*(4*T - Z));

/*
  f1T += (1/24.) * deltaT4 *
         (8*(-116*T2 + 10*T*Z + Z2 + 12*At*(T + Z) + 12*Btt*T*(T + Z)))/
         (T*(-4*T + Z)*(-4*T + Z)*(-4*T + Z)*(-4*T + Z));
  f2T += (1/24.) * deltaT4 *
         (-8*(-92*T3 + 12*At*T*(3*T - 2*Z) + 12*Btt*T2*(3*T - 2*Z) + 98*T2*Z - 
         19*T*Z2 + Z3))/(T2*(-4*T + Z)*(-4*T + Z)*(-4*T + Z)*(-4*T + Z));

  f1T += (1/120.) * deltaT5 *
         (-8*(-1416*T3 + 32*T2*Z + 32*T*Z2 - Z3 + 60*At*T*(2*T + 3*Z) + 
         60*Btt*T2*(2*T + 3*Z)))/
         (T2*(4*T - Z)*(4*T - Z)*(4*T - Z)*(4*T - Z)*(4*T - Z));
  f2T += (1/120.) * deltaT5 *
         (8*(-1176*T4 + 60*At*T2*(6*T - 5*Z) + 60*Btt*T3*(6*T - 5*Z) + 
         1560*T3*Z - 404*T2*Z2 + 43*T*Z3 - 2*Z4))/
         (T3*(4*T - Z)*(4*T - Z)*(4*T - Z)*(4*T - Z)*(4*T - Z));
*/
           
  result = g2 * (
      (4*W-Z)*(W/Z + 5./3. + Z/(12*W))*BWW +
      (4*h*Z - 12*Z2 - h2)/(12*W) * BhZ +
      ((h-2*Z)/(12*W)) * AZ +
      (4*W/Z - 4./3. - Z/(6*W)) * AW +
      (3*Z - h)/(12*W) * Ah +
      4*W2/Z - 4*W/3 + 5*Z/9 + h*Z/(6*W) + Z2/(18*W)) +
    Nc*(auL2 + auR2)*f1T +Nc*2*auL*auR*f2T +
    (Nc*(nQ-1)*auL2 + Nc*nQ*adL2 + Nc*(nu-1)*auR2 + Nc*nd*adR2 +
       nL*(aeL2 + anL2) + ne*aeR2)*Z*(2./9. - (2./3.)*B00);
  
  T = TMSbar;
  T2 = T*T;
  T3 = T*T2;
  T4 = T*T3;
  SMDR_MZ_DoTSILZ(1);

  return result;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Two-loop Z boson pole mass QCD contribution.                        */
/* From 1505.04833 eqs. (2.26)-(2.29)                                  */

SMDR_COMPLEX SMDR_MZ_Z2loopQCD (void)
{
  SMDR_COMPLEX result, F1ttZ, F2ttZ, F100Z;

  F1ttZ = -(8*(T-Z)*(2*T-Z)/3)*Mtttt0 + (16*(Z-T)/3)*Tbar0tt +
          1/(3*Z*(4*T-Z)) * (
            4*(-6*T3 + 6*T2*Z + 5*T*Z2 - 2*Z3)*Btt*Btt +
            -16*(3*T2 + 2*T*Z - 2*Z2)*Btt*At +
            -4*(T-Z)*(12*T2 - 30*T*Z + 7*Z2)*Btt +
            (-24*T + 56*Z + 16*Z2/T)*At*At +
            (-48*T2 + 296*T*Z - 104*Z2)*At +        
            -24*T3 + 220*T2*Z - 141*T*Z2 + 23*Z3);

  F2ttZ = 8*T*(2*T-Z)*Mtttt0 + 16*T*Tbar0tt +
          1/(Z*(4*T-Z)) * (
            4*T*(2*T2 + Z2)*Btt*Btt + 16*T*(T+5*Z)*Btt*At +
            4*T*(4*T2 - 50*T*Z + 9*Z2)*Btt + (8*T+64*Z)*At*At +
            8*T*(2*T-13*Z)*At + 8*T3 - 140*T2*Z + 43*T*Z2);

  F100Z = -(8*Z2/3)*M00000Z - 4*Z*B00 - 31*Z/3;

  result = g32 * ((Nc*Nc-1)/4.) * ((auL2 + auR2)*F1ttZ + 2*auL*auR*F2ttZ +
           ((nQ-1)*auL2 + (nu-1)*auR2 + nQ*adL2 + nd*adR2)*F100Z);

  return (result);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/*  Two-loop non-QCD contribution to the Z boson pole mass.           */ 
/* From 1505.04833 ancillary file "coefficients.txt" and eq. (2.30)   */

SMDR_COMPLEX SMDR_MZ_Z2loopnonQCD (void) 
{
SMDR_COMPLEX result;

SMDR_REAL co1, coAh, coAhAh, coAhAt, coAhAW, coAhAZ, 
 coAhB00, coAhBhZ, coAhBtt, coAhBWW, coAt, coAtAt, coAtAW, coAtAZ, coAtB00, 
 coAtBhZ, coAtBtt, coAtBWW, coAW, coAWAW, coAWAZ, coAWB00, coAWBhZ, coAWBtt, 
 coAWBWW, coAZ, coAZAZ, coAZB00, coAZBhZ, coAZBtt, coAZBWW, coB00, coB00B00, 
 coB00BhZ, coB00Btt, coB00BWW, coBhZ, coBhZBhZ, coBhZBtt, coBhZBWW, coBtt, 
 coBttBtt, coBttBWW, coBWW, coBWWBWW, coI00h, coI00t, coI00W, coI00Z, coI0hW, 
 coI0hZ, coI0tW, coI0WZ, coIhhh, coIhtt, coIhWW, coIhZZ, coIttZ, coIWWZ, 
 coM00000, coM0000W, coM0000Z, coM0t0tW, coM0W0W0, coM0W0Wt, coMhhZZh, 
 coMhtZtt, coMhWZWW, coMhZZhZ, coMtttt0, coMtttth, coMttttZ, coMtWtW0, 
 coMWWWW0, coMWWWWh, coMWWWWZ, coS00h, coS00W, coS0tW, coShhZ, coShtt, coShWW, 
 coSttZ, coSWWZ, coSZZZ, coTbar0tt, coTbar0WW, coTh00, coThhZ, coThtt, coThWW, 
 coTt0W, coTtht, coTttZ, coTW00, coTW0t, coTWhW, coTWWZ, coTZ00, coUhZ00, 
 coUhZhZ, coUhZtt, coUhZWW, coUtt0W, coUttht, coUtttZ, coUWW00, coUWW0t, 
 coUWWhW, coUWWWZ, coUZhhh, coUZhtt, coUZhWW, coUZhZZ;

co1 = ((57595*h3)/4608 - (1687*h2*T)/144 - (929*h*T2)/72 - 14*T3 + (515*h2*W)/16 - 
  (1810*h*T*W)/27 + (1414*T2*W)/81 + (13703*h*W2)/432 - (7486*T*W2)/9 - 
  (32738*W3)/27 + (160*W4)/(9*T) - (22747*h4 + 800*h3*T + 12224*h2*T2 - 
    55296*h*T3 - 59168*h3*W - 56320*h2*T*W - 163840*h*T2*W + 107936*h2*W2 + 
    180224*h*T*W2 + 524288*T2*W2 - 310272*h*W3 - 589824*W4)/(4608*(h - 4*Z)) - 
  (11*h4)/(12*Z) + (64*T4)/(3*Z) - (43*h3*W)/(3*Z) + (400*h2*T*W)/(9*Z) + 
  (160*h*T2*W)/(9*Z) - (160*T3*W)/(3*Z) - (176*h2*W2)/(9*Z) + 
  (1184*h*T*W2)/(27*Z) - (21292*T2*W2)/(27*Z) - (92*h*W3)/(9*Z) + 
  (57244*T*W3)/(81*Z) + (124696*W4)/(81*Z) - (128*W5)/(3*T*Z) - 
  (15773*h2*Z)/1152 + (1997*h*T*Z)/108 + (7045*T2*Z)/54 - (24*T3*Z)/h + 
  (2*h3*Z)/(3*W) - (4*T3*Z)/W - (3026*h*W*Z)/27 + (78064*T*W*Z)/81 + 
  (97963*W2*Z)/108 + (12*W3*Z)/h + (872*W3*Z)/(9*T) - 
  (T2*(16*T - 5*Z)*(2*T + Z)*(2*T + Z))/(3*(T - W)*Z) - (16*h3*W2)/(3*Z2) - 
  (320*h2*T*W2)/(9*Z2) - (128*h*T2*W2)/(9*Z2) + (496*T3*W2)/(3*Z2) + 
  (104*h2*W3)/(3*Z2) + (35552*T2*W3)/(27*Z2) + (80*h*W4)/(3*Z2) - 
  (103504*T*W4)/(81*Z2) + (11488*W5)/(9*Z2) + (64961*h*Z2)/864 - 
  (24496*T*Z2)/81 + (30*T2*Z2)/h - (h2*Z2)/(6*W) + (T2*Z2)/W - 
  (161617*W*Z2)/324 + (W2*Z2)/h - (394*W2*Z2)/(3*T) + 
  (h2*(12*h3 - 55*h2*Z + 17*h*Z2 - 64*Z3))/(48*(h - 4*W)*Z) - 
  (40960*T2*W4)/(81*Z3) + (2048*T*W5)/(3*Z3) + (128*W6)/Z3 + (106841*Z3)/2592 + 
  (560*W*Z3)/(9*T) + (13*Z4)/(2*h) - (107*Z4)/(9*T) + Z4/(2*W) + 
  (21*Z4)/(32*(4*W - Z)) + (2688*h2*W2 + 43776*W4 - 3360*h2*W*Z - 3456*h*W2*Z - 
    106752*W3*Z + 1050*h2*Z2 + 4320*h*W*Z2 + 97876*W2*Z2 - 1350*h*Z3 - 
    39869*W*Z3 + 5950*Z4)/(216*(4*T - Z)))/v4;

coAh = ((-6959*h2)/576 + (179*h*T)/54 - (38*T2)/9 + (1949*h*W)/54 + 
      (1840*T*W)/27 - (440*W2)/9 + (h*(2935*h2 - 864*h*T - 9152*h*W + 
         15488*W2))/(576*(h - 4*Z)) + (11*h3)/(12*Z) + (43*h2*W)/(3*Z) - 
      (560*h*T*W)/(27*Z) + (640*T2*W)/(9*Z) - (964*h*W2)/(27*Z) - 
      (1472*T*W2)/(27*Z) + (880*W3)/(9*Z) - (2027*h*Z)/432 - (782*T*Z)/27 - 
      (12*T2*Z)/h - (2*h2*Z)/(3*W) - (343*W*Z)/18 + (22*W2*Z)/h + 
      ((-8*W + 5*Z)*(-8*W + 5*Z)*(-5*h + 6*Z))/(36*(4*T - Z)) - 
      (h2*(6*h2 - 5*h*Z - 13*Z2))/(24*(h - 4*W)*Z) + (16*h2*W2)/(3*Z2) + 
      (448*h*T*W2)/(27*Z2) - (512*T2*W2)/(9*Z2) - (24*h*W3)/Z2 + (64*W4)/Z2 - 
      (257*Z2)/36 + (11*Z3)/h)/v4;
 
coAt = ((-17*h2)/12 + 20*h*T + 18*T2 - (880*h*W)/27 + (20096*T*W)/81 - 
      368*W2 - (11968*W3)/(27*T) + (51*h3 + 136*h2*T - 432*h*T2 - 480*h2*W - 
        1280*h*T*W + 1536*h*W2 + 4096*T*W2)/(36*(h - 4*Z)) - (32*T3)/(3*Z) - 
      (800*h*T*W)/(9*Z) + (64*T2*W)/Z + (992*h*W2)/(27*Z) - 
      (19328*T*W2)/(27*Z) + (26992*W3)/(81*Z) + (5248*W4)/(27*T*Z) + 
      (829*h*Z)/54 - (3731*T*Z)/54 + (24*T2*Z)/h + (4*T2*Z)/W + 
      (89635*W*Z)/162 + (3592*W2*Z)/(9*T) + (8*T2*(2*T + Z)*(2*T + Z))/
       (3*(T - W)*Z) + (640*h*T*W2)/(9*Z2) - (496*T2*W2)/(3*Z2) + 
      (2576*T*W3)/(3*Z2) - (72032*W4)/(81*Z2) - (137441*Z2)/648 - 
      (4480*W*Z2)/(27*T) - (32768*T*W4)/(81*Z3) + (2048*W5)/(3*Z3) + 
      (856*Z3)/(27*T) + (73728*W5 - 61440*W4*Z - 2304*h*W2*Z2 - 80000*W3*Z2 + 
        2880*h*W*Z3 + 112160*W2*Z3 - 900*h*Z4 - 43708*W*Z4 + 6325*Z5)/
       (216*(4*T - Z)*Z2))/v4;
 
coAW = ((-113*h2)/8 + (310*T2)/9 + 27*h*W - (1712*T*W)/27 + (65606*W2)/27 - 
      (64*W3)/T + (9*h3 + 80*h2*W - 1712*h*W2 - 3072*W3)/(24*(h - 4*Z)) + 
      h3/(2*Z) - (32*T3)/(3*Z) + (86*h2*W)/(3*Z) - (704*T2*W)/(9*Z) - 
      (704*h*W2)/(9*Z) + (7552*T*W2)/(27*Z) - (16928*W3)/(27*Z) + 
      (128*W4)/(3*T*Z) - (37*h*Z)/36 - (581*T*Z)/54 - (4*h2*Z)/(3*W) + 
      (4*T2*Z)/W - (24295*W*Z)/54 - (12*W2*Z)/h + (88*W2*Z)/(3*T) + 
      (8*T2*(2*T + Z)*(2*T + Z))/(3*(T - W)*Z) - (h2*(6*h2 - 5*h*Z - 13*Z2))/
       (12*(h - 4*W)*Z) + (32*h2*W2)/(3*Z2) - (464*T2*W2)/(9*Z2) - 
      (272*h*W3)/(3*Z2) - (26992*T*W3)/(27*Z2) - (4864*W4)/(9*Z2) - 
      (13955*Z2)/216 + ((-8*W + 5*Z)*(-8*W + 5*Z)*(-192*W2 + 92*W*Z + Z2))/
       (72*(4*T - Z)*Z) + (2048*T*W4)/(3*Z3) - (1280*W5)/Z3 - (4*Z3)/(3*W) - 
      (21*Z3)/(4*(4*W - Z)))/v4;
 
coAZ = (-h2/8 + (259*h*T)/54 - 4*T2 - (799*h*W)/27 + (3040*T*W)/81 - 
      (35561*W2)/27 - (5120*W3)/(27*T) + (17*h3 + 380*h2*T - 432*h*T2 + 
        488*h2*W - 2560*h*T*W - 776*h*W2 + 8192*T*W2 - 9216*W3)/
       (72*(h - 4*Z)) + h3/(3*Z) - (640*h*T*W)/(27*Z) + (520*h*W2)/(27*Z) + 
      (1408*T*W2)/(9*Z) + (25024*W3)/(27*Z) + (2048*W4)/(27*T*Z) + 
      (317*h*Z)/27 - (2146*T*Z)/81 + (24*T2*Z)/h + (11236*W*Z)/9 - 
      (28*W2*Z)/h + (1664*W2*Z)/(9*T) + (512*h*T*W2)/(27*Z2) - 
      (64*h*W3)/(3*Z2) - (20480*T*W3)/(81*Z2) - (32320*W4)/(27*Z2) - 
      (20101*Z2)/54 - (2240*W*Z2)/(27*T) - 
      ((-8*W + 5*Z)*(-8*W + 5*Z)*(64*W2 + 3*h*Z - 80*W*Z + Z2))/(54*(4*T - Z)*Z) + 
      (8192*T*W4)/(81*Z3) - (256*W5)/Z3 - (20*Z3)/h + (428*Z3)/(27*T) - 
      (2*Z3)/(3*W))/v4;
 
coB00 = ((-103*h3)/288 - (16*T3)/3 + (25*h2*W)/9 - (16*T2*W)/3 - 
      (560*h*W2)/27 - (7864*T*W2)/9 + (191696*W3)/81 + 
      (h2*(103*h2 - 800*h*W + 2560*W2))/(288*(h - 4*Z)) + 
      (133120*T*W3)/(81*Z) - (171904*W4)/(81*Z) - (103*h2*Z)/72 - 
      (16*T2*Z)/3 + (700*h*W*Z)/27 - (3032*T*W*Z)/27 - (97832*W2*Z)/27 + 
      (4*T2*(2*T + Z)*(2*T + Z))/(3*(T - W)) - (20480*T*W4)/(27*Z2) + 
      (2560*W5)/(3*Z2) - (721*h*Z2)/54 + (12728*T*Z2)/81 + (232972*W*Z2)/81 + 
      ((-8*W + 5*Z)*(-8*W + 5*Z)*(160*W2 - 200*W*Z + 103*Z2))/(54*(4*T - Z)) - 
      (153103*Z3)/162)/v4;
 
coBhZ = ((-77*h3)/192 - (493*h2*T)/54 + 4*h*T2 + (700*h2*W)/27 - 
      (1840*h*T*W)/27 + (1621*h*W2)/27 + (256*T*W2)/9 - 32*W3 + 
      (h*(239*h3 - 544*h2*T - 1120*h2*W + 5120*h*T*W + 1984*h*W2 - 
         16384*T*W2 + 18432*W3))/(576*(h - 4*Z)) + (640*h2*T*W)/(27*Z) - 
      (520*h2*W2)/(27*Z) + (1280*h*T*W2)/(27*Z) - (160*h*W3)/(3*Z) - 
      (5693*h2*Z)/432 + (782*h*T*Z)/27 - 16*T2*Z - (2644*h*W*Z)/27 - 
      (320*T*W*Z)/9 - (748*W2*Z)/3 - (512*h2*T*W2)/(27*Z2) + 
      (64*h2*W3)/(3*Z2) + (4133*h*Z2)/108 + (136*T*Z2)/9 + (48*T2*Z2)/h + 
      (1298*W*Z2)/3 - (40*W2*Z2)/h - (h2*(h2 - 4*h*Z + 12*Z2))/
       (8*(h - 4*W)) + ((-8*W + 5*Z)*(-8*W + 5*Z)*(h2 - 4*h*Z + 12*Z2))/(18*(4*T - Z)) - 
      241*Z3 - (20*Z4)/h)/v4;
 
coBtt = ((-17*h3)/288 - (7*h2*T)/72 + (5*h2*W)/9 - (200*h*T*W)/27 + 
      (2752*T2*W)/9 - (16*h*W2)/27 - (3152*T*W2)/9 - (740*W3)/81 + 
      (h*(17*h3 + 28*h2*T - 160*h2*W - 1280*h*T*W + 512*h*W2 + 4096*T*W2))/
       (288*(h - 4*Z)) + (256*h*T*W2)/(27*Z) - (10016*T2*W2)/(9*Z) + 
      (20192*T*W3)/(81*Z) - (26944*W4)/(81*Z) - (17*h2*Z)/72 + 
      (35*h*T*Z)/54 + (1298*T2*Z)/81 + (20*h*W*Z)/27 + (43838*T*W*Z)/81 - 
      (5086*W2*Z)/27 + (114688*T2*W3)/(81*Z2) - (80384*T*W4)/(81*Z2) + 
      (256*W5)/Z2 - (22*h*Z2)/27 - (11977*T*Z2)/162 + (53255*W*Z2)/162 - 
      (16384*T2*W4)/(27*Z3) + (2048*T*W5)/(3*Z3) - (86069*Z3)/648 + 
      (18432*W5 + 30208*W4*Z + 768*h*W2*Z2 - 131808*W3*Z2 - 960*h*W*Z3 + 
        132944*W2*Z3 + 300*h*Z4 - 55532*W*Z4 + 8825*Z5)/(216*(4*T - Z)*Z))/v4;
 
coBWW = ((45*h3)/64 + (5*h2*W)/2 + (149*h*W2)/9 + (1468*T*W2)/3 + 
      (47608*W3)/27 + (h*(h - 16*W)*(h2 + 80*h*W + 192*W2))/(192*(h - 4*Z)) + 
      h4/(8*Z) + (h3*W)/(2*Z) + (2*h2*W2)/Z + (8*h*W3)/(3*Z) - 
      (22000*T*W3)/(27*Z) + (29152*W4)/(27*Z) - (11*h2*Z)/16 - (52*h*W*Z)/9 - 
      (1664*T*W*Z)/27 + (13778*W2*Z)/27 - (9472*T*W4)/(27*Z2) - (928*W5)/Z2 + 
      (11*h*Z2)/36 - (47*T*Z2)/27 - (21287*W*Z2)/108 + 
      ((8*W - 5*Z)*(8*W - 5*Z)*(4*W - Z)*(12*W2 + 20*W*Z + Z2))/(36*Z*(-4*T + Z)) - 
      (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(24*(h - 4*W)*Z) - (768*W6)/Z3 + 
      (2048*T*W5)/(3*Z3) - (6539*Z3)/432 - (21*Z4)/(16*(4*W - Z)))/v4;
 
coI00h = (-3*h*Z)/v4;
 
coI00t = (-24*T*Z)/v4;
 
coI00W = ((-932*W2)/9 - (758*W3)/(9*Z) + (622*W*Z)/9 + (1184*W4)/(3*Z2))/v4;
 
coI00Z = ((-103*h2)/72 + (100*h*W)/9 + (89446*W2)/27 + 
      (h*(103*h2 - 800*h*W + 2560*W2))/(72*(h - 4*Z)) - (208640*W3)/(81*Z) - 
      (103*h*Z)/18 - (162572*W*Z)/81 + (74240*W4)/(81*Z2) + (40460*Z2)/81)/v4;
 
coI0hW = (6*h*Z + 6*W*Z)/v4;
 
coI0hZ = (3*h*Z + 3*Z2)/v4;
 
coI0tW = ((388*T2)/9 - (124*T*W)/9 - 32*W2 - (64*W3)/T - 
      (608*T2*W)/(9*Z) - (256*T*W2)/(9*Z) - (32*W3)/Z + (128*W4)/(3*T*Z) + 
      (275*T*Z)/18 + (4*T2*Z)/W + (175*W*Z)/18 + (88*W2*Z)/(3*T) + 
      ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(24*(4*T - Z)*Z) - 
      (464*T2*W2)/(9*Z2) - (464*T*W3)/(9*Z2) + (928*W4)/(9*Z2) + (25*Z2)/24)/v4;
 
coI0WZ = (18*W2 + (6*W3)/Z + 18*W*Z + 6*Z2)/v4;
 
coIhhh = ((11*h2)/8 + (9*h3)/(8*(h - 4*Z)) - 6*h*Z)/v4;
 
coIhtt = ((145*h*T)/18 - (182*T2)/9 - (40*h*W)/3 + (32*W2)/3 + 
      (3*h*(h - 4*T)*T)/(2*(h - 4*Z)) - (160*h*T*W)/(9*Z) + 
      (640*T2*W)/(9*Z) + (32*h*W2)/(3*Z) + (43*h*Z)/6 - 24*T*Z + 
      (36*T2*Z)/h - (40*W*Z)/3 + ((-h + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(6*(4*T - Z)) + 
      (128*h*T*W2)/(9*Z2) - (512*T2*W2)/(9*Z2) + (25*Z2)/6)/v4;
 
coIhWW = ((-137*h2)/24 + (71*h*W)/6 - (107*W2)/3 + 
      (h*(h2 - 4*h*W + 12*W2))/(4*(h - 4*Z)) - h3/(4*Z) + (37*h2*W)/(3*Z) - 
      (92*h*W2)/(3*Z) + (176*W3)/(3*Z) - (67*h*Z)/24 - (2*h2*Z)/(3*W) + 
      (19*W*Z)/2 - (18*W2*Z)/h + (h2*(6*h2 - 5*h*Z - 13*Z2))/
       (24*(h - 4*W)*Z) + (16*h2*W2)/(3*Z2) - (64*h*W3)/(3*Z2) + (64*W4)/Z2)/v4;
 
coIhZZ = ((-15*h2)/64 + (77*h3)/(192*(h - 4*Z)) + h3/(3*Z) - 
      (133*h*Z)/48 + (71*Z2)/12 - (9*Z3)/h)/v4;
 
coIttZ = ((-17*h2)/72 + (7*h*T)/18 - (800*T*W)/81 - (5120*W3)/(27*T) + 
      (17*h3 + 28*h2*T - 160*h2*W - 1280*h*T*W + 512*h*W2 + 4096*T*W2)/
       (72*(h - 4*Z)) - (80*h*T*W)/(9*Z) + (16*h*W2)/(9*Z) + 
      (4864*T*W2)/(27*Z) + (2048*W4)/(27*T*Z) - (h*Z)/4 - (352*T*Z)/81 + 
      (1664*W2*Z)/(9*T) + (64*h*T*W2)/(9*Z2) - (20480*T*W3)/(81*Z2) - 4*Z2 - 
      (2240*W*Z2)/(27*T) - ((-8*W + 5*Z)*(-8*W + 5*Z)*(512*W2 - 3*h*Z - 640*W*Z + 
         224*Z2))/(108*(4*T - Z)*Z) + (8192*T*W4)/(81*Z3) + (428*Z3)/(27*T))/v4;
 
coIWWZ = (-h2/48 + (2*h*W)/3 + 536*W2 + 
      ((h - 16*W)*(h2 + 80*h*W + 192*W2))/(48*(h - 4*Z)) - (8*h*W2)/(3*Z) - 
      (608*W3)/(3*Z) - (17*h*Z)/12 - (196*W*Z)/3 - (8*h*W3)/Z2 - 
      (1408*W4)/(3*Z2) - 14*Z2 - (256*W5)/Z3 - (2*Z3)/(3*W))/v4;
 
coM00000 = ((14848*W4)/81 - (3968*W3*Z)/9 + (9920*W2*Z2)/27 - 
      (8896*W*Z3)/81)/v4;
 
coM0000W = (96*W4 + (256*W5)/(9*Z) + (784*W3*Z)/9 - (176*W*Z3)/9)/v4;
 
coM0000Z = ((-59392*W4)/81 + (166912*W3*Z)/81 - (71168*W2*Z2)/27 + 
      (129280*W*Z3)/81 - (32368*Z4)/81)/v4;
 
coM0t0tW = ((-32*T3*W)/9 + (160*T2*W2)/9 - (224*T*W3)/9 + 32*W4 + 
      (64*T3*W2)/(9*Z) - (64*T*W4)/(3*Z) + (128*W5)/(9*Z) - (32*T3*Z)/9 + 
      (64*T2*W*Z)/9 - (128*T*W2*Z)/9 + (176*W3*Z)/9 - (8*T2*Z2)/9 - 
      (32*T*W*Z2)/9 - (16*W*Z3)/9)/v4;
 
coM0W0W0 = (384*W4 + (192*W5)/Z)/v4;
 
coM0W0Wt = ((64*T2*W2)/3 - (128*T*W3)/3 + (160*W4)/3 + (32*T3*W2)/(3*Z) - 
      (32*T*W4)/Z + (64*W5)/(3*Z) - (8*T3*Z)/3 + (16*T2*W*Z)/3 + 8*T*W2*Z + 
      (64*W3*Z)/3 - (8*T2*Z2)/3 + (32*T*W*Z2)/3)/v4;
 
coMhhZZh = (-h4 + 12*h2*Z2 - 24*h*Z3)/v4;
 
coMhtZtt = (8*h2*T2 + (1024*h*T*W2)/9 - (2048*T2*W2)/9 - 32*h*T2*Z - 
      (1280*h*T*W*Z)/9 + (2560*T2*W*Z)/9 - (1024*T*W2*Z)/9 + (400*h*T*Z2)/9 - 
      (224*T2*Z2)/9 + (1280*T*W*Z2)/9 - (544*T*Z3)/9)/v4;
 
coMhWZWW = ((-8*h3*W)/3 + (16*h2*W2)/3 - (736*h*W3)/3 + 256*W4 + 
      (32*h2*W*Z)/3 + (128*h*W2*Z)/3 + (1088*W3*Z)/3 - (4*h2*Z2)/3 - 
      (8*h*W*Z2)/3 - (256*W2*Z2)/3 - (8*h*Z3)/3 - (16*W*Z3)/3)/v4;
 
coMhZZhZ = (-h4/3 - (8*h3*Z)/3 + (8*h2*Z2)/3 + (32*h*Z3)/3 - 16*Z4)/v4;
 
coMtttt0 = ((-3712*T2*W2)/27 + (2048*W4)/81 + (2048*T2*W3)/(9*Z) + 
      (896*T2*W*Z)/81 - (64*T*W2*Z)/3 - (512*W3*Z)/9 - (8192*T2*W4)/(81*Z2) + 
      (64*T*W*Z2)/3 + (1216*W2*Z2)/27 - (1088*W*Z3)/81)/v4;
 
coMtttth = ((80*h2*T*W)/9 - (160*h*T2*W)/3 + (640*T3*W)/9 - (256*T2*W2)/9 - 
      (64*h2*T*W2)/(9*Z) + (128*h*T2*W2)/(3*Z) - (512*T3*W2)/(9*Z) - 
      (16*h2*T*Z)/9 + (26*h*T2*Z)/3 - (56*T3*Z)/9 + (320*T2*W*Z)/9 - 
      (136*T2*Z2)/9)/v4;
 
coMttttZ = ((4864*T2*W2)/27 - (8192*W4)/81 - (20480*T2*W3)/(81*Z) - 
      (2240*T2*W*Z)/81 + (1664*T*W2*Z)/9 + (20480*W3*Z)/81 + 
      (8192*T2*W4)/(81*Z2) - (874*T2*Z2)/81 - (2080*T*W*Z2)/9 - 
      (8704*W2*Z2)/27 + (704*T*Z3)/9 + (16640*W*Z3)/81 - (4112*Z4)/81)/v4;
 
coMtWtW0 = (16*T3*W + 16*T2*W2 - (320*T*W3)/3 + (224*W4)/3 - 
      (64*T3*W2)/(3*Z) - (128*T2*W3)/(3*Z) + (64*T*W4)/(3*Z) + 
      (128*W5)/(3*Z) - (8*T3*Z)/3 + (152*T2*W*Z)/3 - (80*T*W2*Z)/3 - 
      (64*W3*Z)/3)/v4;
 
coMWWWW0 = ((-1504*W4)/3 - (64*W5)/(3*Z) + 304*W3*Z + (256*W6)/Z2 - 
      (104*W2*Z2)/3 - (8*W*Z3)/3)/v4;
 
coMWWWWh = ((8*h3*W)/3 - (52*h2*W2)/3 + 16*h*W3 + (272*W4)/3 - 
      (8*h3*W2)/(3*Z) + (16*h2*W3)/Z - (160*h*W4)/(3*Z) + (64*W5)/Z - 
      (2*h3*Z)/3 + (16*h2*W*Z)/3 - (40*h*W2*Z)/3 - (64*W3*Z)/3 - (h2*Z2)/3 - 
      (4*h*W*Z2)/3 - (4*W2*Z2)/3)/v4;
 
coMWWWWZ = ((2000*W4)/3 - (1408*W5)/(3*Z) - 416*W3*Z - (256*W6)/Z2 - 
      (8*W2*Z2)/3 + (64*W*Z3)/3 - 3*Z4)/v4;
 
coS00h = ((103*h2)/36 + (200*h*W)/27 + (1600*W2)/27 - 
      (h*(103*h2 - 800*h*W + 2560*W2))/(36*(h - 4*Z)) - (640*h*W2)/(27*Z) - 
      (103*h*Z)/27 - (2000*W*Z)/27 + (1030*Z2)/27)/v4;
 
coS00W = ((2872*W2)/9 - (3200*W3)/(9*Z) - (560*W*Z)/9 - (2368*W4)/(3*Z2))/v4;
 
coS0tW = ((-698*T2)/9 + (604*T*W)/9 + (136*W2)/3 + (128*W3)/T + 
      (32*T3)/(3*Z) + (1312*T2*W)/(9*Z) - (320*T*W2)/(9*Z) - (64*W3)/Z - 
      (256*W4)/(3*T*Z) + (55*T*Z)/3 - (8*T2*Z)/W + (137*W*Z)/9 - 
      (176*W2*Z)/(3*T) - (8*T2*(2*T + Z)*(2*T + Z))/(3*(T - W)*Z) - 
      ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(12*(4*T - Z)*Z) + 
      (928*T2*W2)/(9*Z2) + (928*T*W3)/(9*Z2) - (1856*W4)/(9*Z2) - (25*Z2)/12)/v4;
 
coShhZ = ((11*h2)/12 - (11*h3)/(4*(h - 4*Z)) - (2*h3)/(3*Z) + (71*h*Z)/9 + 
      (10*Z2)/3)/v4;
 
coShtt = ((17*h2)/36 - (157*h*T)/9 + (130*T2)/9 + (760*h*W)/27 - 
      (2800*T*W)/27 - (256*W2)/27 - (17*h3 + 28*h2*T - 160*h2*W - 
        1280*h*T*W + 512*h*W2 + 4096*T*W2)/(36*(h - 4*Z)) + 
      (640*h*T*W)/(9*Z) - (1280*T2*W)/(9*Z) - (704*h*W2)/(27*Z) + 
      (1472*T*W2)/(27*Z) - (404*h*Z)/27 + (1406*T*Z)/27 + (320*W*Z)/27 - 
      ((-h + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(3*(4*T - Z)) - (512*h*T*W2)/(9*Z2) + 
      (1024*T2*W2)/(9*Z2) - (55*Z2)/27)/v4;
 
coShWW = ((107*h2)/8 - (320*h*W)/9 + (460*W2)/9 - 
      ((h - 16*W)*(h2 + 80*h*W + 192*W2))/(24*(h - 4*Z)) - (80*h2*W)/(3*Z) + 
      (712*h*W2)/(9*Z) - (1408*W3)/(9*Z) - (49*h*Z)/18 + (4*h2*Z)/(3*W) + 
      (16*W*Z)/9 - (32*h2*W2)/(3*Z2) + (224*h*W3)/(3*Z2) - (128*W4)/Z2 + 
      (5*Z2)/9)/v4;
 
coSttZ = (-(h*T) + 14*T2 + (3040*T*W)/81 - (1792*W2)/9 + 
      (10240*W3)/(27*T) - (3*h*(h - 4*T)*T)/(h - 4*Z) - (9344*T*W2)/(27*Z) + 
      (5120*W3)/(27*Z) - (4096*W4)/(27*T*Z) + (1712*T*Z)/81 - (48*T2*Z)/h + 
      (2720*W*Z)/27 - (3328*W2*Z)/(9*T) + (-8*W + 5*Z)*(-8*W + 5*Z)*(-8*W + 5*Z)*(-8*W + 5*Z)/(9*(4*T - Z)*Z) + 
      (40960*T*W3)/(81*Z2) - (2048*W4)/(27*Z2) - (281*Z2)/27 + 
      (4480*W*Z2)/(27*T) - (16384*T*W4)/(81*Z3) - (856*Z3)/(27*T))/v4;
 
coSWWZ = (-h2/2 - 2*h*W - (10598*W2)/9 - (h*(h2 - 4*h*W + 12*W2))/
       (2*(h - 4*Z)) + (4672*W3)/(9*Z) + (14*h*Z)/3 + 104*W*Z + (24*W2*Z)/h + 
      (11776*W4)/(9*Z2) + (64*Z2)/3 + (512*W5)/Z3 + (4*Z3)/(3*W))/v4;
 
coSZZZ = ((-31*h2)/48 - (3*h3)/(16*(h - 4*Z)) + (25*h*Z)/12 - (23*Z2)/3 + 
      (12*Z3)/h)/v4;
 
coTbar0tt = ((-3712*T*W2)/27 + (1024*W3)/9 + (2048*T*W3)/(9*Z) - 
      (4096*W4)/(81*Z) + (896*T*W*Z)/81 - (2432*W2*Z)/27 - 
      (8192*T*W4)/(81*Z2) + (2176*W*Z2)/81)/v4;
 
coTbar0WW = (-448*W3 + (320*W4)/(3*Z) + 80*W2*Z + (256*W5)/Z2 + 
      (16*W*Z2)/3)/v4;
 
coTh00 = ((103*h3)/48 + (350*h2*W)/27 + (2080*h*W2)/27 - 
      (h2*(103*h2 - 800*h*W + 2560*W2))/(48*(h - 4*Z)) - (640*h2*W2)/(27*Z) - 
      (721*h2*Z)/108 - (2600*h*W*Z)/27 - (640*W2*Z)/9 + (1339*h*Z2)/27 + 
      (800*W*Z2)/9 - (412*Z3)/9)/v4;
 
coThhZ = (3*h3 - (3*h4)/(h - 4*Z) - (4*h4)/(3*Z) + (56*h2*Z)/9 + 
      (64*h*Z2)/9 - 16*Z3)/v4;
 
coThtt = ((17*h3)/48 - (169*h2*T)/12 + (56*h*T2)/9 + (430*h2*W)/27 - 
      (1040*h*T*W)/9 + (640*T2*W)/9 + (704*h*W2)/27 - (512*T*W2)/9 - 
      (h*(17*h3 + 28*h2*T - 160*h2*W - 1280*h*T*W + 512*h*W2 + 4096*T*W2))/
       (48*(h - 4*Z)) + (160*h2*T*W)/(3*Z) - (640*h*T2*W)/(9*Z) - 
      (416*h2*W2)/(27*Z) + (640*h*T*W2)/(9*Z) - (512*T2*W2)/(9*Z) - 
      (893*h2*Z)/108 + (289*h*T*Z)/9 - (56*T2*Z)/9 - (880*h*W*Z)/27 + 
      (640*T*W*Z)/9 - (128*W2*Z)/9 - (h*(h - Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/
       (6*(-4*T + Z)) - (128*h2*T*W2)/(3*Z2) + (512*h*T2*W2)/(9*Z2) + 
      (829*h*Z2)/54 - (164*T*Z2)/9 + (160*W*Z2)/9 - (68*Z3)/9)/v4;
 
coThWW = ((169*h3)/32 + (55*h2*W)/9 - (550*h*W2)/9 + (368*W3)/3 - 
      (h*(h - 16*W)*(h2 + 80*h*W + 192*W2))/(32*(h - 4*Z)) - h4/(4*Z) - 
      (43*h3*W)/(3*Z) + (484*h2*W2)/(9*Z) - (144*h*W3)/Z + (64*W4)/Z - 
      (391*h2*Z)/72 + (2*h3*Z)/(3*W) + (176*h*W*Z)/9 + 24*W2*Z - 
      (16*h3*W2)/(3*Z2) + (160*h2*W3)/(3*Z2) - (64*h*W4)/Z2 + (19*h*Z2)/18 - 
      (2*h2*Z2)/(3*W) - 12*W*Z2 + (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/
       (12*(h - 4*W)*Z) - (2*Z3)/3)/v4;
 
coTt0W = ((-566*T3)/9 - (494*T2*W)/3 - (836*T*W2)/9 + 40*W3 + 
      (16*T4)/(3*Z) + (656*T3*W)/(9*Z) - (1136*T2*W2)/(9*Z) - (256*T*W3)/Z - 
      (64*W4)/(3*Z) + (527*T2*Z)/18 - (4*T3*Z)/W + (19*T*W*Z)/18 - 
      (77*W2*Z)/3 - (4*T3*(2*T + Z)*(2*T + Z))/(3*(T - W)*(T - W)) - 
      ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(32*(4*T - Z)) - 
      (4*T2*(2*T + Z)*(6*T2 - 29*T*Z - 10*Z2))/(9*(T - W)*Z) + 
      (464*T3*W2)/(9*Z2) + (464*T2*W3)/(9*Z2) - (928*T*W4)/(9*Z2) - 
      (481*T*Z2)/72 + (4*T2*Z2)/W - (5*W*Z2)/8 - (25*Z3)/32)/v4;
 
coTtht = ((-83*h2*T)/36 - (175*h*T2)/9 + (184*T3)/9 + (20*h2*W)/9 + 
      (40*h*T*W)/9 - (1600*T2*W)/27 + (160*h*W2)/9 - (2240*T*W2)/27 - 
      (T*(17*h3 + 28*h2*T - 160*h2*W - 1280*h*T*W + 512*h*W2 + 4096*T*W2))/
       (36*(h - 4*Z)) + (80*h2*T*W)/(9*Z) + (640*h*T2*W)/(9*Z) - 
      (1280*T3*W)/(9*Z) - (16*h2*W2)/(9*Z) - (64*h*T*W2)/(9*Z) + 
      (512*T2*W2)/(27*Z) - (25*h2*Z)/36 - (107*h*T*Z)/9 + (1274*T2*Z)/27 - 
      (200*h*W*Z)/9 + (2800*T*W*Z)/27 - 32*W2*Z - (64*h2*T*W2)/(9*Z2) - 
      (512*h*T2*W2)/(9*Z2) + (1024*T3*W2)/(9*Z2) + (125*h*Z2)/18 - 
      (704*T*Z2)/27 + 40*W*Z2 - ((-8*W + 5*Z)*(-8*W + 5*Z)*(h2 - 10*h*Z + 18*Z2))/
       (36*(4*T - Z)) - (25*Z3)/2)/v4;
 
coTttZ = (3*h*T2 + 12*T3 - (1280*T2*W)/81 - (640*T*W2)/27 - (2560*W3)/27 - 
      (3*h*(h - 4*T)*T2)/(h - 4*Z) - (8192*T2*W2)/(27*Z) - 
      (10240*T*W3)/(81*Z) + (1024*W4)/(27*Z) + (1280*T2*Z)/81 - (48*T3*Z)/h + 
      (10400*T*W*Z)/81 + (544*W2*Z)/9 + (-8*W + 5*Z)*(-8*W + 5*Z)*(-8*W + 5*Z)*(-8*W + 5*Z)/(12*(4*T - Z)) + 
      (40960*T2*W3)/(81*Z2) + (4096*T*W4)/(81*Z2) - (2741*T*Z2)/81 + 
      (48*T2*Z2)/h - (40*W*Z2)/27 - (16384*T2*W4)/(81*Z3) - (1223*Z3)/108)/v4;
 
coTW00 = (448*W3 - (32*W4)/(3*Z) - (304*W2*Z)/3 - (1184*W5)/(3*Z2) + 
      (176*W*Z2)/3)/v4;
 
coTW0t = ((256*T3)/9 - 80*T2*W + (452*T*W2)/9 + (232*W3)/3 + 
      (320*W4)/(3*T) + (16*T4)/(3*Z) + (16*T3*W)/(3*Z) + (512*T2*W2)/(9*Z) - 
      (752*T*W3)/(9*Z) + (64*W4)/(3*Z) - (128*W5)/(3*T*Z) + (100*T2*Z)/9 + 
      (187*T*W*Z)/6 - (55*W2*Z)/18 - (280*W3*Z)/(3*T) + 
      (4*T3*(2*T + Z)*(2*T + Z))/(3*(T - W)*(T - W)) - (W*(-4*W + Z)*(8*W + Z)*
        (-8*W + 5*Z)*(-8*W + 5*Z))/(24*(4*T - Z)*Z) + (464*T2*W3)/(9*Z2) + 
      (464*T*W4)/(9*Z2) - (928*W5)/(9*Z2) + (28*T*Z2)/9 + (77*W*Z2)/8 + 
      (88*W2*Z2)/(3*T) - (4*T2*(2*T + Z)*(6*T2 + 35*T*Z + 10*Z2))/
       (9*(T - W)*Z))/v4;
 
coTWhW = ((17*h3)/24 + (199*h2*W)/8 - 90*h*W2 + 144*W3 - 
      ((h - 16*W)*W*(h2 + 80*h*W + 192*W2))/(24*(h - 4*Z)) + h4/(8*Z) + 
      (h3*W)/(2*Z) - (58*h2*W2)/(3*Z) + (112*h*W3)/(3*Z) - (448*W4)/(9*Z) + 
      (4*h2*Z)/3 - (19*h*W*Z)/6 + (184*W2*Z)/3 - (8*h2*W3)/(3*Z2) + 
      (224*h*W4)/(3*Z2) - (128*W5)/Z2 + (h*Z2)/2 + (40*W*Z2)/9 - 
      (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/(24*(h - 4*W)*Z))/v4;
 
coTWWZ = (-h3/8 - (20*h*W2)/3 - (5894*W3)/3 - (h*W*(h2 - 4*h*W + 12*W2))/
       (2*(h - 4*Z)) - (384*W4)/Z + (h2*Z)/2 + (20*h*W*Z)/3 + (2336*W2*Z)/9 + 
      (24*W3*Z)/h + (13696*W5)/(9*Z2) - (3*h*Z2)/2 + (158*W*Z2)/3 - 
      (24*W2*Z2)/h + (h2*(h2 - 4*h*Z + 12*Z2))/(8*(h - 4*W)) + (512*W6)/Z3 + 
      (8*Z3)/3)/v4;
 
coTZ00 = ((-417280*W3)/81 + (148480*W4)/(81*Z) + (177920*W2*Z)/27 - 
      (323200*W*Z2)/81 + (80920*Z3)/81)/v4;
 
coUhZ00 = ((103*h3)/288 - (25*h2*W)/9 - (80*h*W2)/3 - 
      (h2*(103*h2 - 800*h*W + 2560*W2))/(288*(h - 4*Z)) + (103*h2*Z)/72 + 
      (100*h*W*Z)/3 + (1600*W2*Z)/9 - (103*h*Z2)/6 - (2000*W*Z2)/9 + 
      (1030*Z3)/9)/v4;
 
coUhZhZ = ((17*h3)/8 - h4/(8*(h - 4*Z)) - h4/(3*Z) - (29*h2*Z)/6 - 
      (2*h*Z2)/3 + 8*Z3)/v4;
 
coUhZtt = ((17*h3)/288 - (49*h2*T)/72 + (5*h2*W)/3 - (200*h*T*W)/9 + 
      (16*h*W2)/9 + (128*T*W2)/3 - (h*(17*h3 + 28*h2*T - 160*h2*W - 
         1280*h*T*W + 512*h*W2 + 4096*T*W2))/(288*(h - 4*Z)) + 
      (80*h2*T*W)/(9*Z) - (16*h2*W2)/(9*Z) + (128*h*T*W2)/(9*Z) - 
      (11*h2*Z)/24 + (107*h*T*Z)/18 - (20*h*W*Z)/9 - (160*T*W*Z)/3 + 
      (128*W2*Z)/9 - (64*h2*T*W2)/(9*Z2) - (h*Z2)/18 + (14*T*Z2)/3 - 
      (160*W*Z2)/9 - ((-8*W + 5*Z)*(-8*W + 5*Z)*(h2 - 4*h*Z + 12*Z2))/(36*(4*T - Z)) + 
      (95*Z3)/9)/v4;
 
coUhZWW = (h3/192 - (5*h2*W)/3 - (13*h*W2)/3 - 48*W3 - 
      (h*(h - 16*W)*(h2 + 80*h*W + 192*W2))/(192*(h - 4*Z)) + 
      (8*h2*W2)/(3*Z) - (16*h*W3)/Z + (h2*Z)/48 + (28*h*W*Z)/3 - 60*W2*Z + 
      (8*h2*W3)/Z2 - (35*h*Z2)/12 + (88*W*Z2)/3 + (5*Z3)/3)/v4;
 
coUtt0W = ((-352*T2*W)/9 - (992*T*W2)/9 - (176*W3)/9 + (128*W4)/(3*T) + 
      (320*T2*W2)/(9*Z) + (64*T*W3)/Z - (128*W4)/(9*Z) + (77*T2*Z)/9 + 
      (241*T*W*Z)/9 - (122*W2*Z)/9 - (64*W3*Z)/T + 
      ((-4*W + Z)*(8*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/(48*(4*T - Z)) + (331*T*Z2)/36 + 
      (23*W*Z2)/4 + (88*W2*Z2)/(3*T) + (25*Z3)/48)/v4;
 
coUttht = (8*h*T2 + (80*h*T*W)/9 - (640*T2*W)/9 + 16*h*W2 + (320*T*W2)/9 - 
      (64*h*T*W2)/(9*Z) + (512*T2*W2)/(9*Z) - (97*h*T*Z)/9 + (56*T2*Z)/9 - 
      20*h*W*Z - (400*T*W*Z)/9 + (16*W2*Z)/3 + (Z*(-h + Z)*(-8*W + 5*Z)*(-8*W + 5*Z))/
       (12*(4*T - Z)) + (37*h*Z2)/4 + (143*T*Z2)/9 - (20*W*Z2)/3 + 
      (25*Z3)/12)/v4;
 
coUtttZ = (-8*h*T2 - (6400*T*W2)/27 + (25600*W3)/81 + (2048*W4)/(27*T) + 
      (20480*T*W3)/(81*Z) - (10240*W4)/(81*Z) + 16*T2*Z + (8000*T*W*Z)/81 - 
      (10688*W2*Z)/27 - (5120*W3*Z)/(27*T) - (-8*W + 5*Z)*(-8*W + 5*Z)*(-8*W + 5*Z)*(-8*W + 5*Z)/(18*(4*T - Z)) - 
      (8192*T*W4)/(81*Z2) - (764*T*Z2)/81 + (20080*W*Z2)/81 + 
      (1664*W2*Z2)/(9*T) - (10721*Z3)/162 - (2240*W*Z3)/(27*T) + 
      (428*Z4)/(27*T))/v4;
 
coUWW00 = (-384*W3 - (480*W4)/Z + 72*W2*Z)/v4;
 
coUWW0t = (-176*T2*W - (448*T*W2)/3 - 128*W3 + (16*T2*W2)/Z - (48*T*W3)/Z - 
      (160*W4)/Z + 52*T2*Z + (4*T*W*Z)/3 + 24*W2*Z + 4*T*Z2 + (4*T2*Z2)/W)/v4;
 
coUWWhW = ((-17*h3)/24 + (53*h2*W)/2 - (218*h*W2)/3 - (296*W3)/3 - 
      h4/(8*Z) - (h3*W)/(2*Z) - (14*h2*W2)/(3*Z) + (40*h*W3)/(3*Z) - 
      (64*W4)/Z - 8*h2*Z + (64*h*W*Z)/3 + 8*W2*Z + (13*h*Z2)/6 - 
      (2*h2*Z2)/(3*W) + (2*W*Z2)/3 + (h2*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/
       (24*(h - 4*W)*Z))/v4;
 
coUWWWZ = ((8*h2*W)/3 - (16*h*W2)/3 - (2048*W3)/3 + (2368*W4)/(3*Z) - 
      (16*h*W*Z)/3 + (1664*W2*Z)/3 + (256*W5)/Z2 - (88*W*Z2)/3 - (52*Z3)/3 - 
      (2*Z4)/(3*W))/v4;
 
coUZhhh = ((15*h3)/16 + (9*h4)/(16*(h - 4*Z)) - (37*h2*Z)/4 + 9*h*Z2)/v4;
 
coUZhtt = ((13*h2*T)/4 - 17*h*T2 - (1024*T*W2)/9 + 
      (3*h2*(h - 4*T)*T)/(4*(h - 4*Z)) - 15*h*T*Z + 44*T2*Z + 
      (1280*T*W*Z)/9 - (220*T*Z2)/9)/v4;
 
coUZhWW = ((3*h3)/4 - (2*h2*W)/3 + (179*h*W2)/6 + 128*W3 + 
      (h2*(h2 - 4*h*W + 12*W2))/(8*(h - 4*Z)) - 2*h2*Z - 12*h*W*Z - 
      (194*W2*Z)/3 + (35*h*Z2)/6 - (82*W*Z2)/3 + (h2*(h2 - 4*h*Z + 12*Z2))/
       (8*(h - 4*W)))/v4;
 
coUZhZZ = ((103*h3)/96 + (3*h4)/(32*(h - 4*Z)) + (7*h2*Z)/24 + 
      (71*h*Z2)/6 - (82*Z3)/3)/v4;
 
coAhAh = (h/12 - h2/(12*(h - 4*Z)) - (2*h2)/(9*Z) - (5*Z)/9)/v4;
 
coAhAt = ((-92*h)/27 + 6*T - (880*W)/27 - (17*h2 - 160*h*W + 512*W2)/
       (9*(h - 4*Z)) + (400*h*W)/(27*Z) + (320*W2)/(27*Z) + (374*Z)/27 - 
      (12*T*Z)/h - (h*(-8*W + 5*Z)*(-8*W + 5*Z))/(9*(4*T - Z)*Z) - (320*h*W2)/(27*Z2))/v4;
 
coAhAW = ((-91*h)/18 - (34*W)/9 + (h2 - 16*h*W + 192*W2)/(3*(h - 4*Z)) - 
      h2/Z - (52*h*W)/(9*Z) + (80*W2)/(9*Z) - (29*Z)/18 - (2*h*Z)/(3*W) + 
      (6*W*Z)/h + (h*(6*h2 - 5*h*Z - 13*Z2))/(6*(h - 4*W)*Z) + 
      (40*h*W2)/(3*Z2))/v4;
 
coAhAZ = ((139*h)/144 + (13*h2)/(16*(h - 4*Z)) + (7*h2)/(9*Z) - 
      (13*Z)/36 + (3*Z2)/h)/v4;
 
coAhB00 = ((103*h2)/144 + (250*h*W)/27 + (160*W2)/9 - 
      (h*(103*h2 - 800*h*W + 2560*W2))/(144*(h - 4*Z)) - (320*h*W2)/(27*Z) - 
      (515*h*Z)/108 - (200*W*Z)/9 - (640*W2*Z)/(9*h) + (103*Z2)/9 + 
      (800*W*Z2)/(9*h) - (412*Z3)/(9*h))/v4;
 
coAhBhZ = ((149*h2)/72 - (7*h3)/(24*(h - 4*Z)) - (7*h3)/(9*Z) - 
      (29*h*Z)/6 + (14*Z2)/3 - (16*Z3)/h)/v4;
 
coAhBtt = ((17*h2)/144 - (119*h*T)/108 + (110*h*W)/27 - (400*T*W)/9 + 
      (640*T2*W)/(9*h) + (224*W2)/9 - (512*T*W2)/(9*h) - 
      (17*h3 + 28*h2*T - 160*h2*W - 1280*h*T*W + 512*h*W2 + 4096*T*W2)/
       (144*(h - 4*Z)) + (400*h*T*W)/(27*Z) - (112*h*W2)/(27*Z) + 
      (256*T*W2)/(9*Z) - (512*T2*W2)/(9*h*Z) - (40*h*Z)/27 - (T*Z)/9 - 
      (56*T2*Z)/(9*h) - (280*W*Z)/9 + (640*T*W*Z)/(9*h) - (128*W2*Z)/(9*h) - 
      (h*(-8*W + 5*Z)*(-8*W + 5*Z))/(36*(4*T - Z)) - (320*h*T*W2)/(27*Z2) + (119*Z2)/9 - 
      (164*T*Z2)/(9*h) + (160*W*Z2)/(9*h) - (68*Z3)/(9*h))/v4;
 
coAhBWW = ((-45*h2)/32 + (161*h*W)/9 - (214*W2)/3 + (368*W3)/(3*h) - 
      ((h - 16*W)*(h2 + 80*h*W + 192*W2))/(96*(h - 4*Z)) - h3/(4*Z) - 
      (h2*W)/Z + (80*h*W2)/(9*Z) - (112*W3)/(3*Z) + (64*W4)/(h*Z) - 
      (437*h*Z)/72 + (56*W*Z)/3 + (24*W2*Z)/h + (40*h*W3)/(3*Z2) + Z2/2 - 
      (2*h*Z2)/(3*W) - (12*W*Z2)/h + (h*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/
       (12*(h - 4*W)*Z) - (2*Z3)/(3*h))/v4;
 
coAtAt = ((-25*h)/9 - 8*T + (3616*W)/27 - (2560*W2)/(27*T) + 
      (80*h*W)/(9*Z) - (14384*W2)/(27*Z) + (2048*W3)/(27*T*Z) + (76*Z)/81 + 
      (12*T*Z)/h + (4544*W*Z)/(81*T) - (64*h*W2)/(9*Z2) + 
      (57344*W3)/(81*Z2) - (2048*W4)/(81*T*Z2) - (428*Z2)/(27*T) - 
      (8192*W4)/(27*Z3) + (-18432*W4 - 192*h*W2*Z + 39424*W3*Z + 240*h*W*Z2 - 
        20784*W2*Z2 - 75*h*Z3 - 1968*W*Z3 + 2300*Z4)/(27*(4*T - Z)*Z2))/v4;
 
coAtAW = ((190*T)/3 - (1292*W)/27 - (64*W2)/T + (16*T2)/(3*Z) + 
      (64*T*W)/(3*Z) + (6064*W2)/(27*Z) + (128*W3)/(3*T*Z) + (587*Z)/27 + 
      (4*T*Z)/W + (88*W*Z)/(3*T) - (4*T*(2*T + Z)*(2*T + Z))/(3*(T - W)*Z) - 
      (25600*W3)/(27*Z2) + ((-8*W + 5*Z)*(-8*W + 5*Z)*(48*W2 - 40*W*Z + Z2))/
       (9*(4*T - Z)*Z2) + (2048*W4)/(3*Z3))/v4;
 
coAtAZ = ((116*h)/27 + 2*T + (800*W)/27 + (1664*W2)/(9*T) + 
      (2*(17*h2 - 160*h*W + 512*W2))/(9*(h - 4*Z)) - (640*h*W)/(27*Z) + 
      (128*W2)/(27*Z) - (5120*W3)/(27*T*Z) - (232*Z)/27 - (2240*W*Z)/(27*T) + 
      (512*h*W2)/(27*Z2) + (2048*W4)/(27*T*Z2) + (428*Z2)/(27*T) - 
      (2*(-8*W + 5*Z)*(-8*W + 5*Z)*(128*W2 - 3*h*Z - 160*W*Z + 74*Z2))/
       (27*(4*T - Z)*Z2))/v4;
 
coAtB00 = ((-208*T2)/9 - (64*T*W)/3 - (7936*W2)/9 - (32*T*W2)/(9*Z) + 
      (132544*W3)/(81*Z) - (56*T*Z)/3 - (3392*W*Z)/27 - 
      (4*T2*(2*T + Z)*(2*T + Z))/(3*(T - W)*(T - W)) + (8*T*(2*T + Z)*(16*T + 5*Z))/
       (9*(T - W)) - (20480*W4)/(27*Z2) + (12278*Z2)/81 - 
      (2*(-8*W + 5*Z)*(-8*W + 5*Z)*(160*W2 - 200*W*Z + 103*Z2))/(27*(4*T - Z)*Z))/v4;
 
coAtBhZ = ((-385*h2)/54 + 4*h*T - (1840*h*W)/27 - (256*W2)/3 - 
      (h*(17*h2 - 160*h*W + 512*W2))/(18*(h - 4*Z)) + (640*h2*W)/(27*Z) + 
      (1280*h*W2)/(27*Z) + (566*h*Z)/27 - 16*T*Z + (320*W*Z)/3 - 
      (512*h2*W2)/(27*Z2) - (64*Z2)/3 + (48*T*Z2)/h - 
      (2*(-8*W + 5*Z)*(-8*W + 5*Z)*(h2 - 4*h*Z + 12*Z2))/(9*(4*T - Z)*Z))/v4;
 
coAtBtt = ((2144*T*W)/9 - (4552*W2)/9 + (5120*W3)/(27*T) - 
      (9376*T*W2)/(9*Z) + (61696*W3)/(81*Z) - (2048*W4)/(27*T*Z) + 
      (200*T*Z)/81 + (2192*W*Z)/27 - (1664*W2*Z)/(9*T) + 
      (114688*T*W3)/(81*Z2) - (1024*W4)/(3*Z2) + (668*Z2)/81 + 
      (2240*W*Z2)/(27*T) - (16384*T*W4)/(27*Z3) - (428*Z3)/(27*T) + 
      (2*(-512*W4 + 384*W3*Z + 1716*W2*Z2 - 2032*W*Z3 + 525*Z4))/
       (27*(4*T - Z)*Z))/v4;
 
coAtBWW = (-160*T*W + (652*W2)/3 - (16*T*W2)/Z - (27184*W3)/(27*Z) + 
      52*T*Z + (64*W*Z)/27 - (9472*W4)/(27*Z2) + (61*Z2)/27 + (4*T*Z2)/W - 
      ((-4*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z)*(12*W2 + 20*W*Z + Z2))/(9*(4*T - Z)*Z2) + 
      (2048*W5)/(3*Z3))/v4;
 
coAWAW = ((17*h)/2 - (916*W)/9 + h2/Z + (28*h*W)/(3*Z) - (776*W2)/(3*Z) + 
      (206*Z)/9 + (2*h*Z)/(3*W) - (6*W*Z)/h - (h*(6*h2 - 5*h*Z - 13*Z2))/
       (6*(h - 4*W)*Z) + (8*h*W2)/Z2 + (3872*W3)/(9*Z2) + (2*Z2)/(3*W) + 
      (21*Z2)/(2*(4*W - Z)) - (384*W4)/Z3)/v4;
 
coAWAZ = ((-28*h)/9 - (604*W)/9 - (2*(h2 - 16*h*W + 192*W2))/
       (3*(h - 4*Z)) - (32*h*W)/(9*Z) + (1174*W2)/(9*Z) - (46*Z)/9 - 
      (64*h*W2)/(3*Z2) + (3328*W3)/(9*Z2) - (2*Z2)/(3*W))/v4;
 
coAWB00 = ((160*T2)/9 + 16*T*W + (5296*W2)/27 + (32*T*W2)/(9*Z) - 
      (23680*W3)/(27*Z) + (40*T*Z)/3 + (4960*W*Z)/27 + 
      (4*T2*(2*T + Z)*(2*T + Z))/(3*(T - W)*(T - W)) - (4*T*(2*T + Z)*(26*T + 7*Z))/
       (9*(T - W)) + (2560*W4)/(3*Z2) + (2480*Z2)/27)/v4;
 
coAWBhZ = ((22*h2)/9 - (218*h*W)/9 + 96*W2 + (h*(h2 - 16*h*W + 192*W2))/
       (6*(h - 4*Z)) + (32*h2*W)/(9*Z) - (160*h*W2)/(3*Z) - (76*h*Z)/9 + 
      104*W*Z + (64*h2*W2)/(3*Z2) + 18*Z2 - (24*W*Z2)/h + 
      (h*(h2 - 4*h*Z + 12*Z2))/(2*(h - 4*W)))/v4;
 
coAWBtt = ((-2672*T*W)/27 + (2468*W2)/27 + (128*W3)/(3*T) + 
      (8032*T*W2)/(27*Z) - (8768*W3)/(27*Z) + (53*T*Z)/27 - (166*W*Z)/27 - 
      (64*W2*Z)/T - (25600*T*W3)/(27*Z2) + (256*W4)/Z2 + (1483*Z2)/108 + 
      (88*W*Z2)/(3*T) + ((-8*W + 5*Z)*(-8*W + 5*Z)*(48*W2 - 40*W*Z + Z2))/
       (36*(4*T - Z)*Z) + (2048*T*W4)/(3*Z3))/v4;
 
coAWBWW = ((17*h2)/6 - (46*h*W)/3 - (6212*W2)/9 + h3/(2*Z) + (2*h2*W)/Z + 
      (16*h*W2)/(3*Z) + (256*W3)/(3*Z) + 6*h*Z + (673*W*Z)/9 + (32*W4)/Z2 + 
      (781*Z2)/36 + (2*h*Z2)/(3*W) - (h*(h - Z)*(3*h2 + 20*h*Z + 4*Z2))/
       (6*(h - 4*W)*Z) - (768*W5)/Z3 + (2*Z3)/(3*W) + (21*Z3)/(4*(4*W - Z)))/v4;
 
coAZAZ = ((-103*h)/144 - (17*h2)/(16*(h - 4*Z)) - (5*h2)/(9*Z) + 
      (127*Z)/36 - (3*Z2)/h)/v4;
 
coAZB00 = ((-103*h2)/72 - (100*h*W)/27 + (29440*W2)/9 + 
      (h*(103*h2 - 800*h*W + 2560*W2))/(72*(h - 4*Z)) + (320*h*W2)/(27*Z) - 
      (208640*W3)/(81*Z) + (103*h*Z)/54 - (159200*W*Z)/81 + 
      (74240*W4)/(81*Z2) + (39224*Z2)/81)/v4;
 
coAZBhZ = ((-59*h2)/48 + (43*h3)/(48*(h - 4*Z)) + (7*h3)/(9*Z) - 
      (25*h*Z)/36 + (5*Z2)/3 - (12*Z3)/h)/v4;
 
coAZBtt = ((-17*h2)/72 + (49*h*T)/54 - (80*h*W)/27 + (160*T*W)/81 - 
      (128*W2)/27 - (5120*W3)/(27*T) + (17*h3 + 28*h2*T - 160*h2*W - 
        1280*h*T*W + 512*h*W2 + 4096*T*W2)/(72*(h - 4*Z)) - 
      (400*h*T*W)/(27*Z) + (112*h*W2)/(27*Z) + (512*T*W2)/(3*Z) + 
      (2048*W4)/(27*T*Z) + (109*h*Z)/108 - (436*T*Z)/81 + (160*W*Z)/27 + 
      (1664*W2*Z)/(9*T) + (320*h*T*W2)/(27*Z2) - (20480*T*W3)/(81*Z2) - 
      (176*Z2)/27 - (2240*W*Z2)/(27*T) - 
      ((-8*W + 5*Z)*(-8*W + 5*Z)*(512*W2 - 3*h*Z - 640*W*Z + 224*Z2))/
       (108*(4*T - Z)*Z) + (8192*T*W4)/(81*Z3) + (428*Z3)/(27*T))/v4;
 
coAZBWW = (-h2/48 + (22*h*W)/9 + (4960*W2)/9 + 
      ((h - 16*W)*(h2 + 80*h*W + 192*W2))/(48*(h - 4*Z)) - (92*h*W2)/(9*Z) - 
      (192*W3)/Z - (47*h*Z)/36 - (620*W*Z)/9 - (40*h*W3)/(3*Z2) - 
      (1408*W4)/(3*Z2) - (128*Z2)/9 - (256*W5)/Z3 - (2*Z3)/(3*W))/v4;
 
coB00B00 = ((-136544*W3)/81 + (34816*W4)/(81*Z) + (84704*W2*Z)/27 - 
      (203224*W*Z2)/81 + (62666*Z3)/81)/v4;
 
coB00BhZ = ((103*h3)/288 + (325*h2*W)/27 + (560*h*W2)/27 - 
      (h2*(103*h2 - 800*h*W + 2560*W2))/(288*(h - 4*Z)) - 
      (320*h2*W2)/(27*Z) - (1339*h2*Z)/216 - (700*h*W*Z)/27 + (320*W2*Z)/9 + 
      (721*h*Z2)/54 - (400*W*Z2)/9 + (206*Z3)/9)/v4;
 
coB00Btt = ((-32*T2*W)/9 - (9568*T*W2)/9 + (26144*W3)/81 + 
      (64*T2*W2)/(9*Z) + (133696*T*W3)/(81*Z) - (5504*W4)/(27*Z) - 
      (32*T2*Z)/9 + (4240*T*W*Z)/27 + (176*W2*Z)/3 - (20480*T*W4)/(27*Z2) + 
      (1298*T*Z2)/81 - (892*W*Z2)/3 - 
      ((-8*W + 5*Z)*(-8*W + 5*Z)*(160*W2 - 200*W*Z + 103*Z2))/(54*(4*T - Z)) + 
      (20291*Z3)/162)/v4;
 
coB00BWW = (-16*T*W2 - (5296*W3)/27 - (32*T2*W2)/(3*Z) - 
      (32*T*W3)/(3*Z) - (21760*W4)/(27*Z) + (8*T2*Z)/3 - (8*T*W*Z)/3 + 
      (80*W2*Z)/3 + (2560*W5)/(3*Z2) + (4*T*Z2)/3 + (2884*W*Z2)/27 + 
      (412*Z3)/27)/v4;
 
coBhZBhZ = ((217*h3)/144 - h4/(16*(h - 4*Z)) - (2*h4)/(9*Z) - 
      (167*h2*Z)/36 + 5*h*Z2 - 4*Z3)/v4;
 
coBhZBtt = ((17*h3)/288 - (259*h2*T)/216 + (125*h2*W)/27 - 
      (1240*h*T*W)/27 + (304*h*W2)/27 - (128*T*W2)/9 - 
      (h*(17*h3 + 28*h2*T - 160*h2*W - 1280*h*T*W + 512*h*W2 + 4096*T*W2))/
       (288*(h - 4*Z)) + (400*h2*T*W)/(27*Z) - (112*h2*W2)/(27*Z) + 
      (896*h*T*W2)/(27*Z) - (371*h2*Z)/216 + (433*h*T*Z)/54 - 
      (380*h*W*Z)/27 + (160*T*W*Z)/9 - (128*W2*Z)/9 - (320*h2*T*W2)/(27*Z2) + 
      (269*h*Z2)/54 - (14*T*Z2)/9 + (160*W*Z2)/9 - 
      ((-8*W + 5*Z)*(-8*W + 5*Z)*(h2 - 4*h*Z + 12*Z2))/(36*(4*T - Z)) - (41*Z3)/9)/v4;
 
coBhZBWW = (h3/192 - (31*h2*W)/9 - (311*h*W2)/9 + 16*W3 - 
      (h*(h - 16*W)*(h2 + 80*h*W + 192*W2))/(192*(h - 4*Z)) + 
      (92*h2*W2)/(9*Z) - (112*h*W3)/(3*Z) - (13*h2*Z)/144 + (148*h*W*Z)/9 + 
      (92*W2*Z)/3 + (40*h2*W3)/(3*Z2) - (89*h*Z2)/36 + 8*W*Z2 + Z3/3)/v4;
 
coBttBtt = ((-80*h*T*W)/9 + (1232*T2*W)/9 - (9388*T*W2)/27 + (352*W3)/9 + 
      (64*h*T*W2)/(9*Z) - (4816*T2*W2)/(9*Z) + (35200*T*W3)/(81*Z) - 
      (640*W4)/(27*Z) + (16*h*T*Z)/9 - (674*T2*Z)/81 + (9428*T*W*Z)/81 + 
      (607*W2*Z)/9 + (57344*T2*W3)/(81*Z2) - (15872*T*W4)/(81*Z2) - 
      (2561*T*Z2)/162 - (2897*W*Z2)/27 - (8192*T2*W4)/(27*Z3) + 
      (2687*Z3)/72 - (13312*W4 - 33024*W3*Z + 29976*W2*Z2 - 11896*W*Z3 + 
        1875*Z4)/(216*(4*T - Z)))/v4;
 
coBttBWW = (16*T2*W + (1556*T*W2)/3 - (5180*W3)/27 - (64*T2*W2)/(3*Z) - 
      (23008*T*W3)/(27*Z) - (3968*W4)/(27*Z) - (8*T2*Z)/3 - (764*T*W*Z)/27 + 
      97*W2*Z - (9472*T*W4)/(27*Z2) + (256*W5)/Z2 - (137*T*Z2)/27 + 
      (236*W*Z2)/27 - ((-4*W + Z)*(-8*W + 5*Z)*(-8*W + 5*Z)*(12*W2 + 20*W*Z + Z2))/
       (36*(4*T - Z)*Z) + (2048*T*W5)/(3*Z3) + (197*Z3)/108)/v4;
 
coBWWBWW = ((-8*h2*W)/3 + (32*h*W2)/3 + (4058*W3)/9 + (8*h2*W2)/(3*Z) - 
      (32*h*W3)/(3*Z) + (248*W4)/(3*Z) + (2*h2*Z)/3 - (8*h*W*Z)/3 - 
      (51*W2*Z)/2 - (448*W5)/Z2 - (145*W*Z2)/24 - (384*W6)/Z3 + 
      (797*Z3)/288 + (21*Z4)/(32*(4*W - Z)))/v4; 

result = co1 + 
 Ah*coAh + 
 At*coAt + 
 AW*coAW + 
 AZ*coAZ + 
 B00*coB00 + 
 BhZ*coBhZ +  
 Btt*coBtt +  
 BWW*coBWW + 
 Ah*Ah*coAhAh + 
 Ah*At*coAhAt + 
 Ah*AW*coAhAW + 
 Ah*AZ*coAhAZ + 
 Ah*B00*coAhB00 + 
 Ah*BhZ*coAhBhZ + 
 Ah*Btt*coAhBtt + 
 Ah*BWW*coAhBWW + 
 At*At*coAtAt + 
 At*AW*coAtAW + 
 At*AZ*coAtAZ + 
 At*B00*coAtB00 + 
 At*BhZ*coAtBhZ + 
 At*Btt*coAtBtt + 
 At*BWW*coAtBWW + 
 AW*AW*coAWAW + 
 AW*AZ*coAWAZ + 
 AW*B00*coAWB00 + 
 AW*BhZ*coAWBhZ + 
 AW*Btt*coAWBtt + 
 AW*BWW*coAWBWW + 
 AZ*AZ*coAZAZ + 
 AZ*B00*coAZB00 + 
 AZ*BhZ*coAZBhZ + 
 AZ*Btt*coAZBtt + 
 AZ*BWW*coAZBWW + 
 B00*B00*coB00B00 + 
 B00*BhZ*coB00BhZ + 
 B00*Btt*coB00Btt + 
 B00*BWW*coB00BWW + 
 BhZ*BhZ*coBhZBhZ + 
 BhZ*Btt*coBhZBtt + 
 BhZ*BWW*coBhZBWW + 
 Btt*Btt*coBttBtt + 
 Btt*BWW*coBttBWW + 
 BWW*BWW*coBWWBWW +
 I00h*coI00h + 
 I00t*coI00t + 
 I00W*coI00W + 
 I00Z*coI00Z + 
 I0hW*coI0hW + 
 I0hZ*coI0hZ + 
 I0tW*coI0tW + 
 I0WZ*coI0WZ + 
 Ihhh*coIhhh + 
 Ihtt*coIhtt + 
 IhWW*coIhWW + 
 IhZZ*coIhZZ + 
 IttZ*coIttZ + 
 IWWZ*coIWWZ + 
 M00000Z*coM00000 + 
 M0000W*coM0000W + 
 M0000Z*coM0000Z + 
 M0t0tW*coM0t0tW + 
 M0W0W0*coM0W0W0 + 
 M0W0Wt*coM0W0Wt + 
 MhhZZh*coMhhZZh + 
 MhtZtt*coMhtZtt + 
 MhWZWW*coMhWZWW + 
 MhZZhZ*coMhZZhZ + 
 Mtttt0*coMtttt0 + 
 Mtttth*coMtttth + 
 MttttZ*coMttttZ + 
 MtWtW0*coMtWtW0 + 
 MWWWW0*coMWWWW0 + 
 MWWWWh*coMWWWWh + 
 MWWWWZ*coMWWWWZ + 
 S00h*coS00h + 
 S00W*coS00W + 
 S0tW*coS0tW + 
 ShhZ*coShhZ + 
 Shtt*coShtt + 
 ShWW*coShWW + 
 SttZ*coSttZ + 
 SWWZ*coSWWZ + 
 SZZZ*coSZZZ + 
 Tbar0tt*coTbar0tt + 
 Tbar0WW*coTbar0WW + 
 Th00*coTh00 + 
 ThhZ*coThhZ + 
 Thtt*coThtt + 
 ThWW*coThWW + 
 Tt0W*coTt0W + 
 Ttht*coTtht + 
 TttZ*coTttZ + 
 TW00*coTW00 + 
 TW0t*coTW0t + 
 TWhW*coTWhW + 
 TWWZ*coTWWZ + 
 TZ00*coTZ00 + 
 UhZ00*coUhZ00 + 
 UhZhZ*coUhZhZ + 
 UhZtt*coUhZtt + 
 UhZWW*coUhZWW + 
 Utt0W*coUtt0W + 
 Uttht*coUttht + 
 UtttZ*coUtttZ + 
 UWW00*coUWW00 + 
 UWW0t*coUWW0t + 
 UWWhW*coUWWhW + 
 UWWWZ*coUWWWZ + 
 UZhhh*coUZhhh + 
 UZhtt*coUZhtt + 
 UZhWW*coUZhWW + 
 UZhZZ*coUZhZZ; 

  return result;
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Performs all needed basis integral evaluations.                    */

int SMDR_MZ_DoTSILZ (float loopOrder)
{
  TSIL_DATA bar;
  int success = 1;

  /* Do all 1-loop functions: */
  Ah = TSIL_A (h, Q2);
  At = TSIL_A (T, Q2);
  AW = TSIL_A (W, Q2);
  AZ = TSIL_A (Z, Q2);
  Btt = TSIL_B (T, T, Z, Q2);
  B00 = TSIL_B (0, 0, Z, Q2);
  BWW = TSIL_B (W, W, Z, Q2);
  BhZ = TSIL_B (h, Z, Z, Q2);

  Ab = TSIL_A (b, Q2);
  Bbb = TSIL_B (b, b, Z, Q2);

  Atau = TSIL_A (tau, Q2);
  Btautau = TSIL_B (tau, tau, Z, Q2);

  /* Then do 2-loop evals if needed... */

  /* These are needed for the QCD contribution: */
  if (loopOrder > 1.1) {
    success *= TSIL_Manalytic (T, T, T, T, 0, Z, &Mtttt0);
    success *= TSIL_Manalytic (0, 0, 0, 0, 0, Z, &M00000Z);
    success *= TSIL_Tbaranalytic (0, T, T, Z, Q2, &Tbar0tt);
  }

  /* These are needed for nonQCD part: */
  if (loopOrder > 1.6) {
    I00h = TSIL_I2 (0, 0, h, Q2);
    I00t = TSIL_I2 (0, 0, T, Q2);
    I00W = TSIL_I2 (0, 0, W, Q2);
    I00Z = TSIL_I2 (0, 0, Z, Q2);
    I0hW = TSIL_I2 (0, h, W, Q2);
    I0hZ = TSIL_I2 (0, h, Z, Q2);
    I0tW = TSIL_I2 (0, T, W, Q2);
    I0WZ = TSIL_I2 (0, W, Z, Q2);
    Ihhh = TSIL_I2 (h, h, h, Q2);
    Ihtt = TSIL_I2 (h, T, T, Q2);
    IhWW = TSIL_I2 (h, W, W, Q2);
    IhZZ = TSIL_I2 (h, Z, Z, Q2);
    IttZ = TSIL_I2 (T, T, Z, Q2);
    IWWZ = TSIL_I2 (W, W, Z, Q2);

    success *= TSIL_Manalytic (0, 0, 0, 0, W, Z, &M0000W);
    success *= TSIL_Manalytic (0, 0, 0, 0, Z, Z, &M0000Z);
    success *= TSIL_Manalytic (0, W, 0, W, 0, Z, &M0W0W0);
    success *= TSIL_Manalytic (W, W, W, W, 0, Z, &MWWWW0);
    success *= TSIL_Sanalytic (0, 0, h, Z, Q2, &S00h);
    success *= TSIL_Sanalytic (0, 0, W, Z, Q2, &S00W);
    success *= TSIL_Sanalytic (0, T, W, Z, Q2, &S0tW);
    success *= TSIL_Sanalytic (Z, Z, Z, Z, Q2, &SZZZ);
    success *= TSIL_Tanalytic (h, 0, 0, Z, Q2, &Th00);
    success *= TSIL_Tanalytic (T, 0, W, Z, Q2, &Tt0W);
    success *= TSIL_Tanalytic (W, 0, 0, Z, Q2, &TW00);
    success *= TSIL_Tanalytic (W, 0, T, Z, Q2, &TW0t);
    success *= TSIL_Tanalytic (Z, 0, 0, Z, Q2, &TZ00);
    success *= TSIL_Tbaranalytic (0, W, W, Z, Q2, &Tbar0WW);
    success *= TSIL_Uanalytic (h, Z, 0, 0, Z, Q2, &UhZ00);
    success *= TSIL_Uanalytic (W, W, 0, 0, Z, Q2, &UWW00);

    if (1 != success) SMDR_Error ("SMDR_MZ_DoTSILZ",
      "Analytic TSIL computation failed. This can't happen! Please contact the authors with details.", 1);

    TSIL_SetParameters (&bar, 0, T, 0, T, W, Q2);
    TSIL_Evaluate (&bar, Z);
    M0t0tW = TSIL_GetFunction (&bar, "M");
    Utt0W  = TSIL_GetFunction (&bar, "Uyuzv");

    TSIL_SetParameters (&bar, 0, W, 0, W, T, Q2);
    TSIL_Evaluate (&bar, Z);
    M0W0Wt = TSIL_GetFunction (&bar, "M");
    UWW0t  = TSIL_GetFunction (&bar, "Uyuzv");

    TSIL_SetParameters (&bar, h, h, Z, Z, h, Q2);
    TSIL_Evaluate (&bar, Z);
    MhhZZh = TSIL_GetFunction (&bar, "M");
    UhZhZ  = TSIL_GetFunction (&bar, "Uyuzv");
    UZhhh  = TSIL_GetFunction (&bar, "Uuyxv");
    ShhZ  = TSIL_GetFunction (&bar, "Svyz");
    ThhZ  = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, h, T, Z, T, T, Q2);
    TSIL_Evaluate (&bar, Z);
    MhtZtt = TSIL_GetFunction (&bar, "M");
    UZhtt = TSIL_GetFunction (&bar, "Uzxyv");
    UhZtt = TSIL_GetFunction (&bar, "Uxzuv");
    Uttht = TSIL_GetFunction (&bar, "Uuyxv");
    UtttZ = TSIL_GetFunction (&bar, "Uyuzv");
    Shtt  = TSIL_GetFunction (&bar, "Suxv");
    SttZ  = TSIL_GetFunction (&bar, "Svyz");
    Ttht  = TSIL_GetFunction (&bar, "Tuxv");
    Thtt  = TSIL_GetFunction (&bar, "Txuv");
    TttZ  = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, h, W, Z, W, W, Q2);
    TSIL_Evaluate (&bar, Z);
    MhWZWW = TSIL_GetFunction (&bar, "M");
    UZhWW = TSIL_GetFunction (&bar, "Uzxyv");
    UhZWW = TSIL_GetFunction (&bar, "Uxzuv");
    UWWhW = TSIL_GetFunction (&bar, "Uuyxv");
    UWWWZ = TSIL_GetFunction (&bar, "Uyuzv");
    ShWW  = TSIL_GetFunction (&bar, "Suxv");
    SWWZ  = TSIL_GetFunction (&bar, "Svyz");
    TWhW  = TSIL_GetFunction (&bar, "Tuxv");
    ThWW  = TSIL_GetFunction (&bar, "Txuv");
    TWWZ  = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, h, Z, Z, h, Z, Q2);
    TSIL_Evaluate (&bar, Z);
    MhZZhZ = TSIL_GetFunction (&bar, "M");
    UZhZZ = TSIL_GetFunction (&bar, "Uzxyv");

    TSIL_SetParameters (&bar, T, T, T, T, h, Q2);
    TSIL_Evaluate (&bar, Z);
    Mtttth = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, T, T, T, T, Z, Q2);
    TSIL_Evaluate (&bar, Z);
    MttttZ = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, T, W, T, W, 0, Q2);
    TSIL_Evaluate (&bar, Z);
    MtWtW0 = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, W, W, W, W, h, Q2);
    TSIL_Evaluate (&bar, Z);
    MWWWWh = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, W, W, W, W, Z, Q2);
    TSIL_Evaluate (&bar, Z);
    MWWWWZ = TSIL_GetFunction (&bar, "M");
  }

  return 0;
}

