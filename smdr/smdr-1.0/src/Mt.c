/* Top-quark pole mass calculation. In this file:
    method = 0 (expand around tree-level) or
             1 (expand around pole).
    QCDLoopOrder = 0, 1, 2, 3, 4.
    otherLoopOrder = 0 (tree-level) or 
                     1 (1-loop) or
                     1.5 (2-loop mixed QCD/EW) or
                     2 (2-loop full).

    4-loop pure QCD from 1502.01030, updated and perfected in 1606.06754.
    2-loop non-pure-QCD from 1604.01134.
*/
#include "smdr_internal.h"

#define MAXITERS 20

SMDR_REAL PM_TOL;

/* Local functions; not needed elsewhere: */
SMDR_COMPLEX SMDR_Mt_T1loopQCD (void);
SMDR_COMPLEX SMDR_Mt_T2loopQCD (void);
SMDR_COMPLEX SMDR_Mt_T3loopQCD (void);
SMDR_COMPLEX SMDR_Mt_T4loopQCD (void);
SMDR_COMPLEX SMDR_Mt_T1loopnonQCD (void);
SMDR_COMPLEX SMDR_Mt_T2loopmixed (void);
SMDR_COMPLEX SMDR_Mt_T2loopnonQCD (void);
SMDR_COMPLEX SMDR_Mt_t1loopQCD (void);
SMDR_COMPLEX SMDR_Mt_t2loopQCD (void);
SMDR_COMPLEX SMDR_Mt_t3loopQCD (void);
SMDR_COMPLEX SMDR_Mt_t4loopQCD (void);
SMDR_COMPLEX SMDR_Mt_t1loopnonQCD (void);
SMDR_COMPLEX SMDR_Mt_t2loopmixed (void);
SMDR_COMPLEX SMDR_Mt_t2loopnonQCD (void);
int SMDR_Mt_DoTSILt (float);

/* #define global variables used in this file, for convenience with safety. */
#define Ah SMDR_Mt_Ah
#define At SMDR_Mt_At
#define AW SMDR_Mt_AW
#define AZ SMDR_Mt_AZ
#define Ab SMDR_Mt_Ab
#define AT SMDR_Mt_AT
#define Ah2 SMDR_Mt_Ah2
#define At2 SMDR_Mt_At2
#define AW2 SMDR_Mt_AW2
#define AZ2 SMDR_Mt_AZ2
#define BhT SMDR_Mt_BhT
#define BTZ SMDR_Mt_BTZ
#define Bht SMDR_Mt_Bht
#define BtZ SMDR_Mt_BtZ
#define BbW SMDR_Mt_BbW
#define B0W SMDR_Mt_B0W
#define ReB0W SMDR_Mt_ReB0W
#define B0W2 SMDR_Mt_B0W2
#define Bht2 SMDR_Mt_Bht2
#define BtZ2 SMDR_Mt_BtZ2
#define I0hW SMDR_Mt_I0hW
#define I0hZ SMDR_Mt_I0hZ
#define I0tW SMDR_Mt_I0tW
#define I0WZ SMDR_Mt_I0WZ
#define Ihhh SMDR_Mt_Ihhh
#define Ihtt SMDR_Mt_Ihtt
#define IhWW SMDR_Mt_IhWW
#define IhZZ SMDR_Mt_IhZZ
#define IttZ SMDR_Mt_IttZ
#define IWWZ SMDR_Mt_IWWZ
#define S000 SMDR_Mt_S000
#define S0hW SMDR_Mt_S0hW
#define Th0t SMDR_Mt_Th0t
#define Th0W SMDR_Mt_Th0W
#define Thht SMDR_Mt_Thht
#define ThtZ SMDR_Mt_ThtZ
#define TthZ SMDR_Mt_TthZ
#define TW00 SMDR_Mt_TW00
#define TW0h SMDR_Mt_TW0h
#define TW0Z SMDR_Mt_TW0Z
#define TWtW SMDR_Mt_TWtW
#define TZ0t SMDR_Mt_TZ0t
#define TZ0W SMDR_Mt_TZ0W
#define TZht SMDR_Mt_TZht
#define TZtZ SMDR_Mt_TZtZ
#define Tbar00W SMDR_Mt_Tbar00W
#define Tbar0ht SMDR_Mt_Tbar0ht
#define Tbar0tZ SMDR_Mt_Tbar0tZ
#define U0W00 SMDR_Mt_U0W00
#define U0W0t SMDR_Mt_U0W0t
#define U0WhW SMDR_Mt_U0WhW
#define U0WWZ SMDR_Mt_U0WWZ
#define Uht0W SMDR_Mt_Uht0W
#define Uhtht SMDR_Mt_Uhtht
#define UhttZ SMDR_Mt_UhttZ
#define Ut0WW SMDR_Mt_Ut0WW
#define Uthhh SMDR_Mt_Uthhh
#define Uthtt SMDR_Mt_Uthtt
#define UthWW SMDR_Mt_UthWW
#define UthZZ SMDR_Mt_UthZZ
#define UtZ00 SMDR_Mt_UtZ00
#define UtZhZ SMDR_Mt_UtZhZ
#define UtZtt SMDR_Mt_UtZtt
#define UtZWW SMDR_Mt_UtZWW
#define UW00h SMDR_Mt_UW00h
#define UW00Z SMDR_Mt_UW00Z
#define UZt0W SMDR_Mt_UZt0W
#define UZtht SMDR_Mt_UZtht
#define UZttZ SMDR_Mt_UZttZ
#define M00tW0 SMDR_Mt_M00tW0
#define M00WW0 SMDR_Mt_M00WW0
#define M00WWZ SMDR_Mt_M00WWZ
#define M0tt0t SMDR_Mt_M0tt0t
#define Ut000 SMDR_Mt_Ut000
#define Tt00 SMDR_Mt_Tt00
#define Tbar00t SMDR_Mt_Tbar00t
#define M0ttht SMDR_Mt_M0ttht
#define M0ttZt SMDR_Mt_M0ttZt
#define M0tW0W SMDR_Mt_M0tW0W
#define M0tWhW SMDR_Mt_M0tWhW
#define M0tWZW SMDR_Mt_M0tWZW
#define M0ZWt0 SMDR_Mt_M0ZWt0
#define Mhhtth SMDR_Mt_Mhhtth
#define Mhttht SMDR_Mt_Mhttht
#define MhttZt SMDR_Mt_MhttZt
#define MhZttZ SMDR_Mt_MhZttZ
#define MttZZh SMDR_Mt_MttZZh
#define MtZZtt SMDR_Mt_MtZZtt
#define Tpole SMDR_Mt_Tpole
#define lnTpole SMDR_Mt_lnTpole
#define ln2Tpole SMDR_Mt_ln2Tpole
#define ln3Tpole SMDR_Mt_ln3Tpole
#define ln4Tpole SMDR_Mt_ln4Tpole
#define lnT SMDR_Mt_lnT
#define ln2T SMDR_Mt_ln2T
#define ln3T SMDR_Mt_ln3T
#define ln4T SMDR_Mt_ln4T
#define fiveTm2W2 SMDR_Mt_fiveTm2W2
#define fiveTm2W4 SMDR_Mt_fiveTm2W4
#define twoTpZ2 SMDR_Mt_twoTpZ2
#define TmW2 SMDR_Mt_TmW2

SMDR_COMPLEX Ah, At, AW, AZ, Ab, AT, Ah2, At2, AW2, AZ2,
  BhT, BTZ, Bht, BtZ, BbW, B0W, ReB0W, B0W2, Bht2, BtZ2,
  I0hW, I0hZ, I0tW, I0WZ, Ihhh, Ihtt, IhWW, IhZZ, IttZ, 
  IWWZ, S000, S0hW, Th0t, Th0W, Thht, ThtZ, TthZ, TW00, 
  TW0h, TW0Z, TWtW, TZ0t, TZ0W, TZht, TZtZ, Tbar00W, Tbar0ht, Tbar0tZ, 
  U0W00, U0W0t, U0WhW, U0WWZ, Uht0W, Uhtht, UhttZ, Ut0WW, Uthhh, Uthtt, 
  UthWW, UthZZ, UtZ00, UtZhZ, UtZtt, UtZWW, UW00h, UW00Z, UZt0W, UZtht, UZttZ, 
  M00tW0, M00WW0, M00WWZ, M0tt0t, Ut000, Tt00, Tbar00t,
  M0ttht, M0ttZt, M0tW0W, M0tWhW, M0tWZW, M0ZWt0, Mhhtth, Mhttht, 
  MhttZt, MhZttZ, MttZZh, MtZZtt;

SMDR_REAL Tpole, lnTpole, ln2Tpole, ln3Tpole, ln4Tpole;
SMDR_REAL lnT, ln2T, ln3T, ln4T;
SMDR_REAL fiveTm2W2, fiveTm2W4, twoTpZ2, TmW2;

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/*  Computes top-quark Yukawa coupling, given the pole mass. Arguments are:
    MTpoletarget,
    method = 0 (expand around tree-level) or
             1 (expand around pole),
    QCDLoopOrder = 0, 1, 2, 3, 4,
    otherLoopOrder = 0 (tree-level) or 
                     1 (1-loop) or
                     1.5 (2-loop mixed QCD/EW) or
                     2 (2-loop full).
*/
SMDR_REAL SMDR_Eval_yt (SMDR_REAL Q_eval,
                        SMDR_REAL MTpoletarget, 
                        int method, 
                        int QCDLoopOrder, 
                        float otherLoopOrder)
{
  SMDR_REAL yttemp;
  SMDR_REAL ytsave = yt;
  SMDR_COMPLEX CM2T;
  SMDR_REAL targetratio;
  int i;
  char funcname[] = "SMDR_Eval_yt";

  PM_TOL = 1.0e-8;
  if ((SMDR_MTPOLE_TOLERANCE < PM_TOL) && (SMDR_MTPOLE_TOLERANCE > 1.0e-15))
     PM_TOL = SMDR_MTPOLE_TOLERANCE;

  /* If Q_eval is negative, then we just use the current Q.
     Otherwise, we run all the parameters from Q_in to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

  if ( (method != 0) && (method != 1) )
    SMDR_Error (funcname, 
    "Invalid method specified, should be 0 or 1", 3);

  if ( (QCDLoopOrder < 0) || (QCDLoopOrder > 4)) SMDR_Error (funcname, 
    "Invalid QCDLoopOrder specified, should be 0, 1, 2, 3, or 4", 3);

  if ( (TSIL_FABS(otherLoopOrder) > 0.0001) &&
       (TSIL_FABS(otherLoopOrder-1) > 0.0001) &&
       (TSIL_FABS(otherLoopOrder-1.5) > 0.0001) &&
       (TSIL_FABS(otherLoopOrder-2) > 0.0001) )
    SMDR_Error (funcname, 
    "Invalid otherLoopOrder specified, should be 0, 1, 1.5, or 2", 3);

  /* Check input parameters for sanity: */
  SMDR_Check_Mtpole_Range (MTpoletarget, 170, 178);
  SMDR_Check_VEV_Range (199, 301);
  SMDR_Check_g3_Range (0.5, 1.5);
  SMDR_Check_g_Range (0.4, 0.8);
  SMDR_Check_gp_Range (0.25, 0.5);
  SMDR_Check_lambda_Range (0.05, 0.2);
  SMDR_Check_yb_Range (0.05);
  SMDR_Check_ytau_Range (0.02);
  SMDR_Check_Q_Range (49, 501);

  /* Initial guess. */
  yt = 0.93;
  SMDR_Update ();
  Tpole = MTpoletarget * MTpoletarget;

  if (1 == method) {
    T = Tpole;
    T2 = T*T; T3 = T2*T; T4 = T3*T;
    fiveTm2W2 = (5*T - 2*W) * (5*T - 2*W); 
    fiveTm2W4 = fiveTm2W2 * fiveTm2W2;
    twoTpZ2 = (2*T + Z) * (2*T + Z);
    TmW2 = (T - W) * (T - W);
    SMDR_Mt_DoTSILt (otherLoopOrder);
    for (i = 0; i < MAXITERS; i++) {
      CM2T = (yt2*v2)/2.L;
      if (QCDLoopOrder > 0.5) CM2T += ONELOOPFACTOR*SMDR_Mt_T1loopQCD();
      if (QCDLoopOrder > 1.5) CM2T += TWOLOOPFACTOR*SMDR_Mt_T2loopQCD();
      if (QCDLoopOrder > 2.5) CM2T += THREELOOPFACTOR*SMDR_Mt_T3loopQCD();
      if (QCDLoopOrder > 3.5) CM2T += FOURLOOPFACTOR*SMDR_Mt_T4loopQCD();
      if (otherLoopOrder > 0.1) CM2T += ONELOOPFACTOR*SMDR_Mt_T1loopnonQCD();
      if (otherLoopOrder > 1.1) CM2T += TWOLOOPFACTOR*SMDR_Mt_T2loopmixed();
      if (otherLoopOrder > 1.9) CM2T += TWOLOOPFACTOR*SMDR_Mt_T2loopnonQCD();
      targetratio = MTpoletarget/SMDR_SQRT(SMDR_CREAL(CM2T));
      yt = yt * targetratio;
      yt2 = yt*yt;
      /*
      printf("eval_yt targetratio = %Lf,   yt = %Lf\n",targetratio,yt); 
      */
      if ( TSIL_FABS(targetratio - 1.L) < PM_TOL) break;
    }
  }

  if (0 == method) {
    for (i = 0; i < MAXITERS; i++) {
/* Changed first argument below from Q to -1 on 2/18/2019. Check
   that's the right thing?? */
      SMDR_Eval_Mt_pole (-1, method, QCDLoopOrder, otherLoopOrder,
                         &SMDR_Mt_pole, &SMDR_Gammat_pole);
      targetratio = MTpoletarget/SMDR_Mt_pole;
      yt = yt * targetratio; 
      /*
      printf("eval_yt targetratio = %Lf,   yt = %Lf\n",targetratio,yt);
      */ 
      if ( TSIL_FABS(targetratio - 1.L) < PM_TOL) break;
    }
  }

  yttemp = yt;
  yt = ytsave; 
  SMDR_Update();
  return(yttemp);
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/*  Computes top pole mass. Arguments are:
    method = 0 (expand around tree-level) or
             1 (expand around pole).
    QCDLoopOrder = 0, 1, 2, 3, 4.
    otherLoopOrder = 0 (tree-level) or 
                     1 (1-loop) or
                     1.5 (2-loop mixed QCD/EW) or
                     2 (2-loop full).
*/
void SMDR_Eval_Mt_pole (SMDR_REAL Q_eval,
                        int method, 
                        int QCDLoopOrder, 
                        float otherLoopOrder,
                        SMDR_REAL *Mtpoleresult,
                        SMDR_REAL *Gammatpoleresult)
{
  int i; 
  SMDR_COMPLEX CM2T, CM2Ttemp;
  SMDR_REAL Ttemp;
  char funcname[] = "SMDR_Eval_Mt_pole";

  PM_TOL = 1.0e-8;
  if ((SMDR_MTPOLE_TOLERANCE < PM_TOL) && (SMDR_MTPOLE_TOLERANCE > 1.0e-15))
     PM_TOL = SMDR_MTPOLE_TOLERANCE;

  /* If Q_eval is negative, then we just use the current Q.
     Otherwise, we run all the parameters from Q_in to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

  if ( (method != 0) && (method != 1) )
    SMDR_Error (funcname, 
    "Invalid method specified, should be 0 or 1", 3);

  if ( (QCDLoopOrder < 0) || (QCDLoopOrder > 4)) SMDR_Error (funcname, 
    "Invalid QCDLoopOrder specified, should be 0, 1, 2, 3, or 4", 3);

  if ( (TSIL_FABS(otherLoopOrder) > 0.0001) &&
       (TSIL_FABS(otherLoopOrder-1) > 0.0001) &&
       (TSIL_FABS(otherLoopOrder-1.5) > 0.0001) &&
       (TSIL_FABS(otherLoopOrder-2) > 0.0001) )
    SMDR_Error (funcname, 
    "Invalid otherLoopOrder specified, should be 0, 1, 1.5, or 2", 3);

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
  fiveTm2W2 = (5*T - 2*W) * (5*T - 2*W); 
  fiveTm2W4 = fiveTm2W2 * fiveTm2W2;
  twoTpZ2 = (2*T + Z) * (2*T + Z);
  TmW2 = (T - W) * (T - W);

  Tpole = Ttemp = T;
  CM2T = T;

  if (0 == method) {
    SMDR_Mt_DoTSILt (otherLoopOrder);
    if (QCDLoopOrder > 0.5) CM2T += ONELOOPFACTOR*SMDR_Mt_t1loopQCD();
    if (QCDLoopOrder > 1.5) CM2T += TWOLOOPFACTOR*SMDR_Mt_t2loopQCD();
    if (QCDLoopOrder > 2.5) CM2T += THREELOOPFACTOR*SMDR_Mt_t3loopQCD();
    if (QCDLoopOrder > 3.5) CM2T += FOURLOOPFACTOR*SMDR_Mt_t4loopQCD();
    if (otherLoopOrder > 0.1) CM2T += ONELOOPFACTOR*SMDR_Mt_t1loopnonQCD();
    if (otherLoopOrder > 1.1) CM2T += TWOLOOPFACTOR*SMDR_Mt_t2loopmixed();
    if (otherLoopOrder > 1.9) CM2T += TWOLOOPFACTOR*SMDR_Mt_t2loopnonQCD();
  } 

  if (1 == method) {
    for (i = 0; i < MAXITERS; i++) {
      T = Tpole;
      T2 = T*T; T3 = T2*T; T4 = T3*T;
      fiveTm2W2 = (5*T - 2*W) * (5*T - 2*W); 
      fiveTm2W4 = fiveTm2W2 * fiveTm2W2;
      twoTpZ2 = (2*T + Z) * (2*T + Z);
      TmW2 = (T - W) * (T - W);
      CM2Ttemp = CM2T;
      CM2T = Ttemp;
      SMDR_Mt_DoTSILt (otherLoopOrder);
      if (QCDLoopOrder > 0.5) CM2T += ONELOOPFACTOR*SMDR_Mt_T1loopQCD();
      if (QCDLoopOrder > 1.5) CM2T += TWOLOOPFACTOR*SMDR_Mt_T2loopQCD();
      if (QCDLoopOrder > 2.5) CM2T += THREELOOPFACTOR*SMDR_Mt_T3loopQCD();
      if (QCDLoopOrder > 3.5) CM2T += FOURLOOPFACTOR*SMDR_Mt_T4loopQCD();
      if (otherLoopOrder > 0.1) CM2T += ONELOOPFACTOR*SMDR_Mt_T1loopnonQCD();
      if (otherLoopOrder > 1.1) CM2T += TWOLOOPFACTOR*SMDR_Mt_T2loopmixed();
      if (otherLoopOrder > 1.9) CM2T += TWOLOOPFACTOR*SMDR_Mt_T2loopnonQCD();
      Tpole = SMDR_CREAL(CM2T);
      if ( TSIL_CABS((CM2T/CM2Ttemp) - 1.L) < PM_TOL) break;
    }
  }

  *Mtpoleresult = SMDR_SQRT(SMDR_CREAL(CM2T));
  *Gammatpoleresult = -SMDR_CIMAG(CM2T)/(*Mtpoleresult);

  return;
}

/* ------------------------------------------------------------ */
/* 1604.01134 equations (2.1) and (2.2)                         */

SMDR_COMPLEX SMDR_Mt_t1loopQCD ()
{
  SMDR_COMPLEX result;
  result = g32 * T * ((32./3.)  - 8 * lnT);
  return (result);
}

/* ------------------------------------------------------------------ */
/* 1604.01134 equations (2.1) and (2.3).                              */

SMDR_COMPLEX SMDR_Mt_t2loopQCD ()
{
  SMDR_REAL c20 = 292.01441887879963L;
  SMDR_REAL c21 = -204.L;
  SMDR_REAL c22 = 60.L;
  SMDR_REAL c2b = 52.63789013914324L;
  SMDR_COMPLEX result;

  result = g34 * T * (c20 + c21 * lnT + c22 * ln2T + c2b * yb/yt);
  return (result);
}

/* ------------------------------------------------------------ */
/* 1604.01134 equations (2.1), (2.9) and (2.11)-(2.14)          */

SMDR_COMPLEX SMDR_Mt_t3loopQCD ()
{
  SMDR_REAL c30 = 10831.186915635191L;
  SMDR_REAL c31 = -5977.0007973114525L;
  SMDR_REAL c32 = 2300.L;
  SMDR_REAL c33 = -440.L;
  return (g36 * T * (c30 + c31 * lnT + c32 * ln2T + c33 * ln3T));
}

/* ------------------------------------------------------------ */
/* 1604.01134 equations (2.1), (2.10) and (2.15)-(2.19)         */

SMDR_COMPLEX SMDR_Mt_t4loopQCD ()
{
  /* Old and busted from 1502.01030:
  SMDR_REAL c40 = 491247.57540717465L;
  */
  /* New and cool from 1606.06754: */
  SMDR_REAL c40 = 497274.; /* Uncertainty 921.6 */
  SMDR_REAL c41 = -289975.47247110726L;
  SMDR_REAL c42 = 80453.91255559417L;
  SMDR_REAL c43 = -21913.333333333332L;
  SMDR_REAL c44 = 3190.L;
  return (g38 * T * (c40 + c41 * lnT + c42 * ln2T + c43 * ln3T + c44 * ln4T));
}

/* ------------------------------------------------------------ */
/* 1604.01134 equations (2.1) and (2.20)                        */

SMDR_COMPLEX SMDR_Mt_t1loopnonQCD ()
{
  SMDR_COMPLEX BFVbW, BFVTZ, result;

  BbW = TSIL_B (b, W, T, Q2);

  BFVbW = (W - b - T) * BbW + AW - Ab + T + 
          ((W*(T + b) - (b - T)*(b - T)) * BbW + (b - T) * AW)/(2.*W);

  BFVTZ = (Z - T) * BtZ + AZ - At + T; 

  result = (8.L/9.L) * (g2 * gp2/(g2 + gp2)) * (T - 3.L * At) +
    (yt2/2.L) * ((h-4.L*T)*Bht + Ah - 2*At - Ab) - (yb2/2.) * Ab +
    (g2/2.L) * BFVbW + 
    ((g2 + gp2)/4.L -(2.L/3.L)*gp2 +(8.L/9.L) * gp2*gp2/(g2 + gp2)) * BFVTZ +
    ((4.L/3.L) * gp2)*((2.L/3.L)*gp2/(g2 + gp2) - 0.5L) * T * (3.L*BtZ - 2.L);

  return result;
}

/* --------------------------------------------------------------- */
/* 1604.01134 eq. (2.1) and ancillary file delta2mixed_secII.txt   */

SMDR_COMPLEX SMDR_Mt_t2loopmixed ()
{
  SMDR_COMPLEX t2g32;
  #include "includes/t2g32.c"
  return (t2g32 * g32/v2); 
}

/* --------------------------------------------------------------- */
/* 1604.01134 eq. (2.1) and ancillary file delta2nonQCD_secII.txt  */

SMDR_COMPLEX SMDR_Mt_t2loopnonQCD (void) 
{
  SMDR_COMPLEX t2EW;
  #include "includes/t2EW.c"
  t2EW = t2EW/(v2*v2);
  return (t2EW);
}

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */

/* ------------------------------------------------------------ */
/* 1604.01134 equations (3.2) and (3.3)                         */

SMDR_COMPLEX SMDR_Mt_T1loopQCD ()
{
  return (g32 * Tpole * ((32./3.) - 8 * lnTpole));
}

/* ------------------------------------------------------------ */
/* 1604.01134 equations (3.2) and (3.4)                         */

SMDR_COMPLEX SMDR_Mt_T2loopQCD ()
{
  SMDR_REAL a20 = 263.56997443435523L;
  SMDR_REAL a21 = -97.33333333333333L;
  SMDR_REAL a22 = -4.L;
  SMDR_REAL a2b = 52.63789013914324L;

  return (g34 * Tpole * (a20 + a21 * lnTpole + a22 * ln2Tpole + a2b * yb/yt));
}

/* ------------------------------------------------------------ */
/* 1604.01134 equations (3.2), (3.5), and (3.7)-(3.10)          */

SMDR_COMPLEX SMDR_Mt_T3loopQCD ()
{
  SMDR_REAL a30 = 8734.402071325276L;
  SMDR_REAL a31 = -1326.103428583991L;
  SMDR_REAL a32 = -36.L;
  SMDR_REAL a33 = 8.L;
  return (g36 * Tpole * (a30 + a31 * lnTpole + a32 * ln2Tpole + 
          a33 * ln3Tpole));
}

/* ------------------------------------------------------------ */
/* 1604.01134 equations (3.2), (3.6), and (3.11)-(3.15)         */

SMDR_COMPLEX SMDR_Mt_T4loopQCD ()
{
  /* Old and busted from 1502.01030:
  SMDR_REAL a40 = 364091.74916320975L;
  */
  /* New and cool from 1606.06754: */
  SMDR_REAL a40 = 370118.; /* Uncertainty 921.6 */
  SMDR_REAL a41 = -89526.05631921464L;
  SMDR_REAL a42 = 3286.9379578845037L;
  SMDR_REAL a43 = 81.33333333333333L;
  SMDR_REAL a44 = -26.L;
  return (g38 * Tpole * (a40 + a41 * lnTpole + a42 * ln2Tpole +
          a43 * ln3Tpole + a44 * ln4Tpole));
}

/* ------------------------------------------------------------ */
/* 1604.01134 equations (3.2) and (3.16)                        */

SMDR_COMPLEX SMDR_Mt_T1loopnonQCD ()
{
  SMDR_COMPLEX BFVbW, BFVTZ, result;

  BFVbW = (W - b - Tpole) * BbW + AW - Ab + Tpole + 
          ((W*(Tpole + b) - (b - Tpole)*(b - Tpole)) * BbW + 
          (b - Tpole) * AW)/(2.*W);

  BFVTZ = (Z - Tpole)*BTZ + AZ - AT + Tpole; 

  result = (8.L/9.L) * (g2 * gp2/(g2 + gp2)) * (Tpole - 3.L * AT) +
    (yt2/2.L) * ((h-4.L*Tpole)*BhT + Ah - 2*AT - Ab) - (yb2/2.) * Ab +
    (g2/2.L) * BFVbW +
    ((g2 + gp2)/4.L - (2.L/3.L)*gp2 + 
    (8.L/9.L) * gp2 * gp2/(g2 + gp2)) * BFVTZ +
    ((4.L/3.L) * gp2)*((2.L/3.L)*gp2/(g2 + gp2) - 0.5L) * 
    Tpole * (3.L*BTZ - 2.L);

  return result;
}

/* ---------------------------------------------------------------- */
/* 1604.01134 eq. (2.1) and ancillary file Delta2mixed_secIII.txt   */

SMDR_COMPLEX SMDR_Mt_T2loopmixed()
{
  SMDR_COMPLEX T2g32;
  #include "includes/tt2g32.c"
  return (T2g32 * g32/v2);
}

/* ------------------------------------------------------------------ */
/* 1604.01134 eq. (2.1) and ancillary file Delta2nonQCD_secIII.txt    */

SMDR_COMPLEX SMDR_Mt_T2loopnonQCD (void) 
{
  SMDR_COMPLEX T2EW;
  #include "includes/tt2EW.c"
  return (T2EW/v4);
}

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* Performs the needed basis integral evaluations.                    */

int SMDR_Mt_DoTSILt (float otherLoopOrder)
{
  TSIL_DATA bar;
  int success = 1;

  lnTpole = TSIL_CLOG (Tpole/Q2);
  ln2Tpole = lnTpole * lnTpole;
  ln3Tpole = ln2Tpole * lnTpole;
  ln4Tpole = ln3Tpole * lnTpole;

  lnT = TSIL_CLOG (T/Q2);
  ln2T = lnT * lnT;
  ln3T = ln2T * lnT;
  ln4T = ln3T * lnT;

  At = TSIL_A (T, Q2);
  Ab = TSIL_A (b, Q2);
  AW = TSIL_A (W, Q2);
  AZ = TSIL_A (Z, Q2);
  Ah = TSIL_A (h, Q2);
  BtZ = TSIL_B (T, Z, T, Q2);
  Bht = TSIL_B (h, T, T, Q2);
  B0W = TSIL_B (0, W, T, Q2);
  ReB0W = SMDR_CREAL(B0W);

  At2 = At * At;
  AW2 = AW * AW;
  AZ2 = AZ * AZ;
  Ah2 = Ah * Ah;
  B0W2 = B0W * B0W;
  BtZ2 = BtZ * BtZ;
  Bht2 = Bht * Bht;

  AT = TSIL_A (Tpole, Q2);
  BTZ = TSIL_B (Tpole, Z, Tpole, Q2);
  BhT = TSIL_B (h, Tpole, Tpole, Q2);
  BbW = TSIL_B (b, W, Tpole, Q2);

  /* These are needed for the 2-loop mixed QCD contributions: */
  if (otherLoopOrder > 1.1) {

    I0tW = TSIL_I2 (0, T, W, Q2);
    Ihtt = TSIL_I2 (h, T, T, Q2);
    IttZ = TSIL_I2 (T, T, Z, Q2);

    success *= TSIL_Manalytic (0, T, T, 0, T, T, &M0tt0t);
    success *= TSIL_Uanalytic (T, 0, 0, 0, T, Q2, &Ut000);
    success *= TSIL_Tanalytic (T, 0, 0, T, Q2, &Tt00);
    success *= TSIL_Tbaranalytic (0, 0, T, T, Q2, &Tbar00t);

    success *= TSIL_Manalytic (0, 0, T, W, 0, T, &M00tW0);
    success *= TSIL_Uanalytic (0, W, 0, T, T, Q2, &U0W0t);
    success *= TSIL_Tanalytic (W, 0, 0, T, Q2, &TW00);
    success *= TSIL_Tbaranalytic (0, 0, W, T, Q2, &Tbar00W);

    if (1 != success) printf ("hellno!!!\n");

    TSIL_SetParameters (&bar, 0, T, T, h, T, Q2);
    TSIL_Evaluate (&bar, T);
    M0ttht = TSIL_GetFunction (&bar, "M");
    Uthtt  = TSIL_GetFunction (&bar, "Uyuzv");
    Th0t  = TSIL_GetFunction (&bar, "Tuxv");
    Tbar0ht  = TSIL_GetFunction (&bar, "TBARxuv");

    TSIL_SetParameters (&bar, 0, T, T, Z, T, Q2);
    TSIL_Evaluate (&bar, T);
    M0ttZt = TSIL_GetFunction (&bar, "M");
    UtZtt  = TSIL_GetFunction (&bar, "Uyuzv");
    TZ0t  = TSIL_GetFunction (&bar, "Tuxv");
    Tbar0tZ  = TSIL_GetFunction (&bar, "TBARxuv");
  }

  /* These are needed for 2-loop nonQCD part: */
  if (otherLoopOrder > 1.9) {

    I0hW = TSIL_I2 (0, h, W, Q2);
    I0hZ = TSIL_I2 (0, h, Z, Q2);
    I0WZ = TSIL_I2 (0, W, Z, Q2);
    Ihhh = TSIL_I2 (h, h, h, Q2);
    IhWW = TSIL_I2 (h, W, W, Q2);
    IhZZ = TSIL_I2 (h, Z, Z, Q2);
    IWWZ = TSIL_I2 (W, W, Z, Q2);

    success *= TSIL_Manalytic (0, 0, W, W, 0, T, &M00WW0);
    success *= TSIL_Manalytic (0, T, W, 0, W, T, &M0tW0W);
    success *= TSIL_Tanalytic (W, T, W, T, Q2, &TWtW);
    success *= TSIL_Uanalytic (T, 0, W, W, T, Q2, &Ut0WW);
    success *= TSIL_Uanalytic (0, W, 0, 0, T, Q2, &U0W00);
    success *= TSIL_Uanalytic (W, 0, 0, h, T, Q2, &UW00h);
    success *= TSIL_Sanalytic (0, 0, 0, T, Q2, &S000);
    success *= TSIL_Tanalytic (h, h, T, T, Q2, &Thht);

    if (1 != success) printf ("fuckno!!!\n");

    TSIL_SetParameters (&bar, 0, 0, W, W, Z, Q2);
    TSIL_Evaluate (&bar, T);
    M00WWZ = TSIL_GetFunction (&bar, "M");
    TZ0W  = TSIL_GetFunction (&bar, "Tvxu");
    TW0Z  = TSIL_GetFunction (&bar, "Tuxv");
    U0WWZ  = TSIL_GetFunction (&bar, "Uxzuv");
    UW00Z  = TSIL_GetFunction (&bar, "Uzxyv");

    TSIL_SetParameters (&bar, 0, T, W, h, W, Q2);
    TSIL_Evaluate (&bar, T);
    M0tWhW = TSIL_GetFunction (&bar, "M");
    Uht0W  = TSIL_GetFunction (&bar, "Uuyxv");
    UthWW  = TSIL_GetFunction (&bar, "Uyuzv");
    U0WhW  = TSIL_GetFunction (&bar, "Uxzuv");
    S0hW  = TSIL_GetFunction (&bar, "Sxuv");
    Th0W  = TSIL_GetFunction (&bar, "Tuxv");
    TW0h  = TSIL_GetFunction (&bar, "Tvxu");

    TSIL_SetParameters (&bar, 0, T, W, Z, W, Q2);
    TSIL_Evaluate (&bar, T);
    M0tWZW = TSIL_GetFunction (&bar, "M");
    UZt0W  = TSIL_GetFunction (&bar, "Uuyxv");
    UtZWW  = TSIL_GetFunction (&bar, "Uyuzv");

    TSIL_SetParameters (&bar, 0, Z, W, T, 0, Q2);
    TSIL_Evaluate (&bar, T);
    M0ZWt0 = TSIL_GetFunction (&bar, "M");
    UtZ00  = TSIL_GetFunction (&bar, "Uuyxv");

    TSIL_SetParameters (&bar, h, h, T, T, h, Q2);
    TSIL_Evaluate (&bar, T);
    Mhhtth = TSIL_GetFunction (&bar, "M");
    Uhtht  = TSIL_GetFunction (&bar, "Uxzuv");
    Uthhh  = TSIL_GetFunction (&bar, "Uzxyv");

    TSIL_SetParameters (&bar, h, T, T, h, T, Q2);
    TSIL_Evaluate (&bar, T);
    Mhttht = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, h, T, T, Z, T, Q2);
    TSIL_Evaluate (&bar, T);
    MhttZt = TSIL_GetFunction (&bar, "M");
    UhttZ  = TSIL_GetFunction (&bar, "Uxzuv");
    UZtht  = TSIL_GetFunction (&bar, "Uuyxv");
    ThtZ = TSIL_GetFunction (&bar, "Txuv");
    TthZ = TSIL_GetFunction (&bar, "Tvxu");
    TZht = TSIL_GetFunction (&bar, "Tuxv");

    TSIL_SetParameters (&bar, h, Z, T, T, Z, Q2);
    TSIL_Evaluate (&bar, T);
    MhZttZ = TSIL_GetFunction (&bar, "M");
    UthZZ  = TSIL_GetFunction (&bar, "Uzxyv");
    UtZhZ  = TSIL_GetFunction (&bar, "Uuyxv");
    UZttZ  = TSIL_GetFunction (&bar, "Uyuzv");
    TZtZ = TSIL_GetFunction (&bar, "Tvyz");

    TSIL_SetParameters (&bar, T, T, Z, Z, h, Q2);
    TSIL_Evaluate (&bar, T);
    MttZZh = TSIL_GetFunction (&bar, "M");

    TSIL_SetParameters (&bar, T, Z, Z, T, T, Q2);
    TSIL_Evaluate (&bar, T);
    MtZZtt = TSIL_GetFunction (&bar, "M");
  }

  return 0;
}
