/* 
   Minimization of the SM effective potential in Landau gauge at up to 
   full 3-loop order (with 4-loop results at leading order in QCD). 
   For the functions defined in this file:

   loopOrder = 0 = tree-level
   loopOrder = 1 = 1-loop
   loopOrder = 2 = 2-loop 
   loopOrder = 2.5 = 2-loop plus 3-loop contributions in the large g3,yt limit
   loopOrder = 3 = full 3-loop
   loopOrder = 3.5 = full 3-loop plus 4-loop at leading order in QCD
*/

#include "smdr_internal.h"

/* Parameters governing iteration in SMDR_Eval_vev() */
#define MAXITERS_FINDVEV 15
#define TOLERANCE_FINDVEV 1.e-9

/* Local function; not needed elsewhere: */
int SMDR_effpot_Do3VILv (float loopOrder);

/* #define global variables used in this file, for convenience with safety. */
#define lnbarh SMDR_effpot_lnbarh
#define lnbarT SMDR_effpot_lnbarT
#define lnbarT2 SMDR_effpot_lnbarT2
#define lnbarT3 SMDR_effpot_lnbarT3
#define lnbarT4 SMDR_effpot_lnbarT4
#define AW SMDR_effpot_AW
#define AZ SMDR_effpot_AZ
#define Ah SMDR_effpot_Ah
#define AT SMDR_effpot_AT
#define Ab SMDR_effpot_Ab
#define Atau SMDR_effpot_Atau
#define I0hW SMDR_effpot_I0hW
#define I0hZ SMDR_effpot_I0hZ
#define I0TW SMDR_effpot_I0TW
#define I0WZ SMDR_effpot_I0WZ
#define Ihhh SMDR_effpot_Ihhh
#define IhTT SMDR_effpot_IhTT
#define IhWW SMDR_effpot_IhWW
#define IhZZ SMDR_effpot_IhZZ
#define ITTZ SMDR_effpot_ITTZ
#define IWWZ SMDR_effpot_IWWZ
#define H00000h SMDR_effpot_H00000h
#define H00000T SMDR_effpot_H00000T
#define H00000W SMDR_effpot_H00000W
#define H00000Z SMDR_effpot_H00000Z
#define H0000hW SMDR_effpot_H0000hW
#define H0000TT SMDR_effpot_H0000TT
#define H0000TW SMDR_effpot_H0000TW
#define H0000WW SMDR_effpot_H0000WW
#define H0000WZ SMDR_effpot_H0000WZ
#define H000hTT SMDR_effpot_H000hTT
#define H000hWW SMDR_effpot_H000hWW
#define H000hZZ SMDR_effpot_H000hZZ
#define H000TTZ SMDR_effpot_H000TTZ
#define H000WWZ SMDR_effpot_H000WWZ
#define H00h0WW SMDR_effpot_H00h0WW
#define H00hhWW SMDR_effpot_H00hhWW
#define H00hhZZ SMDR_effpot_H00hhZZ
#define H00hW00 SMDR_effpot_H00hW00
#define H00hW0Z SMDR_effpot_H00hW0Z
#define H00hWWZ SMDR_effpot_H00hWWZ
#define H00hZ00 SMDR_effpot_H00hZ00
#define H00hZ0W SMDR_effpot_H00hZ0W
#define H00hZWW SMDR_effpot_H00hZWW
#define H00T0TW SMDR_effpot_H00T0TW
#define H00T0WW SMDR_effpot_H00T0WW
#define H00ThTT SMDR_effpot_H00ThTT
#define H00TW0T SMDR_effpot_H00TW0T
#define H00TZ00 SMDR_effpot_H00TZ00
#define H00TZ0W SMDR_effpot_H00TZ0W
#define H00TZWW SMDR_effpot_H00TZWW
#define H00W0WZ SMDR_effpot_H00W0WZ
#define H00Whhh SMDR_effpot_H00Whhh
#define H00WhZZ SMDR_effpot_H00WhZZ
#define H00WThT SMDR_effpot_H00WThT
#define H00WTTZ SMDR_effpot_H00WTTZ
#define H00WW00 SMDR_effpot_H00WW00
#define H00WW0h SMDR_effpot_H00WW0h
#define H00WW0W SMDR_effpot_H00WW0W
#define H00WW0Z SMDR_effpot_H00WW0Z
#define H00WWhW SMDR_effpot_H00WWhW
#define H00WWWZ SMDR_effpot_H00WWWZ
#define H00WZ00 SMDR_effpot_H00WZ00
#define H00WZ0h SMDR_effpot_H00WZ0h
#define H00WZhZ SMDR_effpot_H00WZhZ
#define H00WZTT SMDR_effpot_H00WZTT
#define H00Zhhh SMDR_effpot_H00Zhhh
#define H00ZhWW SMDR_effpot_H00ZhWW
#define H00ZWhW SMDR_effpot_H00ZWhW
#define H00ZZ00 SMDR_effpot_H00ZZ00
#define H00ZZ0W SMDR_effpot_H00ZZ0W
#define H00ZZWW SMDR_effpot_H00ZZWW
#define H0hhWhW SMDR_effpot_H0hhWhW
#define H0hhZhZ SMDR_effpot_H0hhZhZ
#define H0hTZTT SMDR_effpot_H0hTZTT
#define H0hWWWh SMDR_effpot_H0hWWWh
#define H0hWWWZ SMDR_effpot_H0hWWWZ
#define H0hZWZW SMDR_effpot_H0hZWZW
#define H0hZZZh SMDR_effpot_H0hZZZh
#define H0TTT0T SMDR_effpot_H0TTT0T
#define H0TTThT SMDR_effpot_H0TTThT
#define H0TTTZT SMDR_effpot_H0TTTZT
#define H0TTW0W SMDR_effpot_H0TTW0W
#define H0TTWhW SMDR_effpot_H0TTWhW
#define H0TTWZW SMDR_effpot_H0TTWZW
#define H0WWW0W SMDR_effpot_H0WWW0W
#define H0WWWhW SMDR_effpot_H0WWWhW
#define H0WWWZW SMDR_effpot_H0WWWZW
#define H0WWZhZ SMDR_effpot_H0WWZhZ
#define H0WZZWW SMDR_effpot_H0WZZWW
#define Hhhhhhh SMDR_effpot_Hhhhhhh
#define HhhThTT SMDR_effpot_HhhThTT
#define HhhWhWW SMDR_effpot_HhhWhWW
#define HhhZhZZ SMDR_effpot_HhhZhZZ
#define HhTTThT SMDR_effpot_HhTTThT
#define HhTTTZT SMDR_effpot_HhTTTZT
#define HhTZTTZ SMDR_effpot_HhTZTTZ
#define HhWWWhW SMDR_effpot_HhWWWhW
#define HhWWWZW SMDR_effpot_HhWWWZW
#define HhWZWWZ SMDR_effpot_HhWZWWZ
#define HhZZZhZ SMDR_effpot_HhZZZhZ
#define HTTZZTT SMDR_effpot_HTTZZTT
#define HWWZZWW SMDR_effpot_HWWZZWW
#define G000hT SMDR_effpot_G000hT
#define G000TZ SMDR_effpot_G000TZ
#define G0TTWW SMDR_effpot_G0TTWW
#define Gh000W SMDR_effpot_Gh000W
#define Gh000Z SMDR_effpot_Gh000Z
#define Gh00hh SMDR_effpot_Gh00hh
#define Gh00TT SMDR_effpot_Gh00TT
#define Gh00WW SMDR_effpot_Gh00WW
#define Gh00ZZ SMDR_effpot_Gh00ZZ
#define Gh0W0Z SMDR_effpot_Gh0W0Z
#define Gh0Whh SMDR_effpot_Gh0Whh
#define Gh0WTT SMDR_effpot_Gh0WTT
#define Gh0WWW SMDR_effpot_Gh0WWW
#define Gh0WZZ SMDR_effpot_Gh0WZZ
#define Gh0Zhh SMDR_effpot_Gh0Zhh
#define Gh0ZTT SMDR_effpot_Gh0ZTT
#define Gh0ZWW SMDR_effpot_Gh0ZWW
#define Gh0ZZZ SMDR_effpot_Gh0ZZZ
#define Ghhhhh SMDR_effpot_Ghhhhh
#define GhhhTT SMDR_effpot_GhhhTT
#define GhhhWW SMDR_effpot_GhhhWW
#define GhhhZZ SMDR_effpot_GhhhZZ
#define GhTTWW SMDR_effpot_GhTTWW
#define GhTTZZ SMDR_effpot_GhTTZZ
#define GhWWZZ SMDR_effpot_GhWWZZ
#define GT000W SMDR_effpot_GT000W
#define GT00hT SMDR_effpot_GT00hT
#define GT00TZ SMDR_effpot_GT00TZ
#define GT0WhT SMDR_effpot_GT0WhT
#define GT0WTZ SMDR_effpot_GT0WTZ
#define GThTTZ SMDR_effpot_GThTTZ
#define GW000Z SMDR_effpot_GW000Z
#define GW00hW SMDR_effpot_GW00hW
#define GW00WZ SMDR_effpot_GW00WZ
#define GW0h0T SMDR_effpot_GW0h0T
#define GW0h0Z SMDR_effpot_GW0h0Z
#define GW0hhW SMDR_effpot_GW0hhW
#define GW0hWZ SMDR_effpot_GW0hWZ
#define GW0T0Z SMDR_effpot_GW0T0Z
#define GW0ThW SMDR_effpot_GW0ThW
#define GW0TWZ SMDR_effpot_GW0TWZ
#define GW0ZhW SMDR_effpot_GW0ZhW
#define GW0ZWZ SMDR_effpot_GW0ZWZ
#define GWhWWZ SMDR_effpot_GWhWWZ
#define GZ00hZ SMDR_effpot_GZ00hZ
#define GZ00TT SMDR_effpot_GZ00TT
#define GZ00WW SMDR_effpot_GZ00WW
#define GZ0h0W SMDR_effpot_GZ0h0W
#define GZ0hhZ SMDR_effpot_GZ0hhZ
#define GZ0hTT SMDR_effpot_GZ0hTT
#define GZ0hWW SMDR_effpot_GZ0hWW
#define GZ0WhZ SMDR_effpot_GZ0WhZ
#define GZ0WTT SMDR_effpot_GZ0WTT
#define GZ0WWW SMDR_effpot_GZ0WWW
#define GZhZTT SMDR_effpot_GZhZTT
#define GZhZWW SMDR_effpot_GZhZWW
#define GZTTWW SMDR_effpot_GZTTWW
#define Fh00T SMDR_effpot_Fh00T
#define Fh00W SMDR_effpot_Fh00W
#define Fh00Z SMDR_effpot_Fh00Z
#define Fh0hW SMDR_effpot_Fh0hW
#define Fh0hZ SMDR_effpot_Fh0hZ
#define Fh0TT SMDR_effpot_Fh0TT
#define Fh0TW SMDR_effpot_Fh0TW
#define Fh0WW SMDR_effpot_Fh0WW
#define Fh0WZ SMDR_effpot_Fh0WZ
#define Fh0ZZ SMDR_effpot_Fh0ZZ
#define FhhTT SMDR_effpot_FhhTT
#define FhhWW SMDR_effpot_FhhWW
#define FhhZZ SMDR_effpot_FhhZZ
#define FhTTZ SMDR_effpot_FhTTZ
#define FhWWZ SMDR_effpot_FhWWZ
#define FT00W SMDR_effpot_FT00W
#define FT00Z SMDR_effpot_FT00Z
#define FT0hW SMDR_effpot_FT0hW
#define FT0TW SMDR_effpot_FT0TW
#define FT0WZ SMDR_effpot_FT0WZ
#define FThTZ SMDR_effpot_FThTZ
#define FTTWW SMDR_effpot_FTTWW
#define FTTZZ SMDR_effpot_FTTZZ
#define FW00Z SMDR_effpot_FW00Z
#define FW0hh SMDR_effpot_FW0hh
#define FW0hT SMDR_effpot_FW0hT
#define FW0hZ SMDR_effpot_FW0hZ
#define FW0TT SMDR_effpot_FW0TT
#define FW0TZ SMDR_effpot_FW0TZ
#define FW0WW SMDR_effpot_FW0WW
#define FW0ZZ SMDR_effpot_FW0ZZ
#define FWhWZ SMDR_effpot_FWhWZ
#define FZ0hh SMDR_effpot_FZ0hh
#define FZ0hW SMDR_effpot_FZ0hW
#define FZ0TT SMDR_effpot_FZ0TT
#define FZ0TW SMDR_effpot_FZ0TW
#define FZ0WW SMDR_effpot_FZ0WW
#define FZ0WZ SMDR_effpot_FZ0WZ
#define FZ0ZZ SMDR_effpot_FZ0ZZ
#define FZhTT SMDR_effpot_FZhTT
#define FZhWW SMDR_effpot_FZhWW
#define FBAR00hT SMDR_effpot_FBAR00hT
#define FBAR00hW SMDR_effpot_FBAR00hW
#define FBAR00hZ SMDR_effpot_FBAR00hZ
#define FBAR00TW SMDR_effpot_FBAR00TW
#define FBAR00TZ SMDR_effpot_FBAR00TZ
#define FBAR00WZ SMDR_effpot_FBAR00WZ
#define FBAR0hTT SMDR_effpot_FBAR0hTT
#define FBAR0hWW SMDR_effpot_FBAR0hWW
#define FBAR0TTZ SMDR_effpot_FBAR0TTZ
#define FBAR0WWZ SMDR_effpot_FBAR0WWZ

SMDR_REAL lnbarh, lnbarT, lnbarT2, lnbarT3, lnbarT4;
SMDR_REAL AW, AZ, Ah, AT, Ab, Atau;

SMDR_REAL I0hW, I0hZ, I0TW, I0WZ, Ihhh, IhTT, IhWW, IhZZ, ITTZ, IWWZ;

SMDR_COMPLEX H00000h, H00000T, H00000W, H00000Z, H0000hW, H0000TT,
             H0000TW, H0000WW, H0000WZ, H000hTT, H000hWW, H000hZZ, 
             H000TTZ, H000WWZ, H00h0WW, H00hhWW, H00hhZZ, H00hW00, 
             H00hW0Z, H00hWWZ, H00hZ00, H00hZ0W, H00hZWW, H00T0TW, 
             H00T0WW, H00ThTT, H00TW0T, H00TZ00, H00TZ0W, H00TZWW, 
             H00W0WZ, H00Whhh, H00WhZZ, H00WThT, H00WTTZ, H00WW00, 
             H00WW0h, H00WW0W, H00WW0Z, H00WWhW, H00WWWZ, H00WZ00, 
             H00WZ0h, H00WZhZ, H00WZTT, H00Zhhh, H00ZhWW, H00ZWhW, 
             H00ZZ00, H00ZZ0W, H00ZZWW, H0hhWhW, H0hhZhZ, H0hTZTT, 
             H0hWWWh, H0hWWWZ, H0hZWZW, H0hZZZh, H0TTT0T, H0TTThT, 
             H0TTTZT, H0TTW0W, H0TTWhW, H0TTWZW, H0WWW0W, H0WWWhW, 
             H0WWWZW, H0WWZhZ, H0WZZWW, Hhhhhhh, HhhThTT, HhhWhWW, 
             HhhZhZZ, HhTTThT, HhTTTZT, HhTZTTZ, HhWWWhW, HhWWWZW, 
             HhWZWWZ, HhZZZhZ, HTTZZTT, HWWZZWW;

SMDR_COMPLEX G000hT, G000TZ, G0TTWW, Gh000W, Gh000Z, Gh00hh, Gh00TT, 
             Gh00WW, Gh00ZZ, Gh0W0Z, Gh0Whh, Gh0WTT, Gh0WWW, Gh0WZZ, 
             Gh0Zhh, Gh0ZTT, Gh0ZWW, Gh0ZZZ, Ghhhhh, GhhhTT, GhhhWW, 
             GhhhZZ, GhTTWW, GhTTZZ, GhWWZZ, GT000W, GT00hT, GT00TZ, 
             GT0WhT, GT0WTZ, GThTTZ, GW000Z, GW00hW, GW00WZ, GW0h0T, 
             GW0h0Z, GW0hhW, GW0hWZ, GW0T0Z, GW0ThW, GW0TWZ, GW0ZhW, 
             GW0ZWZ, GWhWWZ, GZ00hZ, GZ00TT, GZ00WW, GZ0h0W, GZ0hhZ, 
             GZ0hTT, GZ0hWW, GZ0WhZ, GZ0WTT, GZ0WWW, GZhZTT, GZhZWW, GZTTWW;

SMDR_COMPLEX Fh00T, Fh00W, Fh00Z, Fh0hW, Fh0hZ, Fh0TT, Fh0TW, Fh0WW, 
             Fh0WZ, Fh0ZZ, FhhTT, FhhWW, FhhZZ, FhTTZ, FhWWZ, FT00W, 
             FT00Z, FT0hW, FT0TW, FT0WZ, FThTZ, FTTWW, FTTZZ, FW00Z, 
             FW0hh, FW0hT, FW0hZ, FW0TT, FW0TZ, FW0WW, FW0ZZ, FWhWZ, 
             FZ0hh, FZ0hW, FZ0TT, FZ0TW, FZ0WW, FZ0WZ, FZ0ZZ, FZhTT, FZhWW;

SMDR_COMPLEX FBAR00hT, FBAR00hW, FBAR00hZ, FBAR00TW, FBAR00TZ, 
             FBAR00WZ, FBAR0hTT, FBAR0hWW, FBAR0TTZ, FBAR0WWZ;

/* ----------------------------------------------------------------- */
/* 
   Returns the minimum value of the effective potential (with Goldstone
   boson resummation). The minimum value returned is the one for which m2 
   is adjusted to assure that the resummed effective potential is minimized
   at the input VEV v.
*/

SMDR_REAL SMDR_Eval_Veffmin (SMDR_REAL Q_eval, float loopOrder) 
{
  SMDR_REAL Veffminresult;
  SMDR_REAL V10, V2hat0, V3hat0, V3hat0leading, V40;
  char funcname[] = "SMDR_Eval_Veffmin";

  /* Coefficients appearing in V3hat0 */
  SMDR_REAL co1, coAh, coAhAh, coAhAhAh, coAhAhAT, coAhAhAW, coAhAhAZ,
    coAhAT, coAhATAT, coAhATAW, coAhATAZ, coAhAW, coAhAWAW, coAhAWAZ,
    coAhAZ, coAhAZAZ, coAhI0hW, coAhI0hZ, coAhI0TW, coAhI0WZ, coAhIhhh,
    coAhIhTT, coAhIhWW, coAhIhZZ, coAhITTZ, coAhIWWZ, coAT, coATAT,
    coATATAT, coATATAW, coATATAZ, coATAW, coATAWAW, coATAWAZ, coATAZ,
    coATAZAZ, coATI0hW, coATI0hZ, coATI0TW, coATI0WZ, coATIhhh,
    coATIhTT, coATIhWW, coATIhZZ, coATITTZ, coATIWWZ, coAW, coAWAW,
    coAWAWAW, coAWAWAZ, coAWAZ, coAWAZAZ, coAWI0hW, coAWI0hZ, coAWI0TW,
    coAWI0WZ, coAWIhhh, coAWIhTT, coAWIhWW, coAWIhZZ, coAWITTZ,
    coAWIWWZ, coAZ, coAZAZ, coAZAZAZ, coAZI0hW, coAZI0hZ, coAZI0TW,
    coAZI0WZ, coAZIhhh, coAZIhTT, coAZIhWW, coAZIhZZ, coAZITTZ,
    coAZIWWZ, coFBAR00TW, coFBAR0hTT, coFBAR0hWW, coFBAR0TTZ, coFBAR0WWZ,
    coFh00W, coFh00Z, coFh0TT, coFh0TW, coFh0WW, coFh0ZZ, coFhhTT, coFhhWW,
    coFhhZZ, coFhTTZ, coFhWWZ, coFT00W, coFT0hW, coFT0WZ, coFThTZ, coFTTWW,
    coFTTZZ, coFW00Z, coFW0hT, coFW0TZ, coFWhWZ, coFZ0TT, coFZ0TW, coFZ0WW,
    coFZhTT, coFZhWW, coG0TTWW, coGhhhhh, coGhhhTT, coGhhhWW, coGhhhZZ,
    coGhTTWW, coGhTTZZ, coGhWWZZ, coGT000W, coGT0WhT, coGT0WTZ, coGThTTZ,
    coGW00hW, coGW00WZ, coGW0ThW, coGW0TWZ, coGWhWWZ, coGZ00hZ, coGZ00TT,
    coGZ00WW, coGZhZTT, coGZhZWW, coGZTTWW, coH00000W, coH00000Z, coH0000WW,
    coH000WWZ, coH00T0TW, coH00T0WW, coH00TZWW, coH00WZ00, coH00WZTT,
    coH00ZZ00, coH0TTT0T, coH0TTThT, coH0TTTZT, coH0TTW0W, coH0TTWhW,
    coH0TTWZW, coH0WWW0W, coH0WWWhW, coH0WWWZW, coHhhhhhh, coHhhThTT,
    coHhhWhWW, coHhhZhZZ, coHhTTThT, coHhTTTZT, coHhTZTTZ, coHhWWWhW,
    coHhWWWZW, coHhWZWWZ, coHhZZZhZ, coHTTZZTT, coHWWZZWW, coI0hW, coI0hZ,
    coI0TW, coI0WZ, coIhhh, coIhTT, coIhWW, coIhZZ, coITTZ, coIWWZ,
    coZeta2, coZeta2Ah, coZeta2AT, coZeta2AW, coZeta2AZ, coZeta3;

  if ( (TSIL_FABS(loopOrder) > 0.0001) &&
       (TSIL_FABS(loopOrder-1) > 0.0001) &&
       (TSIL_FABS(loopOrder-2) > 0.0001) &&
       (TSIL_FABS(loopOrder-2.5) > 0.0001) &&
       (TSIL_FABS(loopOrder-3) > 0.0001) &&
       (TSIL_FABS(loopOrder-3.5) > 0.0001) )
    SMDR_Error (funcname,
    "Invalid loopOrder specified, should be 0, 1, 2, 2.5, 3, or 3.5", 3);

  /* If Q_eval is negative, then we just use the current Q and working
     parameters. Otherwise, we run all the input parameters from Q_in 
     to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

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

  m2 = SMDR_Eval_m2 (-1, loopOrder);

  /* The following Update is called only because m2 has change, 
     and should really only affect G and H. Note that SMDR_Eval_m2() 
     has just evaluated all of the effective potential loop integrals, 
     and changing m2 doesn't affect any of them, because they all involve 
     h, not H, and they don't involve G ever. So it probably isn't actually 
     needed at all. But it isn't the bottleneck, and can't hurt. */
  SMDR_Update (); 

  Veffminresult = Lambda + (m2/2.0) * v2 + (k/4.0) * v4;

  if (loopOrder > 0.999) {   
    V10 = 3*AW*g2*v2/8 + 3*AZ*g2pgp2*v2/16 + Ah*k*v2/2 
          - 3*AT*v2*yt2/2 - 3*Ab*v2*yb2/2 - Atau*v2*ytau2/2 + 
          (v4*(3*g4 + 2*g2*gp2 + gp4 - 64*k2 + 48*yt4))/128;

    Veffminresult += ONELOOPFACTOR * V10;
  }

  /* From 1709.02397, ancillary file SMVresummedGexp.anc */
  if (loopOrder > 1.999) {    
    V2hat0 = (AW*AZ*(15*g4 - 10*g2*gp2 - gp4))/(4*g2pgp2) - 
    (AT*AZ*(9*g4 - 6*g2*gp2 + 17*gp4))/(6*g2pgp2) + 
    AZ*AZ*((5*(3*g4 + 5*gp4))/(6*g2pgp2) + k/2) - (3*Ah*Ah*k)/4 - 2*Ah*AW*k - 
    Ah*AZ*k + AW*AW*((49*g4 + 2*g2*gp2 + gp4)/(8*g2pgp2) + k) + 
    (Ah*(3*g4 + 2*g2*gp2 + gp4)*v2)/8 + 
    ((3*g2 - gp2)*(33*g4 + 22*g2*gp2 + gp4)*IWWZ*v2)/(32*g2pgp2) + 
    IhWW*((-3*g4)/8 + g2*k - 2*k2)*v2 + 
    IhZZ*((-3*g2pgp22)/16 + (g2pgp2*k)/2 - k2)*v2 - 3*Ihhh*k2*v2 + 
    AT*AW*(-3*g2 + 3*yt2) + 
    AT*AT*(24*g32 + (9*g4 + 90*g2*gp2 + 17*gp4)/(12*g2pgp2) + 3*yt2) + 
    AW*v2*(-(g2*(-29*g4 - 104*g2*gp2 - 3*gp4))/(24*g2pgp2) + g2*k - 
    (3*g2*yt2)/2) + ITTZ*v2*((-9*g4 + 6*g2*gp2 - 17*gp4)/48 + 
    ((9*g4 + 66*g2*gp2 - 7*gp4)*yt2)/(24*g2pgp2)) + 
    AT*v2*(-8*g32*yt2 - ((27*g4 + 92*g2*gp2 + gp4)*yt2)/(12*g2pgp2)) + 
    AZ*v2*(-(-30*g6 + 69*g4*gp2 + 84*g2*gp4 + 93*gp6)/(72*g2pgp2) + 
    (g2pgp2*k)/2 - ((9*g4 - 6*g2*gp2 + 17*gp4)*yt2)/(12*g2pgp2)) + 
    I0TW*v2*((-3*g4)/4 + (3*g2*yt2)/4 + (3*yt4)/2) + 
    IhTT*v2*(-3*k*yt2 + 3*yt4) + v4*((279*g8 + 144*g6*gp2 + 28*g4*gp4 + 
    56*g2*gp6 + 37*gp8)/(192*g2pgp2) + 
    ((-3*g4 - 2*g2*gp2 - gp4)*k)/8 + ((3*g2 - gp2)*gp2*yt2)/12 + 8*g32*yt4 + 
    ((18*g4 + 91*g2*gp2 + 9*gp4)*yt4)/(24*g2pgp2)) + 
    ((171*g6 + 69*g4*gp2 + 109*g2*gp4 + 103*gp6)*v4*Zeta2)/192;

    Veffminresult += TWOLOOPFACTOR * V2hat0;
  }

  if ((loopOrder > 2.4) && (loopOrder < 2.6)) {    
    /* From 1406.2355 [equation (4.32) + v^2/2 equation (4.21)]. */
    V3hat0leading = T * T * (
      g34 * (1957.3316481591396 - 1842.2008457655597 * lnbarT 
             + 868. * lnbarT2 - 184. * lnbarT3) + 
      g32 * yt2 * (-780.3020832735587 + 1220.9289095713748 * lnbarT 
                   - 360. * lnbarT2 + 60. * lnbarT3) +
      yt4 * (265.2930347126314 - 696.4100650178304 * lnbarT + 63. * lnbarT2 
             + 36. * lnbarh * lnbarT2 - 29.25 * lnbarT3));

    Veffminresult += THREELOOPFACTOR * V3hat0leading;
  } else
  if (loopOrder > 2.999) {    
    /* From 1709.02397, ancillary file SMVresummedGexp.anc */
    #include "includes/V3hat0.c"

    Veffminresult += THREELOOPFACTOR * V3hat0;
  }

  if (loopOrder > 3.001) {    
    /* From 1508.00912 equation (5.3) */
    V40 = g36 * T * T * (59366.97 - 54056.36 * lnbarT + 
          27699.06 * lnbarT2 - 9144. * lnbarT3 + 1380. * lnbarT4);  
    Veffminresult += FOURLOOPFACTOR * V40;
  }

  return Veffminresult;
}

/* ----------------------------------------------------------------- */

SMDR_REAL SMDR_Eval_m2 (SMDR_REAL Q_eval, float loopOrder) 
{
  char funcname[] = "SMDR_Eval_m2";
  SMDR_REAL Delta;

  if ( (TSIL_FABS(loopOrder) > 0.0001) &&
       (TSIL_FABS(loopOrder-1) > 0.0001) &&
       (TSIL_FABS(loopOrder-2) > 0.0001) &&
       (TSIL_FABS(loopOrder-2.5) > 0.0001) &&
       (TSIL_FABS(loopOrder-3) > 0.0001) &&
       (TSIL_FABS(loopOrder-3.5) > 0.0001) )
    SMDR_Error (funcname,
    "Invalid loopOrder specified, should be 0, 1, 2, 2.5, 3, or 3.5", 3);

  /* If Q_eval is negative, then we just use the current Q and working
     parameters. Otherwise, we run all the input parameters from Q_in
     to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

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

  Delta = SMDR_Eval_vevDelta (-1, loopOrder);

  return (-k * v2 - Delta);
}

/* ----------------------------------------------------------------- */

SMDR_REAL SMDR_Eval_vev (SMDR_REAL Q_eval, float loopOrder)
{
  SMDR_REAL oldv, lastv, newv;
  int i;
  int converged = 0;
  char funcname[] = "SMDR_Eval_vev";

  if ( (TSIL_FABS(loopOrder) > 0.0001) &&
       (TSIL_FABS(loopOrder-1) > 0.0001) &&
       (TSIL_FABS(loopOrder-2) > 0.0001) &&
       (TSIL_FABS(loopOrder-2.5) > 0.0001) &&
       (TSIL_FABS(loopOrder-3) > 0.0001) &&
       (TSIL_FABS(loopOrder-3.5) > 0.0001) )
    SMDR_Error (funcname,
    "Invalid loopOrder specified, should be 0, 1, 2, 2.5, 3, or 3.5", 3);

  /* If Q_eval is negative, then we just use the current Q and working
     parameters. Otherwise, we run all the input parameters from Q_in
     to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

  /* Check input parameters for sanity: */
  SMDR_Check_lambda_Range (0.05, 0.2);
  SMDR_Check_m2_Range (-10000,-2500);
  SMDR_Check_g3_Range (0.5, 1.5);
  SMDR_Check_g_Range (0.4, 0.8);
  SMDR_Check_gp_Range (0.25, 0.5);
  SMDR_Check_yt_Range (0.5, 1.5);
  SMDR_Check_yb_Range (0.05);
  SMDR_Check_ytau_Range (0.02);
  SMDR_Check_Q_Range (49, 501);

  /* Save original value. */
  oldv = v;

  /* Make an initial guess, in case v wasn't set yet. */
  v = 247.;

  /* Then iterate to find v2... */
  for (i = 0; i < MAXITERS_FINDVEV; i++) {

    lastv = v;
    v = newv = TVIL_SQRT(-(m2 + SMDR_Eval_vevDelta (-1, loopOrder))/k);
     
    /*
    printf("SMDR_Eval_vev iteration %d:  v = %.6lf;\n", i, (double) v);
    */

    /* Check for convergence */
    if (TVIL_FABS((v - lastv)/lastv) < TOLERANCE_FINDVEV) {
      converged = 1;
      break;
    }
  }

  if (0 == converged)
    SMDR_Warn (funcname, "Didn't converge after MAXITERS_FINDVEV iterations!");

  /* Restore original values. */
  v = oldv;
  SMDR_Update ();

  return newv;
}

/* ----------------------------------------------------------------- */

SMDR_REAL SMDR_Eval_vevDelta (SMDR_REAL Q_eval, float loopOrder) 
{
  SMDR_REAL Delta1, Delta2, Delta3, Delta3leading, Delta4;
  SMDR_REAL Deltaresult = 0;

  /* Coefficients appearing in Delta3: */  
  SMDR_REAL cod1, codAh, codAhAh, codAhAhAh, codAhAhAT, codAhAhAW, 
    codAhAhAZ, codAhAT, codAhATAT, codAhATAW, codAhATAZ, codAhAW, codAhAWAW, 
    codAhAWAZ, codAhAZ, codAhAZAZ, codAhI0hW, codAhI0hZ, codAhI0TW, 
    codAhI0WZ, codAhIhhh, codAhIhTT, codAhIhWW, codAhIhZZ, codAhITTZ, 
    codAhIWWZ, codAT, codATAT, codATATAT, codATATAW, codATATAZ, codATAW, 
    codATAWAW, codATAWAZ, codATAZ, codATAZAZ, codATI0hW, codATI0hZ, 
    codATI0TW, codATI0WZ, codATIhhh, codATIhTT, codATIhWW, codATIhZZ, 
    codATITTZ, codATIWWZ, codAW, codAWAW, codAWAWAW, codAWAWAZ, codAWAZ, 
    codAWAZAZ, codAWI0hW, codAWI0hZ, codAWI0TW, codAWI0WZ, codAWIhhh, 
    codAWIhTT, codAWIhWW, codAWIhZZ, codAWITTZ, codAWIWWZ, codAZ, 
    codAZAZ, codAZAZAZ, codAZI0hW, codAZI0hZ, codAZI0TW, codAZI0WZ, 
    codAZIhhh, codAZIhTT, codAZIhWW, codAZIhZZ, codAZITTZ, codAZIWWZ, 
    codFBAR00hT, codFBAR00hW, codFBAR00hZ, codFBAR00TW, codFBAR00TZ, 
    codFBAR00WZ, codFBAR0hTT, codFBAR0hWW, codFBAR0TTZ, codFBAR0WWZ, 
    codFh00T, codFh00W, codFh00Z, codFh0hW, codFh0hZ, codFh0TT, codFh0TW, 
    codFh0WW, codFh0WZ, codFh0ZZ, codFhhTT, codFhhWW, codFhhZZ, codFhTTZ, 
    codFhWWZ, codFT00W, codFT00Z, codFT0hW, codFT0TW, codFT0WZ, codFThTZ, 
    codFTTWW, codFTTZZ, codFW00Z, codFW0hh, codFW0hT, codFW0hZ, codFW0TT, 
    codFW0TZ, codFW0WW, codFW0ZZ, codFWhWZ, codFZ0hh, codFZ0hW, codFZ0TT, 
    codFZ0TW, codFZ0WW, codFZ0WZ, codFZ0ZZ, codFZhTT, codFZhWW, codG000hT, 
    codG000TZ, codG0TTWW, codGh000W, codGh000Z, codGh00hh, codGh00TT, 
    codGh00WW, codGh00ZZ, codGh0W0Z, codGh0Whh, codGh0WTT, codGh0WWW, 
    codGh0WZZ, codGh0Zhh, codGh0ZTT, codGh0ZWW, codGh0ZZZ, codGhhhhh, 
    codGhhhTT, codGhhhWW, codGhhhZZ, codGhTTWW, codGhTTZZ, codGhWWZZ, 
    codGT000W, codGT00hT, codGT00TZ, codGT0WhT, codGT0WTZ, codGThTTZ, 
    codGW000Z, codGW00hW, codGW00WZ, codGW0h0T, codGW0h0Z, codGW0hhW, 
    codGW0hWZ, codGW0T0Z, codGW0ThW, codGW0TWZ, codGW0ZhW, codGW0ZWZ, 
    codGWhWWZ, codGZ00hZ, codGZ00TT, codGZ00WW, codGZ0h0W, codGZ0hhZ, 
    codGZ0hTT, codGZ0hWW, codGZ0WhZ, codGZ0WTT, codGZ0WWW, codGZhZTT, 
    codGZhZWW, codGZTTWW, codH00000h, codH00000T, codH00000W, codH00000Z, 
    codH0000hW, codH0000TT, codH0000TW, codH0000WW, codH0000WZ, codH000hTT, 
    codH000hWW, codH000hZZ, codH000TTZ, codH000WWZ, codH00h0WW, codH00hhWW, 
    codH00hhZZ, codH00hW00, codH00hW0Z, codH00hWWZ, codH00hZ00, codH00hZ0W, 
    codH00hZWW, codH00T0TW, codH00T0WW, codH00ThTT, codH00TW0T, codH00TZ00, 
    codH00TZ0W, codH00TZWW, codH00W0WZ, codH00Whhh, codH00WhZZ, codH00WThT, 
    codH00WTTZ, codH00WW00, codH00WW0h, codH00WW0W, codH00WW0Z, codH00WWhW, 
    codH00WWWZ, codH00WZ00, codH00WZ0h, codH00WZhZ, codH00WZTT, codH00Zhhh, 
    codH00ZhWW, codH00ZWhW, codH00ZZ00, codH00ZZ0W, codH00ZZWW, codH0hhWhW, 
    codH0hhZhZ, codH0hTZTT, codH0hWWWh, codH0hWWWZ, codH0hZWZW, codH0hZZZh, 
    codH0TTT0T, codH0TTThT, codH0TTTZT, codH0TTW0W, codH0TTWhW, codH0TTWZW, 
    codH0WWW0W, codH0WWWhW, codH0WWWZW, codH0WWZhZ, codH0WZZWW, codHhhhhhh, 
    codHhhThTT, codHhhWhWW, codHhhZhZZ, codHhTTThT, codHhTTTZT, codHhTZTTZ, 
    codHhWWWhW, codHhWWWZW, codHhWZWWZ, codHhZZZhZ, codHTTZZTT, codHWWZZWW, 
    codI0hW, codI0hZ, codI0TW, codI0WZ, codIhhh, codIhTT, codIhWW, 
    codIhZZ, codITTZ, codIWWZ, codZeta2, codZeta2Ah, codZeta2AT, 
    codZeta2AW, codZeta2AZ, codZeta3;

  /* If Q_eval is negative, then we just use the current Q and working
     parameters. Otherwise, we run all the input parameters from Q_in
     to Q_eval. */
  if (Q_eval > 0) {
    SMDR_RGeval_SM (Q_eval, 5);
  }

  SMDR_Update ();
  SMDR_effpot_Do3VILv (loopOrder);

  if (loopOrder > 0.999) {   
    /* From 1406.2355 eq. (4.19), or SMDeltas.anc of 1709.02397 */
    Delta1 = ((3*g4 + 2*g2*gp2 + gp4)*v2)/8 + 3*k*Ah 
             - 6*yt2*AT - 6*yb2*Ab - 2*ytau2*Atau 
             + (3*g2*AW)/2 + (3*(g2 + gp2)*AZ)/4;

    Deltaresult += ONELOOPFACTOR * Delta1;
  }

  if (loopOrder > 1.999) {    
    /* From 1406.2355 eq. (4.20), or SMDeltas.anc of 1709.02397 */
    Delta2 = 
    (-2*AT*AZ*(9*g4 - 6*g2*gp2 + 17*gp4))/(3*g2pgp2*v2) + 
    (Ah*AW*((3*g2)/2 + (3*g4)/(4*g2m2k) - 8*k))/v2 + 
    (Ah*AZ*((3*g2pgp2)/4 + (3*g2pgp22)/(8*g2pgp2m2k) - 4*k))/v2 + 
    (3*Ah*Ah*k)/(2*v2) + 
    (AW*AW*((-3*g4)/(4*g2m2k) + (49*g4 + 2*g2*gp2 + gp4)/(2*g2pgp2) + 
    ((g4 - 10*g2*gp2 - 8*gp4)*k)/g2pgp22))/v2 + 
    (AW*AZ*((15*g4 - 10*g2*gp2 - gp4)/g2pgp2 + (6*(8*g4 + 8*g2*gp2 + gp4)*k)/
    g2pgp22))/v2 + 
    (AZ*AZ*((-3*g2pgp22)/(8*g2pgp2m2k) + (10*(3*g4 + 5*gp4))/(3*g2pgp2) + 
    ((g4 + 26*g2*gp2 + 7*gp4)*k)/(2*g2pgp22)))/v2 - (9*Ah*AT*yt2)/v2 + 
    (AT*AW*(-12*g2 + 12*yt2))/v2 + 
    (AT*AT*(96*g32 + (9*g4 + 90*g2*gp2 + 17*gp4)/(3*g2pgp2) + 18*k + 15*yt2))/
    v2 + AW*((-3*g6)/(8*g2m2k) - (g2*(-101*g4 - 428*g2*gp2 - 39*gp4))/
    (24*g2pgp2) - (g2*(16*g4 + 14*g2*gp2 + gp4)*k)/(2*g2pgp22) - 
    (21*g2*yt2)/2) + AZ*((-3*g2pgp23)/(16*g2pgp2m2k) - 
    (-195*g6 + 165*g4*gp2 - 633*g2*gp4 - 129*gp6)/(144*g2pgp2) - 
    ((20*g4 + 28*g2*gp2 - gp4)*k)/(2*g2pgp2) - 
    ((63*g4 + 30*g2*gp2 + 95*gp4)*yt2)/(12*g2pgp2)) + 
    AT*(16*g32*yt2 - ((27*g4 + 80*g2*gp2 - 11*gp4)*yt2)/(3*g2pgp2) + 12*k*yt2 - 
    12*yt4) + Ah*((-3*g2pgp23)/(32*g2pgp2m2k) - (3*g6)/(16*g2m2k) + 
    (3*(3*g4 + 2*g2*gp2 + gp4))/2 + 18*k2 - (7*k*threeg2pgp2)/4 + 9*k*yt2 - 
    (21*yt4)/2) + v2*((9*g2pgp24)/(64*g2pgp2m2k) + (9*g8)/(32*g2m2k) + 
    (5571*g8 + 3018*g6*gp2 - 1836*g4*gp4 - 1266*g2*gp6 - 255*gp8)/(576*g2pgp2) + 
    ((81*g8 + 158*g6*gp2 + 110*g4*gp4 + 28*g2*gp6 + gp8)*k)/(16*g2pgp22) - 
    48*k3 + 2*k2*threeg2pgp2 + ((-81*g4 + 138*g2*gp2 - 91*gp4 - 864*k2)*yt2)/
    48 + 24*g32*yt4 + ((18*g4 + 87*g2*gp2 + 5*gp4)*yt4)/(6*g2pgp2) + 12*k*yt4 + 
    9*yt6) + v2*((171*g6 + 69*g4*gp2 + 109*g2*gp4 + 103*gp6)/48 - 
    (3*g2*(3*g6 + 8*g4*gp2 + 2*g2*gp4 - 2*gp6)*k)/(8*g2pgp22) + 12*k3 + 
    6*k*yt4)*Zeta2 + ((3*g2*k)/2 + 12*k2)*I0hW + ((3*g2pgp2*k)/4 + 6*k2)*
    I0hZ + (-3*g4 + 3*g2*yt2 + 6*yt4)*I0TW + 
    (3*k*twog2pgp23*I0WZ)/(2*g2pgp22) - 15*k2*Ihhh + 
    (-18*k*yt2 + (27*yt4)/2)*IhTT + 
    ((-15*g4)/8 + (3*g6)/(16*g2m2k) + (11*g2*k)/2 - 14*k2)*IhWW + 
    ((-15*g2pgp22)/16 + (3*g2pgp23)/(32*g2pgp2m2k) + (11*g2pgp2*k)/4 - 7*k2)*
    IhZZ + ((-9*g4 + 6*g2*gp2 - 17*gp4)/12 + 
    ((9*g4 + 66*g2*gp2 - 7*gp4)*yt2)/(6*g2pgp2))*ITTZ + 
    ((33*g4 + 22*g2*gp2 + gp4)*threeg2mgp2*IWWZ)/(8*g2pgp2);

    Deltaresult += TWOLOOPFACTOR * Delta2;
  }

  if ((loopOrder > 2.4) && (loopOrder < 2.6)) {    
    /* From 1406.2355 equation (4.21). */
    Delta3leading = v2 * yt4 * (
      g34 * (1036.23 - 974.20 * lnbarT + 592 * lnbarT2 - 184 * lnbarT3) + 
      g32 * yt2 * (-169.84 + 860.93 * lnbarT - 270 * lnbarT2 + 60 * lnbarT3) +
      yt4 * (-82.91 - 753.02 * lnbarT + 36 * lnbarh * lnbarT + 
      82.125 * lnbarT2 + 54 * lnbarh * lnbarT2 - 56.25 * lnbarT3));

    /* Should be the same: */
    Delta3leading = v2 *  yt4 * (g34 * (1036.2312252763598258864
      - 974.200845765559649533875*lnbarT + 592.*lnbarT2 - 184.*lnbarT3)
      + g32 * yt2 * (-169.83762848787117086562
      + 860.928909571374931921951*lnbarT - 270.*lnbarT2 + 60.*lnbarT3)
      + yt4 * (-82.91199779628403181268
      - 753.018878221098564951523*lnbarT + 36.*lnbarh*lnbarT + 82.125*lnbarT2
      + 54.*lnbarh*lnbarT2 - 56.25*lnbarT3));

    Deltaresult += THREELOOPFACTOR * Delta3leading;
  } else
  if (loopOrder > 2.999) {    
    /* From 1709.02397, ancillary file SMDeltas.anc. */
    #include "includes/V3SMDelta3.c"

    Deltaresult += THREELOOPFACTOR*Delta3;
  }

  if (loopOrder > 3.001) {    
    /* From 1508.00912 equation (5.5) */
    Delta4 = g36 * yt2 * T * (64677.58 - 52714.59 * lnbarT + 
             27966.13 * lnbarT2 - 12768 * lnbarT3 + 2760 * lnbarT4);  
    Deltaresult += FOURLOOPFACTOR*Delta4;
  }

  return TVIL_CREAL(Deltaresult);
}

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

int SMDR_effpot_Do3VILv (float loopOrder) {
  int success = 1;
  TVIL_RESULT result;

  /* Must be done before perpetrating 3VIL: */
  /* Actually no, this is obsolete?
  TVIL_Setup ();
  */

  if (loopOrder > 0.999) {
    AW = TVIL_CREAL(TVIL_A(W,Q2));
    AZ = TVIL_CREAL(TVIL_A(Z,Q2));
    Ah = TVIL_CREAL(TVIL_A(h,Q2));
    AT = TVIL_CREAL(TVIL_A(T,Q2));
    Ab = TVIL_CREAL(TVIL_A(b,Q2));
    Atau = TVIL_CREAL(TVIL_A(tau,Q2));
    lnbarh = Ah/h + 1;
    lnbarT = AT/T + 1;
    lnbarT2 = lnbarT * lnbarT;
    lnbarT3 = lnbarT2 * lnbarT;
    lnbarT4 = lnbarT3 * lnbarT;
  }

  if (loopOrder > 1.999) {
    I0hW = TVIL_CREAL(TVIL_I2 (0, h, W, Q2));
    I0hZ = TVIL_CREAL(TVIL_I2 (0, h, Z, Q2));
    I0TW = TVIL_CREAL(TVIL_I2 (0, T, W, Q2));
    I0WZ = TVIL_CREAL(TVIL_I2 (0, W, Z, Q2));
    Ihhh = TVIL_CREAL(TVIL_I2 (h, h, h, Q2));
    IhTT = TVIL_CREAL(TVIL_I2 (h, T, T, Q2));
    IhWW = TVIL_CREAL(TVIL_I2 (h, W, W, Q2));
    IhZZ = TVIL_CREAL(TVIL_I2 (h, Z, Z, Q2));
    ITTZ = TVIL_CREAL(TVIL_I2 (T, T, Z, Q2));
    IWWZ = TVIL_CREAL(TVIL_I2 (W, W, Z, Q2));
  }

  if (loopOrder > 2.999) {
  /* Analytic cases. */

  success *= TVIL_Hanalytic (0, 0, 0, 0, 0, h, Q2, &H00000h);
  success *= TVIL_Hanalytic (0, 0, 0, 0, 0, T, Q2, &H00000T);
  success *= TVIL_Hanalytic (0, 0, 0, 0, 0, W, Q2, &H00000W);
  success *= TVIL_Hanalytic (0, 0, 0, 0, 0, Z, Q2, &H00000Z);
  success *= TVIL_Hanalytic (0, 0, 0, 0, h, W, Q2, &H0000hW);
  success *= TVIL_Hanalytic (0, 0, 0, 0, T, T, Q2, &H0000TT);
  success *= TVIL_Hanalytic (0, 0, 0, 0, T, W, Q2, &H0000TW);
  success *= TVIL_Hanalytic (0, 0, 0, 0, W, W, Q2, &H0000WW);
  success *= TVIL_Hanalytic (0, 0, 0, 0, W, Z, Q2, &H0000WZ);
  success *= TVIL_Hanalytic (0, 0, 0, h, T, T, Q2, &H000hTT);
  success *= TVIL_Hanalytic (0, 0, 0, h, W, W, Q2, &H000hWW);
  success *= TVIL_Hanalytic (0, 0, 0, h, Z, Z, Q2, &H000hZZ);
  success *= TVIL_Hanalytic (0, 0, 0, T, T, Z, Q2, &H000TTZ);
  success *= TVIL_Hanalytic (0, 0, 0, W, W, Z, Q2, &H000WWZ);
  success *= TVIL_Hanalytic (0, 0, h, 0, W, W, Q2, &H00h0WW);
  success *= TVIL_Hanalytic (0, 0, h, h, W, W, Q2, &H00hhWW);
  success *= TVIL_Hanalytic (0, 0, h, h, Z, Z, Q2, &H00hhZZ);
  success *= TVIL_Hanalytic (0, 0, h, W, 0, 0, Q2, &H00hW00);
  success *= TVIL_Hanalytic (0, 0, h, Z, 0, 0, Q2, &H00hZ00);
  success *= TVIL_Hanalytic (0, 0, T, 0, T, W, Q2, &H00T0TW);
  success *= TVIL_Hanalytic (0, 0, T, 0, W, W, Q2, &H00T0WW);
  success *= TVIL_Hanalytic (0, 0, T, h, T, T, Q2, &H00ThTT);
  success *= TVIL_Hanalytic (0, 0, T, W, 0, T, Q2, &H00TW0T);
  success *= TVIL_Hanalytic (0, 0, T, Z, 0, 0, Q2, &H00TZ00);
  success *= TVIL_Hanalytic (0, 0, W, 0, W, Z, Q2, &H00W0WZ);
  success *= TVIL_Hanalytic (0, 0, W, h, h, h, Q2, &H00Whhh);
  success *= TVIL_Hanalytic (0, 0, W, W, 0, 0, Q2, &H00WW00);
  success *= TVIL_Hanalytic (0, 0, W, W, 0, h, Q2, &H00WW0h);
  success *= TVIL_Hanalytic (0, 0, W, W, 0, W, Q2, &H00WW0W);
  success *= TVIL_Hanalytic (0, 0, W, W, 0, Z, Q2, &H00WW0Z);
  success *= TVIL_Hanalytic (0, 0, W, W, h, W, Q2, &H00WWhW);
  success *= TVIL_Hanalytic (0, 0, W, W, W, Z, Q2, &H00WWWZ);
  success *= TVIL_Hanalytic (0, 0, W, Z, 0, 0, Q2, &H00WZ00);
  success *= TVIL_Hanalytic (0, 0, Z, h, h, h, Q2, &H00Zhhh);
  success *= TVIL_Hanalytic (0, 0, Z, Z, 0, 0, Q2, &H00ZZ00);
  success *= TVIL_Hanalytic (0, 0, Z, Z, 0, W, Q2, &H00ZZ0W);
  success *= TVIL_Hanalytic (0, 0, Z, Z, W, W, Q2, &H00ZZWW);
  success *= TVIL_Hanalytic (0, h, h, W, h, W, Q2, &H0hhWhW);
  success *= TVIL_Hanalytic (0, h, h, Z, h, Z, Q2, &H0hhZhZ);
  success *= TVIL_Hanalytic (0, h, W, W, W, h, Q2, &H0hWWWh);
  success *= TVIL_Hanalytic (0, h, Z, Z, Z, h, Q2, &H0hZZZh);
  success *= TVIL_Hanalytic (0, T, T, T, 0, T, Q2, &H0TTT0T);
  success *= TVIL_Hanalytic (0, T, T, T, h, T, Q2, &H0TTThT);
  success *= TVIL_Hanalytic (0, T, T, T, Z, T, Q2, &H0TTTZT);
  success *= TVIL_Hanalytic (0, T, T, W, 0, W, Q2, &H0TTW0W);
  success *= TVIL_Hanalytic (0, W, W, W, 0, W, Q2, &H0WWW0W);
  success *= TVIL_Hanalytic (0, W, W, W, h, W, Q2, &H0WWWhW);
  success *= TVIL_Hanalytic (0, W, W, W, Z, W, Q2, &H0WWWZW);
  success *= TVIL_Hanalytic (0, W, Z, Z, W, W, Q2, &H0WZZWW);
  success *= TVIL_Hanalytic (h, h, h, h, h, h, Q2, &Hhhhhhh);

  success *= TVIL_Ganalytic (0, 0, 0, h, T, Q2, &G000hT);
  success *= TVIL_Ganalytic (0, 0, 0, T, Z, Q2, &G000TZ);
  success *= TVIL_Ganalytic (0, T, T, W, W, Q2, &G0TTWW);
  success *= TVIL_Ganalytic (h, 0, 0, h, h, Q2, &Gh00hh);
  success *= TVIL_Ganalytic (h, 0, 0, T, T, Q2, &Gh00TT);
  success *= TVIL_Ganalytic (h, 0, W, h, h, Q2, &Gh0Whh);
  success *= TVIL_Ganalytic (h, 0, Z, h, h, Q2, &Gh0Zhh);
  success *= TVIL_Ganalytic (h, 0, Z, Z, Z, Q2, &Gh0ZZZ);
  success *= TVIL_Ganalytic (h, h, h, h, h, Q2, &Ghhhhh);
  success *= TVIL_Ganalytic (W, 0, h, h, W, Q2, &GW0hhW);
  success *= TVIL_Ganalytic (Z, 0, h, h, Z, Q2, &GZ0hhZ);

  success *= TVIL_Fanalytic (h, 0, h, W, Q2, &Fh0hW);
  success *= TVIL_Fanalytic (h, 0, h, Z, Q2, &Fh0hZ);
  success *= TVIL_Fanalytic (h, 0, Z, Z, Q2, &Fh0ZZ);
  success *= TVIL_Fanalytic (T, 0, 0, W, Q2, &FT00W);
  success *= TVIL_Fanalytic (W, 0, h, h, Q2, &FW0hh);
  success *= TVIL_Fanalytic (Z, 0, h, h, Q2, &FZ0hh);
  success *= TVIL_Fanalytic (Z, 0, Z, Z, Q2, &FZ0ZZ);

  success *= TVIL_FBARanalytic (0, 0, h, T, Q2, &FBAR00hT);
  success *= TVIL_FBARanalytic (0, 0, h, W, Q2, &FBAR00hW);
  success *= TVIL_FBARanalytic (0, 0, h, Z, Q2, &FBAR00hZ);
  success *= TVIL_FBARanalytic (0, 0, T, W, Q2, &FBAR00TW);
  success *= TVIL_FBARanalytic (0, 0, T, Z, Q2, &FBAR00TZ);
  success *= TVIL_FBARanalytic (0, 0, W, Z, Q2, &FBAR00WZ);
  success *= TVIL_FBARanalytic (0, h, T, T, Q2, &FBAR0hTT);
  success *= TVIL_FBARanalytic (0, h, W, W, Q2, &FBAR0hWW);
  success *= TVIL_FBARanalytic (0, T, T, Z, Q2, &FBAR0TTZ);
  success *= TVIL_FBARanalytic (0, W, W, Z, Q2, &FBAR0WWZ);

  if (success != 1) 
    printf("Disaster! Some analytic function evaluation failed.\n");

  /*
  printf("Done with analytic 3-loop functions. success = %d\n",success);
  */

  /* Evaluate non-analytic 3-loop functions: */

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 444, h, W, 0, Z, Q2, &result);

  H00hW0Z = TVIL_GetFunction (&result, "H");
  GZ0h0W = TVIL_GetFunction (&result, "Gzuwxy");
  GW000Z = TVIL_GetFunction (&result, "Gxuvyz");
  Gh000Z = TVIL_GetFunction (&result, "Gwuzvy");
  FZ0hW = TVIL_GetFunction (&result, "Fzvwx");
  FW0hZ = TVIL_GetFunction (&result, "Fxvwz");
  Fh00W = TVIL_GetFunction (&result, "Fwuxy");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, h, Z, 0, W, Q2, &result);

  H00hZ0W = TVIL_GetFunction (&result, "H");
  Gh000W = TVIL_GetFunction (&result, "Gwuzvy");
  GW0h0Z = TVIL_GetFunction (&result, "Gzuwxy");
  Fh00Z = TVIL_GetFunction (&result, "Fwuxy");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, h, W, W, Z, Q2, &result);

  H00hWWZ = TVIL_GetFunction (&result, "H");
  Gh0W0Z = TVIL_GetFunction (&result, "Gwuzvy");
  GW00WZ = TVIL_GetFunction (&result, "Gxuvyz");
  GZ0hWW = TVIL_GetFunction (&result, "Gzuwxy");
  GW0hWZ = TVIL_GetFunction (&result, "Gyvwxz");
  Fh0WW = TVIL_GetFunction (&result, "Fwuxy");
  FW00Z = TVIL_GetFunction (&result, "Fyuvz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, h, Z, W, W, Q2, &result);

  H00hZWW = TVIL_GetFunction (&result, "H");
  GZ00WW = TVIL_GetFunction (&result, "Gxuvyz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, T, Z, 0, W, Q2, &result);

  H00TZ0W = TVIL_GetFunction (&result, "H");
  GT000W = TVIL_GetFunction (&result, "Gwuzvy");
  GW0T0Z = TVIL_GetFunction (&result, "Gzuwxy");
  FT00Z = TVIL_GetFunction (&result, "Fwuxy");
  FT0WZ = TVIL_GetFunction (&result, "Fwvxz");
  FW0TZ = TVIL_GetFunction (&result, "Fzvwx");
  FZ0TW = TVIL_GetFunction (&result, "Fxvwz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, T, Z, W, W, Q2, &result);

  H00TZWW = TVIL_GetFunction (&result, "H");
  GW0TWZ = TVIL_GetFunction (&result, "Gzuwxy");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, W, h, Z, Z, Q2, &result);

  H00WhZZ = TVIL_GetFunction (&result, "H");
  GZ0WhZ = TVIL_GetFunction (&result, "Gzuwxy");
  Gh00ZZ = TVIL_GetFunction (&result, "Gxuvyz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, W, T, h, T, Q2, &result);

  H00WThT = TVIL_GetFunction (&result, "H");
  GW0h0T = TVIL_GetFunction (&result, "Gwuzvy");
  GT00hT = TVIL_GetFunction (&result, "Gxuvyz");
  GT0WhT = TVIL_GetFunction (&result, "Gzuwxy");
  Gh0WTT = TVIL_GetFunction (&result, "Gyvwxz");
  FW0hT = TVIL_GetFunction (&result, "Fwuxy");
  FW0TT = TVIL_GetFunction (&result, "Fwvxz");
  FT0TW = TVIL_GetFunction (&result, "Fzvwx");
  Fh00T = TVIL_GetFunction (&result, "Fyuvz");
  Fh0TW = TVIL_GetFunction (&result, "Fyuwx");
  FT0hW = TVIL_GetFunction (&result, "Fxuwy");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, W, T, T, Z, Q2, &result);

  H00WTTZ = TVIL_GetFunction (&result, "H");
  GT00TZ = TVIL_GetFunction (&result, "Gxuvyz");
  GT0WTZ = TVIL_GetFunction (&result, "Gyvwxz");
  GZ0WTT = TVIL_GetFunction (&result, "Gzuwxy");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, W, Z, 0, h, Q2, &result);

  H00WZ0h = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, W, Z, h, Z, Q2, &result);

  H00WZhZ = TVIL_GetFunction (&result, "H");
  GZ00hZ = TVIL_GetFunction (&result, "Gxuvyz");
  Gh0WZZ = TVIL_GetFunction (&result, "Gyvwxz");
  FZ0WZ = TVIL_GetFunction (&result, "Fxvwz");
  FW0ZZ = TVIL_GetFunction (&result, "Fwvxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, W, Z, T, T, Q2, &result);

  H00WZTT = TVIL_GetFunction (&result, "H");
  GZ00TT = TVIL_GetFunction (&result, "Gxuvyz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, Z, h, W, W, Q2, &result);

  H00ZhWW = TVIL_GetFunction (&result, "H");
  Gh00WW = TVIL_GetFunction (&result, "Gxuvyz");
  GW0ZhW = TVIL_GetFunction (&result, "Gyvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, 0, Z, W, h, W, Q2, &result);

  H00ZWhW = TVIL_GetFunction (&result, "H");
  Gh0ZWW = TVIL_GetFunction (&result, "Gyvwxz");
  GW00hW = TVIL_GetFunction (&result, "Gxuvyz");
  FZ0WW = TVIL_GetFunction (&result, "Fwvxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, h, T, Z, T, T, Q2, &result);

  H0hTZTT = TVIL_GetFunction (&result, "H");
  Gh0ZTT = TVIL_GetFunction (&result, "Gvuxwy");
  GZ0hTT = TVIL_GetFunction (&result, "Gxuvyz");
  GThTTZ = TVIL_GetFunction (&result, "Gyvwxz");
  FThTZ = TVIL_GetFunction (&result, "Fwvxz");
  FhTTZ = TVIL_GetFunction (&result, "Fvwxz");
  FZhTT = TVIL_GetFunction (&result, "Fxwvz");
  Fh0TT = TVIL_GetFunction (&result, "Fvuyz");
  FZ0TT = TVIL_GetFunction (&result, "Fxuwy");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, h, W, W, W, Z, Q2, &result);

  H0hWWWZ = TVIL_GetFunction (&result, "H");
  Gh0WWW = TVIL_GetFunction (&result, "Gvuxwy");
  GWhWWZ = TVIL_GetFunction (&result, "Gyvwxz");
  GZ0WWW = TVIL_GetFunction (&result, "Gzuwxy");
  FWhWZ = TVIL_GetFunction (&result, "Fwvxz");
  FhWWZ = TVIL_GetFunction (&result, "Fvwxz");
  FZhWW = TVIL_GetFunction (&result, "Fzvwx");
  FW0WW = TVIL_GetFunction (&result, "Fxuwy");
  Fh0WZ = TVIL_GetFunction (&result, "Fvuyz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, h, Z, W, Z, W, Q2, &result);

  H0hZWZW = TVIL_GetFunction (&result, "H");
  GW0ZWZ = TVIL_GetFunction (&result, "Gzuwxy");
  GZhZWW = TVIL_GetFunction (&result, "Gyvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, T, T, W, h, W, Q2, &result);

  H0TTWhW = TVIL_GetFunction (&result, "H");
  GhTTWW = TVIL_GetFunction (&result, "Gyvwxz");
  GW0ThW = TVIL_GetFunction (&result, "Gzuwxy");
  FTTWW = TVIL_GetFunction (&result, "Fvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, T, T, W, Z, W, Q2, &result);

  H0TTWZW = TVIL_GetFunction (&result, "H");
  GZTTWW = TVIL_GetFunction (&result, "Gyvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (0, W, W, Z, h, Z, Q2, &result);

  H0WWZhZ = TVIL_GetFunction (&result, "H");
  GhWWZZ = TVIL_GetFunction (&result, "Gyvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, h, T, h, T, T, Q2, &result);

  HhhThTT = TVIL_GetFunction (&result, "H");
  GhhhTT = TVIL_GetFunction (&result, "Gvuxwy");
  FhhTT = TVIL_GetFunction (&result, "Fvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, h, W, h, W, W, Q2, &result);

  HhhWhWW = TVIL_GetFunction (&result, "H");
  GhhhWW = TVIL_GetFunction (&result, "Gvuxwy");
  FhhWW = TVIL_GetFunction (&result, "Fvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, h, Z, h, Z, Z, Q2, &result);

  HhhZhZZ = TVIL_GetFunction (&result, "H");
  GhhhZZ = TVIL_GetFunction (&result, "Gvuxwy");
  FhhZZ = TVIL_GetFunction (&result, "Fvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, T, T, T, h, T, Q2, &result);

  HhTTThT = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, T, T, T, Z, T, Q2, &result);

  HhTTTZT = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, T, Z, T, T, Z, Q2, &result);

  HhTZTTZ = TVIL_GetFunction (&result, "H");
  GhTTZZ = TVIL_GetFunction (&result, "Guvxwz");
  GZhZTT = TVIL_GetFunction (&result, "Gzuwxy");
  FTTZZ = TVIL_GetFunction (&result, "Fvwxz");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, W, W, W, h, W, Q2, &result);

  HhWWWhW = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, W, W, W, Z, W, Q2, &result);

  HhWWWZW = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, W, Z, W, W, Z, Q2, &result);

  HhWZWWZ = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (h, Z, Z, Z, h, Z, Q2, &result);

  HhZZZhZ = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (T, T, Z, Z, T, T, Q2, &result);

  HTTZZTT = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  TVIL_Evaluate (W, W, Z, Z, W, W, Q2, &result);

  HWWZZWW = TVIL_GetFunction (&result, "H");

  /* --------------------------------------------------------------- */

  /*
  printf("Done with non-analytic 3-loop functions.\n");
  */
  }

  return success;
}

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
