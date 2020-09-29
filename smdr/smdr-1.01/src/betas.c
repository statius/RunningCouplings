/* 
   Standard Model Beta functions at 3 loops and beyond. 
   Sources:
   1210.6873, 1212.6829, 1303.4364, Bednyakov, Pikelner, Velizhanin
   Checked in part with
   1205.2892, 1303.2890, Chetyrkin and Zoller,
   and 
   0111190 Ford, Jack, Jones.

   Included at 4-loop order:
     betak, leading in QCD (SPM 1508.00912 and Chetyrkin+Zoller 1604.00853)
     yt leading in QCD (Chetyrkin 9703278, Vermaseren, Larin, van Ritbergen
      9703284) 

   4-loop beta functions for pure QCD: 
     9701390 van Ritbergen, Vermaseren, Larin
     0405193 Czakon

   4-loop beta functions for QCD coupling, but including also contributions 
   from yt and lambda:
     1508.02680 Bednyakov and Pikelner
     1508.03624 Zoller 

   Remaining 4-loop beta functions for SM gauge couplings, from
   1912.07624 "Gauge coupling beta functions to four-loop
   order in the Standard Model", by Joshua Davies, Florian Herren, Colin
   Poole, Matthias Steinhauser, Anders Eller Thomsen.

   Five-loop pure QCD beta functions are given in 
   arXiv:1606.08659 Baikov, Chetyrkin, Kuhn.
   arXiv:1606.08662 Luthe, Maier, Marquard, Schroeder
   arXiv:1701.01404 Herzog, Ruijl, Ueda, Vermaseren, Vogt

   Five-loop QCD contributions to yt beta function are given in 
   arXiv:1402.6611

   kappa = ONELOOPFACTOR = 1/(16 Pi^2)
   Zeta3 = Zeta[3] = 1.2020569031595943

   Running parameters are:
   g3, g, gp = gauge couplings
   yt, yb, ytau = Yukawa couplings
   v  = VEV (in normalization where v is roughly 246 GeV)
   k = lambda = Higgs^4 coupling in normalization where mh^2 = 2 k v^2
   m2 = negative Higgs mass^2 parameter.

   BetaX = Q dX/dQ
   gamma = -Q d[ln(v)]/dQ

   Betag3 = kappa betag31 + kappa^2 betag32 + kappa^3 betag33 + ...;
   Betag = kappa betag1 + kappa^2 betag2 + kappa^3 betag3;
   Betagp = kappa betagp1 + kappa^2 betagp2 + kappa^3 betagp3;
   Betayt = kappa betayt1 + kappa^2 betayt2 + kappa^3 betayt3 + ...;
   Betayb = kappa betayb1 + kappa^2 betayb2 + kappa^3 betayb3;
   Betaytau = kappa betaytau1 + kappa^2 betaytau2 + kappa^3 betaytau3;
   Betak = kappa betak1 + kappa^2 betak2 + kappa^3 betak3 + ...;
   Betam2 = m2 (kappa betam21 + kappa^2 betam22 + kappa^3 betam23);
   gamma = kappa gamma1 + kappa^2 gamma2 + kappa^3 gamma3;
   BetaLambda = kappa betaLambda1 + kappa^2 betaLambda2 + kappa^3 betaLambda3;
*/

#include "smdr_internal.h"

/* See 1508.02680 and 1508.03624. The parameter Rgammafive reflects a
   gamma_5 ambiguity at 4-loop order. The value 3 apparently is
   apparently the correct one, as argued in 1508.02680 and later by
   Poole and Thomsen, in 1901.02749.
*/
#define Rgammafive 3

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* 
   Full 3-loop SM beta functions, with known 4-loop (partly
   implemented) and 5-loop pure QCD (not implemented yet but will be).
*/

int SMDR_Betas (int loopOrder)
{ 
  SMDR_REAL betag31, betag32, betag33, betag34, betag35;
  SMDR_REAL betag1, betag2, betag3, betag4;
  SMDR_REAL betagp1, betagp2, betagp3, betagp4;
  SMDR_REAL betayt1, betayt2, betayt3, betayt4, betayt5;
  SMDR_REAL betayb1, betayb2, betayb3;
  SMDR_REAL betaytau1, betaytau2, betaytau3;
  SMDR_REAL betak1, betak2, betak3, betak4;
  SMDR_REAL betam21, betam22, betam23;
  SMDR_REAL gamma, gamma1, gamma2, gamma3;
  SMDR_REAL betaLambda1,  betaLambda2,  betaLambda3;
  SMDR_REAL betayc1, betayu1, betays1, betayd1, betaymu1, betaye1;
  SMDR_REAL betayc2, betayu2, betays2, betayd2, betaymu2, betaye2;
  SMDR_REAL betayc3, betayu3, betays3, betayd3;
  SMDR_REAL Y2S, Y4S, chi4S;

  char funcname[] = "SMDR_Betas";

  if ( (loopOrder < 1) || (loopOrder > 5) )
    SMDR_Error (funcname,
    "Invalid loop order specified, should be 1, 2, 3, 4, or 5.", 3);

  SMDR_update_for_betas ();

  Y2S = 3 * (yt2 + yb2 + yc*yc + ys*ys + yu*yu + yd*yd) + 
        ytau2 + ymu*ymu + ye*ye;

  Y4S = yt2 * (8*g32 + (9*g2)/4 + (17*gp2)/12) +
        yb2 * (8*g32 + (9*g2)/4 + (5*gp2)/12) + 
        ytau2 * ((3*g2)/4 + (5*gp2)/4); 

  chi4S = (9./4.) * (3 * yt4 + 3 * yb4 + ytau4 - (2 * yt2 * yb2)/3);

  betag31 = -7*g3*g32;

  betag1 = (-19*g*g2)/6;

  betagp1 = (41*gp*gp2)/6;

  betayt1 = yt * (-8*g32 - (9*g2)/4 - (17*gp2)/12 + 3.*(yt2 - yb2)/2. + Y2S);

  betayb1 = yb * (-8*g32 - (9*g2)/4 - (5*gp2)/12 + 3.*(yb2 - yt2)/2. + Y2S);
 
  betaytau1 = ytau * (-(9*g2)/4 - (15*gp2)/4 + (3.*ytau2)/2. + Y2S);
 
  betak1 = (9*g4)/8 + (3*g2*gp2)/4 + (3*gp4)/8 - 9*g2*k - 3*gp2*k +
    24*k2 - 6*yb4 - 6*yt4 + 4*k*Y2S - 2*ytau4;
 
  betam21 = (-9*g2)/2 - (3*gp2)/2 + 12*k + 2*Y2S;
 
  betaLambda1 = 2;

  gamma1 = (-9*g2)/4 - (3*gp2)/4 + Y2S;

  /* One-loop beta functions for light quark/lepton Yukawa added 2/12/2019 */
  betayc1 = yc * (-8*g32 - (9*g2)/4 - (17*gp2)/12 + 3*(yc*yc - ys*ys)/2 + Y2S);
  betays1 = ys * (-8*g32 - (9*g2)/4 - (5*gp2)/12 + 3*(ys*ys - yc*yc)/2 + Y2S); 
  betaymu1 = ymu * (-(9*g2)/4 - (15*gp2)/4 + (3*ymu*ymu)/2 + Y2S); 
  betayu1 = yu * (-8*g32 - (9*g2)/4 - (17*gp2)/12 + Y2S);
  betayd1 = yd * (-8*g32 - (9*g2)/4 - (5*gp2)/12 + Y2S); 
  betaye1 = ye * (-(9*g2)/4 - (15*gp2)/4 + Y2S); 

  beta_g3 = ONELOOPFACTOR * betag31;
  beta_g  = ONELOOPFACTOR * betag1;
  beta_gp = ONELOOPFACTOR * betagp1;
  beta_yt = ONELOOPFACTOR * betayt1;
  beta_yb = ONELOOPFACTOR * betayb1;
  beta_ytau = ONELOOPFACTOR * betaytau1; 
  beta_k = ONELOOPFACTOR * betak1;
  beta_m2 = ONELOOPFACTOR * betam21;
  beta_Lambda = ONELOOPFACTOR * betaLambda1;
  gamma = ONELOOPFACTOR * gamma1;
  beta_yc = ONELOOPFACTOR * betayc1;
  beta_yu = ONELOOPFACTOR * betayu1;
  beta_ys = ONELOOPFACTOR * betays1;
  beta_yd = ONELOOPFACTOR * betayd1;
  beta_ye = ONELOOPFACTOR * betaye1;
  beta_ymu = ONELOOPFACTOR * betaymu1;

  if (loopOrder > 1) {

    betag32 = g3 * g32 * ((9*g2)/2 - 26*g32 + (11*gp2)/6 - 2*yb2 - 2*yt2);

    betag2 = g * g2 * ((35*g2)/6 + 12*g32 + (3*gp2)/2 - (3*yb2)/2 -
      (3*yt2)/2 - (ytau2)/2);

    betagp2 = gp * gp2 * ((9*g2)/2 + (44*g32)/3 + (199*gp2)/18 -
      (5*yb2)/6 - (17*yt2)/6 - (5*ytau2)/2);

    betayt2 = yt * ((-23*g4)/4 + 9*g2*g32 - 108*g34 - (3*g2*gp2)/4 +
      (19*g32*gp2)/9 + (1187*gp4)/216 + 6*k2 +
      (99*g2*yb2)/16 + 4*g32*yb2 + (7*gp2*yb2)/48 -
      yb4/4 + (225*g2*yt2)/16 + 36*g32*yt2 + (131*gp2*yt2)/16 -
      12*k*yt2 - (11*yb2*yt2)/4 - 12*yt4 + (15*g2*ytau2)/8 +
      (25*gp2*ytau2)/8 + (5*yb2*ytau2)/4 - (9*yt2*ytau2)/4 -
      (9*ytau4)/4);

    betayb2 = yb * ((-23*g4)/4 + 9*g2*g32 - 108*g34 - (9*g2*gp2)/4 + 
       (31*g32*gp2)/9 - (127*gp4)/216 + 6*k2 + (225*g2*yb2)/16 + 
       36*g32*yb2 + (79*gp2*yb2)/16 - 12*k*yb2 - 12*yb4 + 
       (99*g2*yt2)/16 + 4*g32*yt2 + (91*gp2*yt2)/48 - 
       (11*yb2*yt2)/4 - yt4/4 + (15*g2*ytau2)/8 + 
       (25*gp2*ytau2)/8 - (9*yb2*ytau2)/4 + (5*yt2*ytau2)/4 - 
       (9*ytau4)/4); 

    betaytau2 = ytau * ((-23*g4)/4 + (9*g2*gp2)/4 + (457*gp4)/24 + 
       6*k2 + (45*g2*yb2)/8 + 20*g32*yb2 + 
       (25*gp2*yb2)/24 - (27*yb4)/4 + (45*g2*yt2)/8 + 
       20*g32*yt2 + (85*gp2*yt2)/24 + (3*yb2*yt2)/2 - 
       (27*yt4)/4 + (165*g2*ytau2)/16 + (179*gp2*ytau2)/16 - 
       12*k*ytau2 - (27*yb2*ytau2)/4 - (27*yt2*ytau2)/4 - 3*ytau4); 

    /* 2-loop beta functions for light quark/lepton Yukawa added 2/12/2019 */

    betayc2 = yc * ((-23*g4)/4 + 9*g2*g32 - 108*g34 - (3*g2*gp2)/4 +
      (19*g32*gp2)/9 + (1187*gp4)/216 + 6*k2 - chi4S + (5*Y4S)/2);

    betayu2 = yu * ((-23*g4)/4 + 9*g2*g32 - 108*g34 - (3*g2*gp2)/4 +
      (19*g32*gp2)/9 + (1187*gp4)/216 + 6*k2 - chi4S + (5*Y4S)/2);

    betays2 = ys * ((-23*g4)/4 + 9*g2*g32 - 108*g34 - (9*g2*gp2)/4 + 
       (31*g32*gp2)/9 - (127*gp4)/216 + 6*k2 - chi4S + (5*Y4S)/2);

    betayd2 = yd * ((-23*g4)/4 + 9*g2*g32 - 108*g34 - (9*g2*gp2)/4 + 
       (31*g32*gp2)/9 - (127*gp4)/216 + 6*k2 - chi4S + (5*Y4S)/2);

    betaymu2 = ymu * ((-23*g4)/4 + (9*g2*gp2)/4 + (457*gp4)/24 + 
       6*k2 - chi4S + (5*Y4S)/2);

    betaye2 = ye * ((-23*g4)/4 + (9*g2*gp2)/4 + (457*gp4)/24 + 
       6*k2 - chi4S + (5*Y4S)/2);

    betak2 = (305*g6)/16 - (289*g4*gp2)/48 - (559*g2*gp4)/48 -
      (379*gp6)/48 - (73*g4*k)/8 + (39*g2*gp2*k)/4 + (629*gp4*k)/24 +
      108*g2*k2 + 36*gp2*k2 - 312*k3 - (9*g4*yb2)/4 +
      (9*g2*gp2*yb2)/2 + (5*gp4*yb2)/4 + (45*g2*k*yb2)/2 +
      80*g32*k*yb2 + (25*gp2*k*yb2)/6 - 144*k2*yb2 - 32*g32*yb4 +
      (4*gp2*yb4)/3 - 3*k*yb4 + 30*yb6 - (9*g4*yt2)/4 +
      (21*g2*gp2*yt2)/2 - (19*gp4*yt2)/4 + (45*g2*k*yt2)/2 +
      80*g32*k*yt2 + (85*gp2*k*yt2)/6 - 144*k2*yt2 - 42*k*yb2*yt2 -
      6*yb4*yt2 - 32*g32*yt4 - (8*gp2*yt4)/3 - 3*k*yt4 - 6*yb2*yt4 +
      30*yt6 - (3*g4*ytau2)/4 + (11*g2*gp2*ytau2)/2 -
      (25*gp4*ytau2)/4 + (15*g2*k*ytau2)/2 + (25*gp2*k*ytau2)/2 -
      48*k2*ytau2 - 4*gp2*ytau4 - k*ytau4 + 10*ytau6;
 
    betam22 = (-145*g4)/16 + (15*g2*gp2)/8 + (557*gp4)/48 + 72*g2*k +
      24*gp2*k - 60*k2 + (45*g2*yb2)/4 + 40*g32*yb2 +
      (25*gp2*yb2)/12 - 72*k*yb2 - (27*yb4)/2 + (45*g2*yt2)/4 +
      40*g32*yt2 + (85*gp2*yt2)/12 - 72*k*yt2 - 21*yb2*yt2 -
      (27*yt4)/2 + (15*g2*ytau2)/4 + (25*gp2*ytau2)/4 - 24*k*ytau2 -
      (9*ytau4)/2;
 
    betaLambda2 = 12 * g2 + 4 * gp2 - 4 * Y2S;

    gamma2 = (-271*g4)/32 + (9*g2*gp2)/16 + (431*gp4)/96 + 6*k2 +
      (45*g2*yb2)/8 + 20*g32*yb2 + (25*gp2*yb2)/24 - (27*yb4)/4 +
      (45*g2*yt2)/8 + 20*g32*yt2 + (85*gp2*yt2)/24 + (3*yb2*yt2)/2 -
      (27*yt4)/4 + (15*g2*ytau2)/8 + (25*gp2*ytau2)/8 - (9*ytau4)/4;

    beta_g3 += TWOLOOPFACTOR * betag32;
    beta_g += TWOLOOPFACTOR * betag2;
    beta_gp += TWOLOOPFACTOR * betagp2;
    beta_yt += TWOLOOPFACTOR * betayt2;
    beta_yb += TWOLOOPFACTOR * betayb2; 
    beta_ytau += TWOLOOPFACTOR * betaytau2; 
    beta_yc += TWOLOOPFACTOR * betayc2; 
    beta_yu += TWOLOOPFACTOR * betayu2; 
    beta_ys += TWOLOOPFACTOR * betays2; 
    beta_yd += TWOLOOPFACTOR * betayd2; 
    beta_ymu += TWOLOOPFACTOR * betaymu2; 
    beta_ye += TWOLOOPFACTOR * betaye2; 
    beta_k += TWOLOOPFACTOR * betak2;
    beta_m2 += TWOLOOPFACTOR * betam22;
    beta_Lambda += TWOLOOPFACTOR * betaLambda2;
    gamma += TWOLOOPFACTOR * gamma2;
  }

  if (loopOrder > 2) {

    betag33 = g3 * g32 * ((109*g4)/8 + 21*g2*g32 + (65*g34)/2 - (g2*gp2)/8 +
      (77*g32*gp2)/9 - (2615*gp4)/216 - (93*g2*yb2)/8 -
      40*g32*yb2 - (89*gp2*yb2)/24 + 15*yb4 -
      (93*g2*yt2)/8 - 40*g32*yt2 - (101*gp2*yt2)/24 +
      18*yb2*yt2 + 15*yt4 + (7*yb2*ytau2)/2 +
      (7*yt2*ytau2)/2);

    betag3 = g * g2 * ((324953*g4)/1728 + 39*g2*g32 + 81*g34 + (291*g2*gp2)/32 -
      (g32*gp2)/3 - (5597*gp4)/576 + (3*g2*k)/2 + (gp2*k)/2 -
      3*k2 - (729*g2*yb2)/32 - 7*g32*yb2 -
      (533*gp2*yb2)/96 + (147*yb4)/16 - (729*g2*yt2)/32 -
      7*g32*yt2 - (593*gp2*yt2)/96 + (117*yb2*yt2)/8 +
      (147*yt4)/16 - (243*g2*ytau2)/32 - (85*gp2*ytau2)/32 +
      (15*yb2*ytau2)/4 + (15*yt2*ytau2)/4 + (29*ytau4)/16);

    betagp3 = gp * gp2 * ((1315*g4)/64 - g2*g32 + 99*g34 +
      (205*g2*gp2)/96 - (137*g32*gp2)/27 - (388613*gp4)/5184 +
      (3*g2*k)/2 + (3*gp2*k)/2 - 3*k2 - (437*g2*yb2)/32 -
      (17*g32*yb2)/3 - (1267*gp2*yb2)/288 + (95*yb4)/16 -
      (785*g2*yt2)/32 - (29*g32*yt2)/3 - (2827*gp2*yt2)/288 +
      (123*yb2*yt2)/8 + (315*yt4)/16 - (543*g2*ytau2)/32 -
      (281*gp2*ytau2)/32 + (157*yb2*ytau2)/12 +
      (199*yt2*ytau2)/12 + (153*ytau4)/16);

    betayt3 = yt * (
 (987*g2*g34)/4 - (4166*g36)/3 + (435*g32*g4)/4 + (455*g6)/576 - 
 (107*g2*g32*gp2)/4 + (1633*g34*gp2)/36 + (273*g4*gp2)/64 + 
 (2045*g2*gp4)/192 + (2047*g32*gp4)/54 + (763523*gp6)/5184 - (171*g4*k)/16 + 
 (39*g2*gp2*k)/8 - (121*gp4*k)/16 + 45*g2*k2 + 15*gp2*k2 - 36*k3 - 
 (27*g2*g32*yb2)/2 - (305*g34*yb2)/2 + (10341*g4*yb2)/256 + 
 (1245*g2*gp2*yb2)/128 - (457*g32*gp2*yb2)/18 - (40673*gp4*yb2)/6912 - 
 (291*k2*yb2)/4 - (2283*g2*yb4)/32 + 82*g32*yb4 - (959*gp2*yb4)/96 + 
 15*k*yb4 + (477*yb6)/16 - 168*g2*g32*yt2 + (3827*g34*yt2)/6 + 
 (32391*g4*yt2)/256 + (2699*g2*gp2*yt2)/128 - 42*g32*gp2*yt2 - 
 (458179*gp4*yt2)/6912 - (135*g2*k*yt2)/2 + 16*g32*k*yt2 - 
 (127*gp2*k*yt2)/6 + (15*k2*yt2)/4 - (2307*g2*yb2*yt2)/32 + 27*g32*yb2*yt2 - 
 (461*gp2*yb2*yt2)/32 + 93*k*yb2*yt2 + (825*yb4*yt2)/8 - (1593*g2*yt4)/16 - 
 157*g32*yt4 - (2437*gp2*yt4)/48 + 198*k*yt4 + (739*yb2*yt4)/16 + 
 (339*yt6)/8 + (1113*g4*ytau2)/128 - (347*g2*gp2*ytau2)/64 - 
 (48295*gp4*ytau2)/1152 - (45*k2*ytau2)/2 - (153*g2*yb2*ytau2)/8 - 
 (43*g32*yb2*ytau2)/6 + (491*gp2*yb2*ytau2)/72 + 22*yb4*ytau2 - 
 (81*g2*yt2*ytau2)/4 + (5*g32*yt2*ytau2)/2 - 21*gp2*yt2*ytau2 + 
 30*k*yt2*ytau2 + (7*yb2*yt2*ytau2)/2 + (21*yt4*ytau2)/2 - 
 (315*g2*ytau4)/16 - (45*gp2*ytau4)/16 + 15*k*ytau4 + (53*yb2*ytau4)/4 + 
 (207*yt2*ytau4)/8 + (71*ytau6)/16 - 144*g2*g34*Zeta3 + 640*g36*Zeta3 - 
 108*g32*g4*Zeta3 + (1125*g6*Zeta3)/8 - (176*g34*gp2*Zeta3)/3 - 
 (81*g4*gp2*Zeta3)/8 - (153*g2*gp4*Zeta3)/8 - (748*g32*gp4*Zeta3)/9 - 
 (13073*gp6*Zeta3)/216 - 108*g2*g32*yb2*Zeta3 - 44*g34*yb2*Zeta3 - 
 (225*g4*yb2*Zeta3)/8 + (9*g2*gp2*yb2*Zeta3)/2 - (28*g32*gp2*yb2*Zeta3)/3 - 
 (199*gp4*yb2*Zeta3)/72 + (63*g2*yb4*Zeta3)/2 - 64*g32*yb4*Zeta3 + 
 (19*gp2*yb4*Zeta3)/6 + (9*yb6*Zeta3)/2 + 180*g2*g32*yt2*Zeta3 - 
 228*g34*yt2*Zeta3 - (729*g4*yt2*Zeta3)/8 + (123*g2*gp2*yt2*Zeta3)/4 + 
 60*g32*gp2*yt2*Zeta3 - (31*gp4*yt2*Zeta3)/24 - (9*g2*yb2*yt2*Zeta3)/2 - 
 32*g32*yb2*yt2*Zeta3 + (5*gp2*yb2*yt2*Zeta3)/6 - 48*yb4*yt2*Zeta3 + 
 (27*yt6*Zeta3)/2 - (81*g4*ytau2*Zeta3)/4 - 3*g2*gp2*ytau2*Zeta3 - 
 (269*gp4*ytau2*Zeta3)/12 + 9*g2*yb2*ytau2*Zeta3 - 9*gp2*yb2*ytau2*Zeta3 - 
 9*g2*yt2*ytau2*Zeta3 + 8*gp2*yt2*ytau2*Zeta3 + 9*g2*ytau4*Zeta3 - 
 9*gp2*ytau4*Zeta3 + 3*ytau6*Zeta3);

  betayb3 = yb * (
(987*g2*g34)/4 - (4166*g36)/3 + (435*g32*g4)/4 + (455*g6)/576 - 
 (51*g2*g32*gp2)/4 + (4165*g34*gp2)/36 - (211*g4*gp2)/64 + 
 (3565*g2*gp4)/192 - (337*g32*gp4)/27 + (93241*gp6)/1728 - (171*g4*k)/16 - 
 (9*g2*gp2*k)/8 - (25*gp4*k)/16 + 45*g2*k2 + 15*gp2*k2 - 36*k3 - 
 168*g2*g32*yb2 + (3827*g34*yb2)/6 + (32391*g4*yb2)/256 + 
 (2831*g2*gp2*yb2)/128 - 30*g32*gp2*yb2 - (209659*gp4*yb2)/6912 - 
 (135*g2*k*yb2)/2 + 16*g32*k*yb2 - (139*gp2*k*yb2)/6 + (15*k2*yb2)/4 - 
 (1593*g2*yb4)/16 - 157*g32*yb4 - (1981*gp2*yb4)/48 + 198*k*yb4 + 
 (339*yb6)/8 - (27*g2*g32*yt2)/2 - (305*g34*yt2)/2 + (10341*g4*yt2)/256 + 
 (1089*g2*gp2*yt2)/128 - (805*g32*gp2*yt2)/18 - (104729*gp4*yt2)/6912 - 
 (291*k2*yt2)/4 - (2307*g2*yb2*yt2)/32 + 27*g32*yb2*yt2 - 
 (1401*gp2*yb2*yt2)/32 + 93*k*yb2*yt2 + (739*yb4*yt2)/16 - (2283*g2*yt4)/32 + 
 82*g32*yt4 - (121*gp2*yt4)/32 + 15*k*yt4 + (825*yb2*yt4)/8 + (477*yt6)/16 + 
 (3033*g4*ytau2)/128 - (411*g2*gp2*ytau2)/64 - (45239*gp4*ytau2)/1152 - 
 (45*k2*ytau2)/2 - (81*g2*yb2*ytau2)/4 + (5*g32*yb2*ytau2)/2 - 
 (137*gp2*yb2*ytau2)/6 + 30*k*yb2*ytau2 + (21*yb4*ytau2)/2 - 
 (153*g2*yt2*ytau2)/8 - (43*g32*yt2*ytau2)/6 + (293*gp2*yt2*ytau2)/72 + 
 (7*yb2*yt2*ytau2)/2 + 22*yt4*ytau2 - (315*g2*ytau4)/16 - (45*gp2*ytau4)/16 + 
 15*k*ytau4 + (207*yb2*ytau4)/8 + (53*yt2*ytau4)/4 + (71*ytau6)/16 - 
 144*g2*g34*Zeta3 + 640*g36*Zeta3 - 108*g32*g4*Zeta3 + (1125*g6*Zeta3)/8 - 
 (176*g34*gp2*Zeta3)/3 - (81*g4*gp2*Zeta3)/8 - (45*g2*gp4*Zeta3)/8 - 
 (220*g32*gp4*Zeta3)/9 - (3845*gp6*Zeta3)/216 + 180*g2*g32*yb2*Zeta3 - 
 228*g34*yb2*Zeta3 - (729*g4*yb2*Zeta3)/8 - 3*g2*gp2*yb2*Zeta3 + 
 44*g32*gp2*yb2*Zeta3 - (19*gp4*yb2*Zeta3)/8 + (27*yb6*Zeta3)/2 - 
 108*g2*g32*yt2*Zeta3 - 44*g34*yt2*Zeta3 - (225*g4*yt2*Zeta3)/8 + 
 (63*g2*gp2*yt2*Zeta3)/4 + (20*g32*gp2*yt2*Zeta3)/3 - (13*gp4*yt2*Zeta3)/72 - 
 (9*g2*yb2*yt2*Zeta3)/2 - 32*g32*yb2*yt2*Zeta3 + (77*gp2*yb2*yt2*Zeta3)/6 + 
 (63*g2*yt4*Zeta3)/2 - 64*g32*yt4*Zeta3 - (17*gp2*yt4*Zeta3)/6 - 
 48*yb2*yt4*Zeta3 + (9*yt6*Zeta3)/2 - (99*g4*ytau2*Zeta3)/4 + 
 18*g2*gp2*ytau2*Zeta3 + (39*gp4*ytau2*Zeta3)/4 - 9*g2*yb2*ytau2*Zeta3 + 
 7*gp2*yb2*ytau2*Zeta3 + 9*g2*yt2*ytau2*Zeta3 - 6*gp2*yt2*ytau2*Zeta3 + 
 9*g2*ytau4*Zeta3 - 9*gp2*ytau4*Zeta3 + 3*ytau6*Zeta3);

betayc3 = yc * g36 * (640 * Zeta3 - 4166./3.);
betays3 = ys * g36 * (640 * Zeta3 - 4166./3.);
betayu3 = yu * g36 * (640 * Zeta3 - 4166./3.);
betayd3 = yd * g36 * (640 * Zeta3 - 4166./3.);

betaytau3 = ytau * (
 (351*g32*g4)/4 + (455*g6)/576 + (981*g4*gp2)/64 + (771*g2*gp4)/64 + 
 (803*g32*gp4)/4 + (607261*gp6)/1728 - (171*g4*k)/16 + (87*g2*gp2*k)/8 - 
 (345*gp4*k)/16 + 45*g2*k2 + 15*gp2*k2 - 36*k3 - (489*g2*g32*yb2)/4 + 
 (622*g34*yb2)/3 + (9099*g4*yb2)/128 + (1823*g2*gp2*yb2)/64 - 
 (991*g32*gp2*yb2)/36 - (54511*gp4*yb2)/3456 - (135*k2*yb2)/2 - 
 (1161*g2*yb4)/16 + (15*g32*yb4)/2 - (411*gp2*yb4)/16 + 45*k*yb4 + 
 (789*yb6)/16 - (489*g2*g32*yt2)/4 + (622*g34*yt2)/3 + (3339*g4*yt2)/128 + 
 (1139*g2*gp2*yt2)/64 - (2419*g32*gp2*yt2)/36 - (167815*gp4*yt2)/3456 - 
 (135*k2*yt2)/2 - (387*g2*yb2*yt2)/8 + 57*g32*yb2*yt2 - 
 (139*gp2*yb2*yt2)/8 + (831*yb4*yt2)/16 - (1161*g2*yt4)/16 + (15*g32*yt4)/2 - 
 (319*gp2*yt4)/16 + 45*k*yt4 + (831*yb2*yt4)/16 + (789*yt6)/16 + 
 (20259*g4*ytau2)/256 - (909*g2*gp2*ytau2)/128 - (18805*gp4*ytau2)/256 - 
 (135*g2*k*ytau2)/2 - (33*gp2*k*ytau2)/2 + (195*k2*ytau2)/4 - 
 (135*g2*yb2*ytau2)/4 - 96*g32*yb2*ytau2 - (29*gp2*yb2*ytau2)/4 + 
 90*k*yb2*ytau2 + (279*yb4*ytau2)/8 - (135*g2*yt2*ytau2)/4 - 
 96*g32*yt2*ytau2 - (77*gp2*yt2*ytau2)/4 + 90*k*yt2*ytau2 - 
 (87*yb2*yt2*ytau2)/4 + (279*yt4*ytau2)/8 - (531*g2*ytau4)/16 - 
 (765*gp2*ytau4)/16 + 108*k*ytau4 + 9*yb2*ytau4 + 9*yt2*ytau4 - 10*ytau6 - 
 108*g32*g4*Zeta3 + (1125*g6*Zeta3)/8 - (81*g4*gp2*Zeta3)/8 - 
 (405*g2*gp4*Zeta3)/8 - 220*g32*gp4*Zeta3 - (3845*gp6*Zeta3)/24 + 
 108*g2*g32*yb2*Zeta3 - 24*g34*yb2*Zeta3 - (297*g4*yb2*Zeta3)/4 - 
 9*g2*gp2*yb2*Zeta3 + 20*g32*gp2*yb2*Zeta3 - (29*gp4*yb2*Zeta3)/12 + 
 27*g2*yb4*Zeta3 - 72*g32*yb4*Zeta3 + 9*gp2*yb4*Zeta3 + 9*yb6*Zeta3 + 
 108*g2*g32*yt2*Zeta3 - 24*g34*yt2*Zeta3 - (243*g4*yt2*Zeta3)/4 - 
 (99*g2*gp2*yt2*Zeta3)/2 + 68*g32*gp2*yt2*Zeta3 - (1157*gp4*yt2*Zeta3)/12 - 
 48*g32*yb2*yt2*Zeta3 + 8*gp2*yb2*yt2*Zeta3 + 27*g2*yt4*Zeta3 - 
 72*g32*yt4*Zeta3 - 3*gp2*yt4*Zeta3 + 9*yt6*Zeta3 - (333*g4*ytau2*Zeta3)/8 + 
 (117*g2*gp2*ytau2*Zeta3)/2 - (75*gp4*ytau2*Zeta3)/8 - 
 27*g2*yb2*ytau2*Zeta3 + 72*g32*yb2*ytau2*Zeta3 - 3*gp2*yb2*ytau2*Zeta3 - 
 27*g2*yt2*ytau2*Zeta3 + 72*g32*yt2*ytau2*Zeta3 + 6*gp2*yt2*ytau2*Zeta3 + 
 (15*ytau6*Zeta3)/2
);

    betak3 = (228259*g8)/1536 - (459*g6*g32)/4 - (165665*g6*gp2)/1728 -
      (153*g4*g32*gp2)/4 - (81509*g4*gp4)/1728 - (187*g2*g32*gp4)/4 -
      (237787*g2*gp6)/3456 - (187*g32*gp6)/4 - (51845*gp8)/512 +
      (58031*g6*k)/144 + 405*g4*g32*k + (6137*g4*gp2*k)/16 +
      (1549*g2*gp4*k)/4 + 165*g32*gp4*k + (88639*gp6*k)/216 -
      (1389*g4*k2)/4 - 666*g2*gp2*k2 - 836*gp4*k2 - 948*g2*k3 -
      316*gp2*k3 + 7176*k4 - (6849*g6*yb2)/128 + (651*g4*g32*yb2)/4 +
      (6099*g4*gp2*yb2)/128 + (233*g2*g32*gp2*yb2)/2 +
      (7603*g2*gp4*yb2)/128 + (683*g32*gp4*yb2)/12 +
      (48523*gp6*yb2)/1152 - (6957*g4*k*yb2)/32 - 489*g2*g32*k*yb2 +
      (2488*g34*k*yb2)/3 - (3009*g2*gp2*k*yb2)/16 -
      (991*g32*gp2*k*yb2)/9 - (149623*gp4*k*yb2)/864 +
      (639*g2*k2*yb2)/2 - 2448*g32*k2*yb2 + (417*gp2*k2*yb2)/2 +
      1746*k3*yb2 + (9909*g4*yb4)/64 - 31*g2*g32*yb4 -
      (532*g34*yb4)/3 - (3239*g2*gp2*yb4)/96 - (641*g32*gp2*yb4)/9 -
      (104383*gp4*yb4)/1728 - (4977*g2*k*yb4)/4 + 1790*g32*k*yb4 -
      (5737*gp2*k*yb4)/12 + 1719*k2*yb4 + (3411*g2*yb6)/16 -
      76*g32*yb6 + (5111*gp2*yb6)/48 + (117*k*yb6)/4 - (1599*yb8)/4 -
      (6849*g6*yt2)/128 + (651*g4*g32*yt2)/4 + (3487*g4*gp2*yt2)/128 +
      (249*g2*g32*gp2*yt2)/2 + (25441*g2*gp4*yt2)/384 +
      (587*g32*gp4*yt2)/12 + (125503*gp6*yt2)/1152 -
      (6957*g4*k*yt2)/32 - 489*g2*g32*k*yt2 + (2488*g34*k*yt2)/3 -
      (6509*g2*gp2*k*yt2)/16 - (2419*g32*gp2*k*yt2)/9 -
      (203887*gp4*k*yt2)/864 + (639*g2*k2*yt2)/2 - 2448*g32*k2*yt2 -
      (195*gp2*k2*yt2)/2 + 1746*k3*yt2 - (2655*g4*yb2*yt2)/32 -
      16*g2*g32*yb2*yt2 + 384*g34*yb2*yt2 + (1001*g2*gp2*yb2*yt2)/
      48 - (709*gp4*yb2*yt2)/32 - (531*g2*k*yb2*yt2)/2 +
      164*g32*k*yb2*yt2 - (929*gp2*k*yb2*yt2)/6 + 234*k2*yb2*yt2 +
      (477*g2*yb4*yt2)/16 - 4*g32*yb4*yt2 - (2299*gp2*yb4*yt2)/48 +
      (6399*k*yb4*yt2)/4 - (717*yb6*yt2)/4 + (9909*g4*yt4)/64 -
      31*g2*g32*yt4 - (532*g34*yt4)/3 - (1079*g2*gp2*yt4)/96 +
      (931*g32*gp2*yt4)/9 + (67793*gp4*yt4)/1728 - (4977*g2*k*yt4)/4 +
      1790*g32*k*yt4 - (2485*gp2*k*yt4)/12 + 1719*k2*yt4 +
      (477*g2*yb2*yt4)/16 - 4*g32*yb2*yt4 + (1337*gp2*yb2*yt4)/48 +
      (6399*k*yb2*yt4)/4 + (3411*g2*yt6)/16 - 76*g32*yt6 +
      (3467*gp2*yt6)/48 + (117*k*yt6)/4 - (717*yb2*yt6)/4 -
      (1599*yt8)/4 - (2283*g6*ytau2)/128 + (1449*g4*gp2*ytau2)/128 +
      (6017*g2*gp4*ytau2)/128 + (10969*gp6*ytau2)/128 -
      (2319*g4*k*ytau2)/32 - (3771*g2*gp2*k*ytau2)/16 -
      (4903*gp4*k*ytau2)/32 + (213*g2*k2*ytau2)/2 -
      (541*gp2*k2*ytau2)/2 + 582*k3*ytau2 + (9*g4*yb2*ytau2)/4 -
      (5*g2*gp2*yb2*ytau2)/2 + (41*gp4*yb2*ytau2)/12 -
      54*g2*k*yb2*ytau2 - 18*gp2*k*yb2*ytau2 - 432*k2*yb2*ytau2 +
      480*k*yb4*ytau2 - (297*yb6*ytau2)/4 + (9*g4*yt2*ytau2)/4 +
      (29*g2*gp2*yt2*ytau2)/2 + (701*gp4*yt2*ytau2)/12 -
      54*g2*k*yt2*ytau2 - 18*gp2*k*yt2*ytau2 - 432*k2*yt2*ytau2 +
      42*k*yb2*yt2*ytau2 + (45*yb4*yt2*ytau2)/4 + 480*k*yt4*ytau2 +
      (45*yb2*yt4*ytau2)/4 - (297*yt6*ytau2)/4 + (3255*g4*ytau4)/64 -
      (15*g2*gp2*ytau4)/32 + (7777*gp4*ytau4)/64 -
      (1587*g2*k*ytau4)/4 + (507*gp2*k*ytau4)/4 + 717*k2*ytau4 +
      480*k*yb2*ytau4 - 144*yb4*ytau4 + 480*k*yt2*ytau4 +
      24*yb2*yt2*ytau4 - 144*yt4*ytau4 + (1137*g2*ytau6)/16 +
      (135*gp2*ytau6)/16 - (1241*k*ytau6)/4 - (297*yb2*ytau6)/4 -
      (297*yt2*ytau6)/4 - (143*ytau8)/4 - (20061*g8*Zeta3)/64 +
      108*g6*g32*Zeta3 - (405*g6*gp2*Zeta3)/16 + 36*g4*g32*gp2*Zeta3 +
      (2217*g4*gp4*Zeta3)/32 + 44*g2*g32*gp4*Zeta3 +
      (2177*g2*gp6*Zeta3)/48 + 44*g32*gp6*Zeta3 + (12457*gp8*Zeta3)/192 +
      (4419*g6*k*Zeta3)/4 - 432*g4*g32*k*Zeta3 - (393*g4*gp2*k*Zeta3)/4 -
      (147*g2*gp4*k*Zeta3)/4 - 176*g32*gp4*k*Zeta3 -
      (1493*gp6*k*Zeta3)/12 - 1026*g4*k2*Zeta3 - 324*g2*gp2*k2*Zeta3 -
      162*gp4*k2*Zeta3 + 144*g2*k3*Zeta3 + 48*gp2*k3*Zeta3 +
      4032*k4*Zeta3 + (297*g6*yb2*Zeta3)/2 - 108*g4*g32*yb2*Zeta3 +
      18*g4*gp2*yb2*Zeta3 - 72*g2*g32*gp2*yb2*Zeta3 +
      9*g2*gp4*yb2*Zeta3 - 36*g32*gp4*yb2*Zeta3 +
      (5*gp6*yb2*Zeta3)/2 - 351*g4*k*yb2*Zeta3 +
      432*g2*g32*k*yb2*Zeta3 - 96*g34*k*yb2*Zeta3 +
      24*g2*gp2*k*yb2*Zeta3 + 80*g32*gp2*k*yb2*Zeta3 -
      (47*gp4*k*yb2*Zeta3)/3 - 864*g2*k2*yb2*Zeta3 +
      2304*g32*k2*yb2*Zeta3 - 384*gp2*k2*yb2*Zeta3 -
      (819*g4*yb4*Zeta3)/8 + 48*g2*g32*yb4*Zeta3 + 64*g34*yb4*Zeta3 -
      (311*g2*gp2*yb4*Zeta3)/4 + (272*g32*gp2*yb4*Zeta3)/3 -
      (2035*gp4*yb4*Zeta3)/72 + 1026*g2*k*yb4*Zeta3 -
      2592*g32*k*yb4*Zeta3 + 498*gp2*k*yb4*Zeta3 + 1512*k2*yb4*Zeta3 -
      54*g2*yb6*Zeta3 + 480*g32*yb6*Zeta3 - 50*gp2*yb6*Zeta3 -
      396*k*yb6*Zeta3 - 72*yb8*Zeta3 + (297*g6*yt2*Zeta3)/2 -
      108*g4*g32*yt2*Zeta3 + (27*g4*gp2*yt2*Zeta3)/2 -
      72*g2*g32*gp2*yt2*Zeta3 - 6*g2*gp4*yt2*Zeta3 -
      36*g32*gp4*yt2*Zeta3 - 5*gp6*yt2*Zeta3 - 351*g4*k*yt2*Zeta3 +
      432*g2*g32*k*yt2*Zeta3 - 96*g34*k*yt2*Zeta3 +
      354*g2*gp2*k*yt2*Zeta3 + 272*g32*gp2*k*yt2*Zeta3 -
      (449*gp4*k*yt2*Zeta3)/3 - 864*g2*k2*yt2*Zeta3 +
      2304*g32*k2*yt2*Zeta3 - 96*gp2*k2*yt2*Zeta3 +
      117*g4*yb2*yt2*Zeta3 + 192*g2*g32*yb2*yt2*Zeta3 +
      31*g2*gp2*yb2*yt2*Zeta3 - 2*gp4*yb2*yt2*Zeta3 +
      108*g2*k*yb2*yt2*Zeta3 - 192*g32*k*yb2*yt2*Zeta3 -
      4*gp2*k*yb2*yt2*Zeta3 - 1728*k2*yb2*yt2*Zeta3 -
      96*g32*yb4*yt2*Zeta3 + 52*gp2*yb4*yt2*Zeta3 +
      288*k*yb4*yt2*Zeta3 - 72*yb6*yt2*Zeta3 - (819*g4*yt4*Zeta3)/8 +
      48*g2*g32*yt4*Zeta3 + 64*g34*yt4*Zeta3 - (743*g2*gp2*yt4*Zeta3)/
      4 - (112*g32*gp2*yt4*Zeta3)/3 + (2957*gp4*yt4*Zeta3)/72 +
      1026*g2*k*yt4*Zeta3 - 2592*g32*k*yt4*Zeta3 + 114*gp2*k*yt4*Zeta3 +
      1512*k2*yt4*Zeta3 - 96*g32*yb2*yt4*Zeta3 -
      56*gp2*yb2*yt4*Zeta3 + 288*k*yb2*yt4*Zeta3 + 144*yb4*yt4*Zeta3 -
      54*g2*yt6*Zeta3 + 480*g32*yt6*Zeta3 + 34*gp2*yt6*Zeta3 -
      396*k*yt6*Zeta3 - 72*yb2*yt6*Zeta3 - 72*yt8*Zeta3 +
      (99*g6*ytau2*Zeta3)/2 - 3*g4*gp2*ytau2*Zeta3 -
      15*g2*gp4*ytau2*Zeta3 - (15*gp6*ytau2*Zeta3)/2 -
      117*g4*k*ytau2*Zeta3 + 252*g2*gp2*k*ytau2*Zeta3 -
      123*gp4*k*ytau2*Zeta3 - 288*g2*k2*ytau2*Zeta3 +
      192*gp2*k2*ytau2*Zeta3 - (273*g4*ytau4*Zeta3)/8 -
      (381*g2*gp2*ytau4*Zeta3)/4 + (375*gp4*ytau4*Zeta3)/8 +
      342*g2*k*ytau4*Zeta3 - 234*gp2*k*ytau4*Zeta3 +
      504*k2*ytau4*Zeta3 - 18*g2*ytau6*Zeta3 + 66*gp2*ytau6*Zeta3 -
      132*k*ytau6*Zeta3 - 24*ytau8*Zeta3;

    betam23 = (49553*g6)/288 + (405*g4*g32)/2 + (3699*g4*gp2)/32 +
      (1713*g2*gp4)/16 + (165*g32*gp4)/2 + (4357*gp6)/27 +
      (4167*g4*k)/16 - (1701*g2*gp2*k)/8 - (5157*gp4*k)/16 - 126*g2*k2 -
      42*gp2*k2 + 2052*k3 - (3789*g4*yb2)/64 - (489*g2*g32*yb2)/2 +
      (1244*g34*yb2)/3 - (865*g2*gp2*yb2)/32 - (991*g32*gp2*yb2)/18 -
      (101527*gp4*yb2)/1728 + (567*g2*k*yb2)/4 - 1224*g32*k*yb2 +
      (393*gp2*k*yb2)/4 + 297*k2*yb2 - (3177*g2*yb4)/8 + 447*g32*yb4 -
      (1067*gp2*yb4)/8 + (351*k*yb4)/2 + (1605*yb6)/8 -
      (3789*g4*yt2)/64 - (489*g2*g32*yt2)/2 + (1244*g34*yt2)/3 -
      (3277*g2*gp2*yt2)/32 - (2419*g32*gp2*yt2)/18 -
      (214543*gp4*yt2)/1728 + (567*g2*k*yt2)/4 - 1224*g32*k*yt2 -
      (219*gp2*k*yt2)/4 + 297*k2*yt2 - (243*g2*yb2*yt2)/4 +
      82*g32*yb2*yt2 - (929*gp2*yb2*yt2)/12 - 315*k*yb2*yt2 +
      (4047*yb4*yt2)/8 - (3177*g2*yt4)/8 + 447*g32*yt4 -
      (431*gp2*yt4)/8 + (351*k*yt4)/2 + (4047*yb2*yt4)/8 +
      (1605*yt6)/8 - (1263*g4*ytau2)/64 - (2331*g2*gp2*ytau2)/32 -
      (6727*gp4*ytau2)/64 + (189*g2*k*ytau2)/4 - (549*gp2*k*ytau2)/4 +
      99*k2*ytau2 - 27*g2*yb2*ytau2 - 9*gp2*yb2*ytau2 -
      216*k*yb2*ytau2 + 144*yb4*ytau2 - 27*g2*yt2*ytau2 -
      9*gp2*yt2*ytau2 - 216*k*yt2*ytau2 + 21*yb2*yt2*ytau2 +
      144*yt4*ytau2 - (987*g2*ytau4)/8 + (291*gp2*ytau4)/8 +
      (261*k*ytau4)/2 + 144*yb2*ytau4 + 144*yt2*ytau4 - (233*ytau6)/8 +
      (2871*g6*Zeta3)/8 - 216*g4*g32*Zeta3 - (549*g4*gp2*Zeta3)/8 -
      (351*g2*gp4*Zeta3)/8 - 88*g32*gp4*Zeta3 - (1673*gp6*Zeta3)/24 -
      324*g4*k*Zeta3 + 72*g2*gp2*k*Zeta3 - 36*gp4*k*Zeta3 -
      216*g2*k2*Zeta3 - 72*gp2*k2*Zeta3 - (243*g4*yb2*Zeta3)/2 +
      216*g2*g32*yb2*Zeta3 - 48*g34*yb2*Zeta3 + 40*g32*gp2*yb2*Zeta3 -
      (35*gp4*yb2*Zeta3)/6 - 648*g2*k*yb2*Zeta3 + 1152*g32*k*yb2*Zeta3 -
      264*gp2*k*yb2*Zeta3 + 324*g2*yb4*Zeta3 - 720*g32*yb4*Zeta3 +
      144*gp2*yb4*Zeta3 + 432*k*yb4*Zeta3 + 90*yb6*Zeta3 -
      (243*g4*yt2*Zeta3)/2 + 216*g2*g32*yt2*Zeta3 - 48*g34*yt2*Zeta3 +
      117*g2*gp2*yt2*Zeta3 + 136*g32*gp2*yt2*Zeta3 -
      (149*gp4*yt2*Zeta3)/6 - 648*g2*k*yt2*Zeta3 +
      1152*g32*k*yt2*Zeta3 - 120*gp2*k*yt2*Zeta3 -
      54*g2*yb2*yt2*Zeta3 - 96*g32*yb2*yt2*Zeta3 -
      2*gp2*yb2*yt2*Zeta3 - 432*k*yb2*yt2*Zeta3 + 72*yb4*yt2*Zeta3 +
      324*g2*yt4*Zeta3 - 720*g32*yt4*Zeta3 + 24*gp2*yt4*Zeta3 +
      432*k*yt4*Zeta3 + 72*yb2*yt4*Zeta3 + 90*yt6*Zeta3 -
      (81*g4*ytau2*Zeta3)/2 + 90*g2*gp2*ytau2*Zeta3 -
      (15*gp4*ytau2*Zeta3)/2 - 216*g2*k*ytau2*Zeta3 +
      72*gp2*k*ytau2*Zeta3 + 108*g2*ytau4*Zeta3 - 72*gp2*ytau4*Zeta3 +
      144*k*ytau4*Zeta3 + 30*ytau6*Zeta3;

    betaLambda3 = (192 * Zeta3 - 204) * yt2 * g32 + 38.25 * yt4 
      + (189./8. - 108 * Zeta3) * yt2 * g2 - (73./8. + 20 * Zeta3) * yt2 * gp2
      + (1767./32. - 18 * Zeta3) * g4 + (36 * Zeta3 - 441./16.) * g2 * gp2 
      + (6 * Zeta3 - 1593./32.) * gp4;
    
    gamma3 = (6785*g6)/576 + (405*g4*g32)/4 + (543*g4*gp2)/64 +
      (153*g2*gp4)/32 + (165*g32*gp4)/4 + (27053*gp6)/432 +
      (117*g4*k)/16 + (39*g2*gp2*k)/8 + (39*gp4*k)/16 + 45*g2*k2 +
      15*gp2*k2 - 36*k3 + (4275*g4*yb2)/128 - (489*g2*g32*yb2)/4 +
      (622*g34*yb2)/3 + (671*g2*gp2*yb2)/64 - (991*g32*gp2*yb2)/36 -
      (27799*gp4*yb2)/3456 - (135*k2*yb2)/2 - (1161*g2*yb4)/16 +
      (15*g32*yb4)/2 - (411*gp2*yb4)/16 + 45*k*yb4 + (789*yb6)/16 +
      (4275*g4*yt2)/128 - (489*g2*g32*yt2)/4 + (622*g34*yt2)/3 +
      (371*g2*gp2*yt2)/64 - (2419*g32*gp2*yt2)/36 -
      (144271*gp4*yt2)/3456 - (135*k2*yt2)/2 - (387*g2*yb2*yt2)/8 +
      57*g32*yb2*yt2 - (139*gp2*yb2*yt2)/8 + (831*yb4*yt2)/16 -
      (1161*g2*yt4)/16 + (15*g32*yt4)/2 - (319*gp2*yt4)/16 + 45*k*yt4 +
      (831*yb2*yt4)/16 + (789*yt6)/16 + (1425*g4*ytau2)/128 -
      (411*g2*gp2*ytau2)/64 - (5959*gp4*ytau2)/128 - (45*k2*ytau2)/2 -
      (27*g2*yb2*ytau2)/2 - (9*gp2*yb2*ytau2)/2 + 18*yb4*ytau2 -
      (27*g2*yt2*ytau2)/2 - (9*gp2*yt2*ytau2)/2 -
      (3*yb2*yt2*ytau2)/2 + 18*yt4*ytau2 - (315*g2*ytau4)/16 -
      (45*gp2*ytau4)/16 + 15*k*ytau4 + 18*yb2*ytau4 + 18*yt2*ytau4 +
      (71*ytau6)/16 + (1953*g6*Zeta3)/16 - 108*g4*g32*Zeta3 -
      (153*g4*gp2*Zeta3)/16 - (135*g2*gp4*Zeta3)/16 - 44*g32*gp4*Zeta3 -
      (1529*gp6*Zeta3)/48 - (27*g4*k*Zeta3)/2 - 9*g2*gp2*k*Zeta3 -
      (9*gp4*k*Zeta3)/2 - (189*g4*yb2*Zeta3)/4 + 108*g2*g32*yb2*Zeta3 -
      24*g34*yb2*Zeta3 - 9*g2*gp2*yb2*Zeta3 + 20*g32*gp2*yb2*Zeta3 -
      (29*gp4*yb2*Zeta3)/12 + 27*g2*yb4*Zeta3 - 72*g32*yb4*Zeta3 +
      9*gp2*yb4*Zeta3 + 9*yb6*Zeta3 - (189*g4*yt2*Zeta3)/4 +
      108*g2*g32*yt2*Zeta3 - 24*g34*yt2*Zeta3 +
      (27*g2*gp2*yt2*Zeta3)/2 + 68*g32*gp2*yt2*Zeta3 +
      (gp4*yt2*Zeta3)/12 - 48*g32*yb2*yt2*Zeta3 +
      8*gp2*yb2*yt2*Zeta3 + 27*g2*yt4*Zeta3 - 72*g32*yt4*Zeta3 -
      3*gp2*yt4*Zeta3 + 9*yt6*Zeta3 - (63*g4*ytau2*Zeta3)/4 +
      18*g2*gp2*ytau2*Zeta3 + (39*gp4*ytau2*Zeta3)/4 +
      9*g2*ytau4*Zeta3 - 9*gp2*ytau4*Zeta3 + 3*ytau6*Zeta3;

    beta_g3 += THREELOOPFACTOR * betag33;
    beta_g  += THREELOOPFACTOR * betag3;
    beta_gp += THREELOOPFACTOR * betagp3;
    beta_yt += THREELOOPFACTOR * betayt3;
    beta_yb += THREELOOPFACTOR * betayb3; 
    beta_yc += THREELOOPFACTOR * betayc3; 
    beta_ys += THREELOOPFACTOR * betays3; 
    beta_yu += THREELOOPFACTOR * betayu3; 
    beta_yd += THREELOOPFACTOR * betayd3; 
    beta_ytau += THREELOOPFACTOR * betaytau3; 
    beta_k  += THREELOOPFACTOR * betak3;
    beta_m2 += THREELOOPFACTOR * betam23;
    beta_Lambda += THREELOOPFACTOR * betaLambda3;
    gamma   += THREELOOPFACTOR * gamma3;
  }

  /* Known partial results at 4-loop order */
  if (loopOrder > 3) {

    /* See 1402.6611 equation (4.2), each coefficient of a power of as^n 
       should be multipled by 2 4^n */
    betayt4 = yt*g38*2308.18;
    beta_yt += FOURLOOPFACTOR * betayt4;

    /* From 1508.00912 */
    betak4 = 8308.17*yt4*g36;
    beta_k += FOURLOOPFACTOR * betak4;

    /* Old, v1.0. See 1508.03624 and 1508.02680. 
    betag34 = g3*g32*(
      g36 * (63559./18. - (44948.*Zeta3)/9.) +
      g34 * yt2 * (-6709./9. + 272.*Zeta3) +
      g32 * yt4 * (423. - 120.*Zeta3 + Rgammafive*(4./3. + 8.*Zeta3)) + 
      yt6 * (-423./4. - 6.*Zeta3) +
      -30. * k * yt4 + 36. * yt2 * k2);
    */

    /* New, from 1912.07624 */
    betag34 = g3*g32*( 
  (-5969*g2*g34)/12 + (63559*g36)/18 + (953*g32*g4)/9 - (176815*g6)/1728 + 
  (23*g2*g32*gp2)/4 - (57739*g34*gp2)/324 - (37597*g4*gp2)/1728 - 
  (46951*g2*gp4)/1728 - (88855*g32*gp4)/486 - (6085099*gp6)/46656 - 
  (473*g2*g32*yb2)/4 - (6709*g34*yb2)/9 - (12887*g4*yb2)/192 - 
  (775*g2*gp2*yb2)/96 - (1487*g32*gp2*yb2)/36 + (210847*gp4*yb2)/5184 + 
  36*k2*yb2 + (3201*g2*yb4)/32 + 427*g32*yb4 + (2869*gp2*yb4)/96 - 30*k*yb4 - 
  (423*yb6)/4 - (473*g2*g32*yt2)/4 - (6709*g34*yt2)/9 - (12887*g4*yt2)/192 + 
  (77*g2*gp2*yt2)/96 - (1283*g32*gp2*yt2)/36 + (362287*gp4*yt2)/5184 + 
  36*k2*yt2 + (1895*g2*yb2*yt2)/16 + (4282*g32*yb2*yt2)/9 + 
  (19033*gp2*yb2*yt2)/432 - (1171*yb4*yt2)/12 + (3201*g2*yt4)/32 + 427*g32*yt4 + 
  (3641*gp2*yt4)/96 - 30*k*yt4 - (1171*yb2*yt4)/12 - (423*yt6)/4 + 
  (63*g4*ytau2)/8 + (253*gp4*ytau2)/24 + (295*g2*yb2*ytau2)/16 + 
  (133*g32*yb2*ytau2)/3 + (4201*gp2*yb2*ytau2)/144 - (181*yb4*ytau2)/8 + 
  (295*g2*yt2*ytau2)/16 + (133*g32*yt2*ytau2)/3 + (3649*gp2*yt2*ytau2)/144 - 
  (515*yb2*yt2*ytau2)/36 - (181*yt4*ytau2)/8 - (135*yb2*ytau4)/8 - 
  (135*yt2*ytau4)/8 + 869*g2*g34*Zeta3 - (44948*g36*Zeta3)/9 - 
  (475*g32*g4*Zeta3)/6 - (935*g6*Zeta3)/4 + (8119*g34*gp2*Zeta3)/27 + 
  (691*g4*gp2*Zeta3)/36 + (973*g2*gp4*Zeta3)/36 + (11275*g32*gp4*Zeta3)/162 + 
  (87365*gp6*Zeta3)/972 - 72*g2*g32*yb2*Zeta3 + 272*g34*yb2*Zeta3 + 
  (117*g4*yb2*Zeta3)/4 - (15*g2*gp2*yb2*Zeta3)/2 - (104*g32*gp2*yb2*Zeta3)/3 - 
  (7*gp4*yb2*Zeta3)/36 + (45*g2*yb4*Zeta3)/2 - 96*g32*yb4*Zeta3 + 
  (27*gp2*yb4*Zeta3)/2 - 6*yb6*Zeta3 - 72*g2*g32*yt2*Zeta3 + 272*g34*yt2*Zeta3 + 
  (117*g4*yt2*Zeta3)/4 - (45*g2*gp2*yt2*Zeta3)/2 - (8*g32*gp2*yt2*Zeta3)/3 - 
  (19*gp4*yt2*Zeta3)/36 + 3*g2*yb2*yt2*Zeta3 - (400*g32*yb2*yt2*Zeta3)/3 - 
  (25*gp2*yb2*yt2*Zeta3)/9 - 8*yb4*yt2*Zeta3 + (45*g2*yt4*Zeta3)/2 - 
  96*g32*yt4*Zeta3 + (7*gp2*yt4*Zeta3)/2 - 8*yb2*yt4*Zeta3 - 6*yt6*Zeta3 + 
  9*g2*yb2*ytau2*Zeta3 - 9*gp2*yb2*ytau2*Zeta3 + 9*g2*yt2*ytau2*Zeta3 - 
  9*gp2*yt2*ytau2*Zeta3 - (16*yb2*yt2*ytau2*Zeta3)/3
  );

    beta_g3 += FOURLOOPFACTOR * betag34;

    /* New, from 1912.07624 */
    betag4 = g*g2*(
 (2587*g2*g34)/3 + (257*g36)/3 - (72881*g32*g4)/72 + (124660945*g6)/62208 + 
 (161*g2*g32*gp2)/12 - (2185*g34*gp2)/27 - (375767*g4*gp2)/6912 - 
 (787709*g2*gp4)/6912 - (52297*g32*gp4)/648 - (6418229*gp6)/62208 + 
 (2905*g4*k)/48 + (115*g2*gp2*k)/8 + (457*gp4*k)/144 - (363*g2*k2)/4 - 
 (109*gp2*k2)/4 + 52*k3 - (2143*yb6)/32 - (361*g2*g32*yb2)/3 - 
 (307*g34*yb2)/6 - (500665*g4*yb2)/2304 - (90029*g2*gp2*yb2)/1152 + 
 (5*g32*gp2*yb2)/3 + (89515*gp4*yb2)/2304 - (75*g2*k*yb2)/2 - 
 (17*gp2*k*yb2)/2 + 75*k2*yb2 + (30213*g2*yb4)/128 + (239*g32*yb4)/4 + 
 (15937*gp2*yb4)/384 - (39*k*yb4)/2 - (2143*yt6)/32 - (361*g2*g32*yt2)/3 - 
 (307*g34*yt2)/6 - (500665*g4*yt2)/2304 - (102497*g2*gp2*yt2)/1152 + 
 (199*g32*gp2*yt2)/9 + (465089*gp4*yt2)/6912 - (75*g2*k*yt2)/2 - 
 (9*gp2*k*yt2)/2 + 75*k2*yt2 + (16735*g2*yb2*yt2)/64 + (739*g32*yb2*yt2)/6 + 
 (33959*gp2*yb2*yt2)/576 - (3265*yb4*yt2)/32 + (30213*g2*yt4)/128 + 
 (239*g32*yt4)/4 + (15805*gp2*yt4)/384 - (39*k*yt4)/2 - (3265*yb2*yt4)/32 - 
 (269*ytau6)/32 - (500665*g4*ytau2)/6912 - (15341*g2*gp2*ytau2)/384 + 
 (302123*gp4*ytau2)/6912 - (25*g2*k*ytau2)/2 - (gp2*k*ytau2)/6 + 25*k2*ytau2 + 
 (8927*g2*yb2*ytau2)/96 + (33*g32*yb2*ytau2)/2 + (5471*gp2*yb2*ytau2)/288 - 
 (91*yb4*ytau2)/4 + (7007*g2*yt2*ytau2)/96 + (33*g32*yt2*ytau2)/2 + 
 (5401*gp2*yt2*ytau2)/288 - (113*yb2*yt2*ytau2)/6 - (91*yt4*ytau2)/4 + 
 (54931*g2*ytau4)/1152 + (10309*gp2*ytau4)/1152 - (13*k*ytau4)/2 - 
 (41*yb2*ytau4)/2 - (41*yt2*ytau4)/2 - 640*g2*g34*Zeta3 - (1760*g36*Zeta3)/3 + 
 (4108*g32*g4*Zeta3)/3 - (78803*g6*Zeta3)/36 + (736*g34*gp2*Zeta3)/9 + 
 (4631*g4*gp2*Zeta3)/36 + (659*g2*gp4*Zeta3)/36 + (2540*g32*gp4*Zeta3)/27 + 
 (21173*gp6*Zeta3)/324 - (9*yb6*Zeta3)/2 - 14*g2*g32*yb2*Zeta3 + 
 84*g34*yb2*Zeta3 + (239*g4*yb2*Zeta3)/6 - (5*g2*gp2*yb2*Zeta3)/3 - 
 26*g32*gp2*yb2*Zeta3 - (143*gp4*yb2*Zeta3)/18 - (63*g2*yb4*Zeta3)/4 - 
 36*g32*yb4*Zeta3 + (33*gp2*yb4*Zeta3)/4 - (9*yt6*Zeta3)/2 - 
 14*g2*g32*yt2*Zeta3 + 84*g34*yt2*Zeta3 + (239*g4*yt2*Zeta3)/6 - 
 (35*g2*gp2*yt2*Zeta3)/3 - (94*g32*gp2*yt2*Zeta3)/3 - (83*gp4*yt2*Zeta3)/6 + 
 (57*g2*yb2*yt2*Zeta3)/2 - 152*g32*yb2*yt2*Zeta3 + (40*gp2*yb2*yt2*Zeta3)/3 + 
 6*yb4*yt2*Zeta3 - (63*g2*yt4*Zeta3)/4 - 36*g32*yt4*Zeta3 + 
 (51*gp2*yt4*Zeta3)/4 + 6*yb2*yt4*Zeta3 - (3*ytau6*Zeta3)/2 + 
 (239*g4*ytau2*Zeta3)/18 - 8*g2*gp2*ytau2*Zeta3 - (97*gp4*ytau2*Zeta3)/18 - 
 g2*yb2*ytau2*Zeta3 - 16*g32*yb2*ytau2*Zeta3 + (7*gp2*yb2*ytau2*Zeta3)/3 + 
 5*g2*yt2*ytau2*Zeta3 - 16*g32*yt2*ytau2*Zeta3 - 12*gp2*yt2*ytau2*Zeta3 + 
 4*yb2*yt2*ytau2*Zeta3 - (59*g2*ytau4*Zeta3)/12 + (49*gp2*ytau4*Zeta3)/12 
  );

    beta_g += FOURLOOPFACTOR * betag4;

    /* New, from 1912.07624 */
    betagp4 = gp*gp2*(
    (-2185*g2*g34)/9 + (1529*g36)/9 - (41971*g32*g4)/216 - (117923*g6)/6912 - 
 (23*g2*g32*gp2)/12 + (83389*g34*gp2)/243 + (572059*g4*gp2)/20736 - 
 (3819731*g2*gp4)/20736 - (3629273*g32*gp4)/5832 - (143035709*gp6)/559872 + 
 (889*g4*k)/48 + (213*g2*gp2*k)/8 + (403*gp4*k)/48 - (327*g2*k2)/4 - 
 (141*gp2*k2)/4 + 52*k3 - (1595*yb6)/32 - 23*g2*g32*yb2 - (415*g34*yb2)/54 - 
 (80311*g4*yb2)/768 - (4145*g2*gp2*yb2)/128 - (395*g32*gp2*yb2)/27 + 
 (3876629*gp4*yb2)/62208 - (51*g2*k*yb2)/2 - (59*gp2*k*yb2)/2 + 63*k2*yb2 + 
 (14977*g2*yb4)/128 + (497*g32*yb4)/12 + (47975*gp2*yb4)/1152 - (19*k*yb4)/2 - 
 (4551*yt6)/32 + (367*g2*g32*yt2)/3 - (5731*g34*yt2)/54 - 
 (439841*g4*yt2)/2304 - (42841*g2*gp2*yt2)/1152 - (503*g32*gp2*yt2)/27 + 
 (8978897*gp4*yt2)/62208 - (27*g2*k*yt2)/2 - (107*gp2*k*yt2)/2 + 99*k2*yt2 + 
 (30655*g2*yb2*yt2)/192 + (6901*g32*yb2*yt2)/54 + (306401*gp2*yb2*yt2)/5184 - 
 (32269*yb4*yt2)/288 + (23821*g2*yt4)/128 + (1429*g32*yt4)/12 + 
 (145295*gp2*yt4)/1152 - (79*k*yt4)/2 - (29281*yb2*yt4)/288 - 
 (4123*ytau6)/96 - (305839*g4*ytau2)/2304 + (2901*g2*gp2*ytau2)/128 + 
 (1082567*gp4*ytau2)/6912 - (g2*k*ytau2)/2 - (73*gp2*k*ytau2)/2 + 
 61*k2*ytau2 + (2869*g2*yb2*ytau2)/32 + (3613*g32*yb2*ytau2)/18 + 
 (11837*gp2*yb2*ytau2)/864 - 80*yb4*ytau2 + (9713*g2*yt2*ytau2)/96 + 
 (3061*g32*yt2*ytau2)/18 + (70529*gp2*yt2*ytau2)/864 - 
 (5933*yb2*yt2*ytau2)/216 - (2371*yt4*ytau2)/24 + (31141*g2*ytau4)/384 + 
 (89161*gp2*ytau4)/1152 - (73*k*ytau4)/2 - (249*yb2*ytau4)/4 - 
 (633*yt2*ytau4)/8 + (736*g2*g34*Zeta3)/3 - (23200*g36*Zeta3)/27 + 
 (1868*g32*g4*Zeta3)/9 - (3109*g6*Zeta3)/12 - (68656*g34*gp2*Zeta3)/81 - 
 (6751*g4*gp2*Zeta3)/108 + (16529*g2*gp4*Zeta3)/108 + 
 (180076*g32*gp4*Zeta3)/243 - (1638851*gp6*Zeta3)/2916 - (5*yb6*Zeta3)/2 - 
 46*g2*g32*yb2*Zeta3 - (20*g34*yb2*Zeta3)/3 + (62*g4*yb2*Zeta3)/3 - 
 (5*g2*gp2*yb2*Zeta3)/6 - 10*g32*gp2*yb2*Zeta3 - (83*gp4*yb2*Zeta3)/18 + 
 (33*g2*yb4*Zeta3)/4 - 4*g32*yb4*Zeta3 + (5*gp2*yb4*Zeta3)/4 - 
 (17*yt6*Zeta3)/2 - 158*g2*g32*yt2*Zeta3 + (796*g34*yt2*Zeta3)/3 + 
 (154*g4*yt2*Zeta3)/3 - (187*g2*gp2*yt2*Zeta3)/6 + (34*g32*gp2*yt2*Zeta3)/3 + 
 (433*gp4*yt2*Zeta3)/18 + (127*g2*yb2*yt2*Zeta3)/2 - 
 (1208*g32*yb2*yt2*Zeta3)/9 - (125*gp2*yb2*yt2*Zeta3)/27 + 
 (14*yb4*yt2*Zeta3)/3 + (213*g2*yt4*Zeta3)/4 - 100*g32*yt4*Zeta3 - 
 (119*gp2*yt4*Zeta3)/12 + (14*yb2*yt4*Zeta3)/3 - (15*ytau6*Zeta3)/2 + 
 (121*g4*ytau2*Zeta3)/3 - (135*g2*gp2*ytau2*Zeta3)/2 + 
 (989*gp4*ytau2*Zeta3)/18 + 47*g2*yb2*ytau2*Zeta3 - 112*g32*yb2*ytau2*Zeta3 + 
 (25*gp2*yb2*ytau2*Zeta3)/9 + 19*g2*yt2*ytau2*Zeta3 - 112*g32*yt2*ytau2*Zeta3 - 
 (1370*gp2*yt2*ytau2*Zeta3)/9 + (28*yb2*yt2*ytau2*Zeta3)/9 + 
 (55*g2*ytau4*Zeta3)/4 - (155*gp2*ytau4*Zeta3)/12
  );

    beta_gp += FOURLOOPFACTOR * betagp4;
  }

  /* Known partial results at 5-loop order */
  if (loopOrder > 4) {

    /* See 1402.6611 equation (4.2), each coefficient of a power of as^n 
       should be multipled by 2 4^n */
    betayt5 = yt*g32*g38*19988.8;
    beta_yt += ONELOOPFACTOR*FOURLOOPFACTOR * betayt5;

    /* See 1701.01404 and 1606.08659, but beware of roundoff error in 
       eq (6) of latter! */
    betag35 = g3*g34*g36*(-271.428);
    beta_g3 += ONELOOPFACTOR*FOURLOOPFACTOR * betag35;
  }
    
  beta_m2 *= m2;
  beta_Lambda *= (m2 * m2);
  beta_v = -v*gamma;  

  return 0;
}  
