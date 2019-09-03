#include "smdr_internal.h"

/* 
  Sources: Schroder and Steinhauser, 0512058 for 4-loop QCD decoupling
  relations for QCD. (In paper, earlier <4-loop contributions should be cited.)

  Equation (27) of 0004189 for Rundec (Chetyrkin, Kuehn, Steinhauser)
  for  3-loop decoupling of running quark masses. Original source is
  9708255 Chetyrkin, Kniehl Steinhauser, equation (20).

  QED contributions at 2-loop order from 1212.5144 Grozin. 

  The direct source for the code below is from 1812.04100.
*/


/* ------------------------------------------------------------------- */
SMDR_REAL z_mass_corr (SMDR_REAL r, int nq, SMDR_REAL lnbarF,
                       SMDR_REAL alpS, int loopOrder)
{
  SMDR_REAL lnr, Delta2r, Delta3r, result;

  if (r < SMDR_TOL) return 0;
 
  lnr = SMDR_LOG(r);

  /* 1812.04100 eq. (3.12) */
  Delta2r = r*(8.*lnr/15. - 76./75.) + r*r*(9.*lnr/70. - 1389./9800.);

  /* 1812.04100 eq. (3.11) */
  if (loopOrder > 1.999) {
    result = alpS * alpS * Delta2r/(12. * PI * PI);
  }

  /* 1812.04100 eqs. (3.13) and (3.15) */
  if (loopOrder > 2.999) {
    Delta3r = (8./9.)*(2.*nq - 31.) * lnbarF * Delta2r +
              r * ((64. * nq/135. - 451./81.) * lnr * lnr +
                  (84887./7290. - 128.*nq/135.) * lnr + 
                  2.77670 - 0.22452 * nq) +
              r * r * ((4. * nq/35. - 239./270.) * lnr * lnr +
                  (580157./396900. - 6.*nq/35.) * lnr + 
                  0.52092 + 0.03556 * nq);

    result += alpS * alpS * alpS * Delta3r/(64. * PI * PI * PI);
  }

  return(result);
};

/* ------------------------------------------------------------------- */

void SMDR_Decouple_bottom (int loopOrder)
{
  SMDR_REAL zc, zs, zu, zd, ze;
  SMDR_REAL rcharm, rstrange, lnbarb;

  if (mb_53 > 1.0) {
    SMDR_QCDQEDmatch ("d", Q_53, mb_53, alphaS_53, alpha_53, 4, loopOrder,
                      &alphaS_43, &alpha_43, &zu, &zd, &ze); 
    rcharm = (mc_53 * mc_53)/(mb_53 * mb_53);
    rstrange = (ms_53 * ms_53)/(mb_53 * mb_53);
    lnbarb = SMDR_LOG(mb_53 * mb_53/(Q_53 * Q_53));
    zc = zu + z_mass_corr(rcharm, 4, lnbarb, alphaS_53, 4);
    zs = zd + z_mass_corr(rstrange, 4, lnbarb, alphaS_53, 4);
  } else {
    SMDR_QCDQEDmatch ("d", Q_53, Q_53, alphaS_53, alpha_53, 4, loopOrder,
                      &alphaS_43, &alpha_43, &zu, &zd, &ze); 
    zc = zu;
    zs = zd;
  }

  Q_43 = Q_53;
  mc_43 = mc_53 * zc;
  ms_43 = ms_53 * zs;
  md_43 = md_53 * zd;
  mu_43 = mu_53 * zu;
  mtau_43 = mtau_53 * ze;
  mmuon_43 = mmuon_53 * ze;
  melectron_43 = melectron_53 * ze;

  return;
}


/* ------------------------------------------------------------------------ */

void SMDR_Decouple_tau (int loopOrder)
{
  SMDR_REAL zu, zd, ze;

  if (mtau_43 > 1.0) {
    SMDR_QCDQEDmatch ("e", Q_43, mtau_43, alphaS_43, alpha_43, 4, loopOrder,
                      &alphaS_42, &alpha_42, &zu, &zd, &ze);
  } else {
    SMDR_QCDQEDmatch ("e", Q_43, Q_43, alphaS_43, alpha_43, 4, loopOrder,
                      &alphaS_42, &alpha_42, &zu, &zd, &ze);
  }

  Q_42 = Q_43;
  mc_42 = mc_43 * zu; 
  ms_42 = ms_43 * zd;
  md_42 = md_43 * zd;
  mu_42 = mu_43 * zu;
  mmuon_42 = mmuon_43 * ze;
  melectron_42 = melectron_43 * ze;

  return;
}

/* ------------------------------------------------------------------------ */

void SMDR_Decouple_charm (int loopOrder)
{
  SMDR_REAL zs, zu, zd, ze;
  SMDR_REAL rstrange, lnbarc;

  if (mc_42 > 0.3) {
    SMDR_QCDQEDmatch ("u", Q_42, mc_42, alphaS_42, alpha_42, 4, loopOrder,
                      &alphaS_32, &alpha_32, &zu, &zd, &ze);
    rstrange = (ms_42 * ms_42)/(mc_42 * mc_42);
    lnbarc = SMDR_LOG(mc_42 * mc_42/(Q_42 * Q_42));
    zs = zd + z_mass_corr(rstrange, 3, lnbarc, alphaS_42, 4);
  } else {
    SMDR_QCDQEDmatch ("u", Q_42, Q_42, alphaS_42, alpha_42, 4, loopOrder,
                      &alphaS_32, &alpha_32, &zu, &zd, &ze);
    zs = zd;
  }

  Q_32 = Q_42;
  ms_32 = ms_42 * zs;
  md_32 = md_42 * zd;
  mu_32 = mu_42 * zu;
  mmuon_32 = mmuon_42 * ze;
  melectron_32 = melectron_42 * ze;

  return;
}

/* ------------------------------------------------------------------------ 
   The input arguments are: 

     Fermion_type = "u", "d", or "e" = type of heavy fermion being decoupled
     Q_match   = MSbar decoupling scale.
     m_Fermion = MSbar mass of heavy fermion being decoupled, in full theory,
                 at Q_match                     
     alphaS_hi  = strong coupling in full theory, at scale Q_match.
     alpha_hi   = EM coupling in full theory, at scale Q_match
     nqlight = number of light quarks, not including the fermion being 
               decoupled. Note that the number of light leptons doesn't matter
               in the approximation being used. (It would matter if 3-loop
               EM contributions were included.)
     loopOrder = 0, 1, 2, 3, 4. EM contributions are neglected beyond 2 loops.

   The outputs are:

     alphaS_lo = strong coupling in decoupled theory, at scale Q_match.
     alpha_lo  = EM coupling in decoupled theory, at scale Q_match.
     zu = ratio of (MSbar light up-type mass in decoupled theory)/
                   (MSbar light up-type mass in full theory). 
     zd = ratio of (MSbar light down-type mass in decoupled theory)/
                   (MSbar light down-type mass in full theory). 
     ze = ratio of (MSbar light charged lepton mass in decoupled theory)/
                   (MSbar light charged lepton mass in full theory). 

     Small mass corrections are not included here, but are added in 
     the functions SMDR_Decouple_b and SMDR_Decouple_c.
*/

void SMDR_QCDQEDmatch (const char *Fermion_type,
                       SMDR_REAL Q_match,
                       SMDR_REAL m_Fermion,
                       SMDR_REAL alphaS_hi,
                       SMDR_REAL alpha_hi,
                       int nqlight,
                       int loopOrder,
                       SMDR_REAL *alphaS_lo,
                       SMDR_REAL *alpha_lo,
                       SMDR_REAL *zu,
                       SMDR_REAL *zd,
                       SMDR_REAL *ze) 
{
  int nqlight2 = nqlight * nqlight;
  SMDR_REAL ZalphaS = 1; /* Decoupling coeff for alphaS */
  SMDR_REAL Zalpha = 1; /* Decoupling coeff for alpha */
  SMDR_REAL Zu = 1; /* Decoupling coeff for light up-type quark masses */
  SMDR_REAL Zd = 1; /* Decoupling coeff for light down-type quark masses */
  SMDR_REAL Ze = 1; /* Decoupling coeff for light charged lepton masses */

  SMDR_REAL lnbarF, lnbar2F, lnbar3F, lnbar4F;
  SMDR_REAL aS, aS2, aS3, aS4, a, a2;
  SMDR_REAL Zfcommonfactor, Zquark3loop, Zquark4loop;
  SMDR_REAL QF2, QF4; /* Powers of electric charge of decoupled heavy fermion */
  int quark; /* 1 (0) if the decoupled heavy fermion is a quark (lepton) */
  int NcF;   /* 3 (1) if the decoupled heavy fermion is a quark (lepton) */
  SMDR_REAL TF, CF;
  char funcname[] = "QCDQEDmatch";

  if ( !strncmp(Fermion_type, "u", 1) ) 
     {QF2 = 4./9.; quark = 1;}
  else if ( !strncmp(Fermion_type, "d", 1) ) 
     {QF2 = 1./9.; quark = 1;}
  else if ( !strncmp(Fermion_type, "e", 1) ) 
     {QF2 = 1.; quark = 0;}
  else SMDR_Error (funcname,
                   "Invalid heavy fermion type, should be u, d, or e.", 1);

  if ((loopOrder < 0) || (loopOrder > 4)) SMDR_Error (funcname,
    "Invalid loopOrder, should be 0, 1, 2, 3, or 4", 1);
    
  NcF = 1 + 2 * quark;
  TF = quark/2.;
  CF = quark * 4./3.;
  QF4 = QF2 * QF2;

  lnbarF = SMDR_LOG ((m_Fermion * m_Fermion)/(Q_match * Q_match));
  lnbar2F = lnbarF * lnbarF; 
  lnbar3F = lnbar2F * lnbarF; 
  lnbar4F = lnbar3F * lnbarF; 

  aS = alphaS_hi/(4.*PI);
  aS2 = aS * aS;
  aS3 = aS2 * aS;
  aS4 = aS3 * aS;

  a = alpha_hi/(4.*PI);
  a2 = a * a;

  /* 1812.04100 eqs. (3.4) and (3.6) */
  if (loopOrder > 0) {
    ZalphaS += (4./3.) * TF * aS * lnbarF;
    Zalpha += (4./3.) * NcF * QF2 * a * lnbarF;
  }

  /* 1812.04100 eqs. (3.5) and (3.5) and (3.11) */
  if (loopOrder > 1.999) {
    ZalphaS += (16./9.) * TF * TF * aS2 * lnbar2F 
               - TF * aS * (CF * aS + QF2 * a) * (4. * lnbarF + 13./3.)
               + TF * aS2 * (20.*lnbarF + 32./3.);
 
    Zalpha += (16./9.) * NcF * NcF * QF4 * a2 * lnbar2F 
              - NcF * QF2 * a * (CF * aS + QF2 * a) * (4. * lnbarF + 13./3.);

    Zfcommonfactor = 2. * lnbar2F + 10.*lnbarF/3. + 89./18.;
    Zu += (TF * (4./3.) * aS2 + NcF * QF2 * (4./9.) * a2) * Zfcommonfactor;
    Zd += (TF * (4./3.) * aS2 + NcF * QF2 * (1./9.) * a2) * Zfcommonfactor;
    Ze += (NcF * QF2 * a2) * Zfcommonfactor;
  }

  
  /* 1812.04100 eqs. (3.8) and (3.13) */
  if ((1 == quark) && (loopOrder > 2.999)) {
    ZalphaS += aS3 * (8.*lnbar3F/27. 
               + (53./9. - 16.*nqlight/9.) * lnbar2F
               + (955./9. - 67.*nqlight/9.) * lnbarF 
               + 62.211628 - 2633.*nqlight/486.);

    Zquark3loop = aS3 * ( (16.*nqlight/27. - 232./27.) * lnbar3F 
                  + (700./27.) * lnbar2F
                  + (212.*nqlight/27. + 71.788714) * lnbarF 
                  + 118.248112 + 1.582567 * nqlight);
    Zu += Zquark3loop;
    Zd += Zquark3loop;
  }

  /* 1812.04100 eqs. (3.9) and (3.14) */
  if ((1 == quark) && (loopOrder > 3.999)) {
    ZalphaS += aS4 * (16.*lnbar4F/81.
      + (3766./81. + 508.*nqlight/81. - 64.*nqlight2/81.) * lnbar3F
      + (4354./27. - 2966.*nqlight/81. - 77.*nqlight2/81.) * lnbar2F
      + (2157.863053 - 335.316171 * nqlight - 6865.*nqlight2/729.) * lnbarF
      + 1323.608830 - 258.542470 * nqlight - 5.62646 * nqlight2);

    Zquark4loop = aS4 * ( 
      (8.*nqlight2/27. - 80.*nqlight/9. + 610./9.) * lnbar4F 
      + (184.*nqlight/9. - 19264./81.) * lnbar3F
      + (496.*nqlight2/81. - 15650./81.*nqlight + 269.583577) * lnbar2F
      + (286.364218 + 39.625147 * nqlight - 1.284061 * nqlight2) * lnbarF
      + 14.375890 * nqlight2 - 375.221169 * nqlight + 1753.616640);

    Zu += Zquark4loop;
    Zd += Zquark4loop;
  }

  *alphaS_lo = alphaS_hi * ZalphaS;
  *alpha_lo = alpha_hi * Zalpha;
  *zu = Zu;
  *zd = Zd;
  *ze = Ze;
  
  return;
}

