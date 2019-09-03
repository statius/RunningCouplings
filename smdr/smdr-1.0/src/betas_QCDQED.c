/*
   Sources: 

   For 4-loop and 5-loop alpha_S beta function:
     Baikov, Chetyrkin, Kuehn, 1606.08659 and 
     Herzog, Ruijl, Ueda, Vermaseren and Vogt 1701.01404.

   For 3-loop alpha and QED contributions to alpha_S, see Appendix A
   of SPM 0608026 and Appendix A of Bednyakov, Kniehl, Pikelner and
   Veretin 1612.00660, which has several typos.

   For three-loop QCD+QED contributions to the beta functions for
   fermion masses, see Appendix A of 0608026, and Appendix A of
   1612.00660, which has several typos.

   For 4-loop and 5-loop pure QCD contributions to the beta functions
   for fermion masses, see 1402.6611.
*/

/* 
   Beta functions for an effective SU(3)_C x U(1)_EM theory,
   applicable when the top quark, W boson, Z boson, and Higgs boson
   have been decoupled, along with possibly some subset of fermions b,
   c, s, u, d, tau, mu, e. The matter content of the theory,
   consisting of the remaining light fermions, is assumed to be:

     nu up-type quarks, 
     nd down-type quarks, 
     ne charged leptons, 

   For example, if only the t, W, Z, h have been decoupled, then
     nu=2; nd=3; ne=3.

   If b, t, W, Z, h have been decoupled, then
     nu=2; nd=2; ne=3.

   If tau, b, t, W, Z, h have been decoupled, then
     nu=2; nd=2; ne=2.

   etc.

   The beta functions here consist of:

     beta_alphaS, for MSbar alpha_S = g3^2/(4 Pi)
     beta_alpha,  for MSbar alpha = e^2/(4 Pi)
     beta_lnmu = (beta_mu)/mu, for the running MSbar up-type quark masses 
     beta_lnmd = (beta_md)/md, for the running MSbar down-type quark masses
     beta_lnme = (beta_me)/me, for the running MSbar charged lepton masses

   Note that there are no Yukawa couplings or other non-gauge couplings in
   this effective theory, so the masses run homogenously, and so: 
   beta_lnmu are the same for the non-decoupled up and charm, 
   beta_lnmd are the same for the non-decoupled down, strange and bottom, and 
   beta_lnme are the same for the non-decoupled e, mu, and tau.

   The beta functions include all QCD and QED contributions at full
   3-loop order, as well as the pure 4-loop and 5-loop order QCD
   contributions for alphaS.
*/

#include "smdr_internal.h"

/* Because we work in terms of alpha's, the loop factor in this file
   is 1/(4 Pi).
*/
#define ONELOOPFAC (1/(4.*PI)) 
#define TWOLOOPFAC (1/(16.*PI*PI)) 
#define THREELOOPFAC (1/(64.*PI*PI*PI)) 
#define FOURLOOPFAC (1/(256.*PI*PI*PI*PI)) 
#define FIVELOOPFAC (1/(1024.*PI*PI*PI*PI*PI)) 
#define nq (nd + nu)
#define nq2 nq*nq
#define nq3 nq2*nq
#define nq4 nq3*nq

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

void SMDR_Betas_QCDQED (int loopOrder, int nu, int nd, int ne)
{ 
  SMDR_REAL alphaS_run2, alphaS_run3, alphaS_run4, alphaS_run5, alphaS_run6;
  SMDR_REAL alpha_run2, alpha_run3;
  SMDR_REAL beta_lnmq4, beta_lnmq5;
  char funcname[] = "SMDR_Betas_QCDQED";

  if ( (loopOrder < 1) || (loopOrder > 5) )
    SMDR_Error (funcname,
    "Invalid loop order specified, should be 1, 2, 3, 4, or 5.", 3);

  alphaS_run2 = alphaS_run * alphaS_run;
  alphaS_run3 = alphaS_run2 * alphaS_run;
  alphaS_run4 = alphaS_run3 * alphaS_run;
  alphaS_run5 = alphaS_run4 * alphaS_run;
  alphaS_run6 = alphaS_run5 * alphaS_run;
  alpha_run2 = alpha_run * alpha_run;
  alpha_run3 = alpha_run2 * alpha_run;
  /* alpha_run4 = alpha_run3 * alpha_run; */

  beta_alphaS = ONELOOPFAC * alphaS_run2 * (-22. + 4.*nq/3.);

  beta_alpha = ONELOOPFAC * alpha_run2 * (8.*(nd + 3.*ne + 4.*nu))/9.; 

  beta_lnmu = ONELOOPFAC * (-8.*alphaS_run - 8.*alpha_run/3.);
 
  beta_lnmd = ONELOOPFAC * (-8.*alphaS_run - 2.*alpha_run/3.);

  beta_lnme = ONELOOPFAC * (-6.*alpha_run);

  if (loopOrder > 1) {

    beta_alphaS += TWOLOOPFAC * alphaS_run2 * (
                   alphaS_run * (-204. + 76.*nq/3.) +
                   4. * alpha_run * (nd + 4.*nu)/9. );

    beta_alpha += TWOLOOPFAC * alpha_run2 * ( 
                  32.*alphaS_run*(nd + 4.*nu)/9. +
                  8. * alpha_run * (nd + 27.*ne + 16.*nu)/27. ); 

    beta_lnmu += TWOLOOPFAC * ( 
                 alphaS_run2*(-404./3. + 40.*nq/9.) 
                 - 32.*alpha_run*alphaS_run/9. 
                 + alpha_run2*(-16./27. + 80.*(nd + 3.*ne + 4.*nu)/81.) );

    beta_lnmd += TWOLOOPFAC * ( 
                 alphaS_run2*(-404./3. + 40.*nq/9.) 
                 - 8.*alpha_run*alphaS_run/9. 
                 + alpha_run2*(-1./27. + 20.*(nd + 3.*ne + 4.*nu)/81.));

    beta_lnme += TWOLOOPFAC * (alpha_run2*(-3. + 20.*(nd + 3.*ne + 4.*nu)/9.) );
  }

  if (loopOrder > 2) {

    beta_alphaS += THREELOOPFAC * alphaS_run2 * (
      alphaS_run2*(-2857. + 5033.*nq/9. - 325.*nq2/27.) +
      56.*alpha_run*alphaS_run*(nd + 4.*nu)/27. +
      alpha_run2*(-44.*(nd+4.*nu)*(nd+3.*ne+4.*nu)/243. - 2.*(nd+16.*nu)/81.) );

    beta_alpha += THREELOOPFAC * alpha_run2 * (
      alphaS_run2*(1000.*(nd + 4.*nu)/27. - 176.*nq*(nd + 4.*nu)/81.) 
      -32.*alpha_run*alphaS_run*(nd + 16.*nu)/81. 
      + alpha_run2*(-88.*(nd + 3.*ne + 4.*nu)*(nd + 27.*ne + 16.*nu)/729. -
         4.*(nd + 243.*ne + 64.*nu)/243.));
 
    beta_lnmu += THREELOOPFAC * (
      alphaS_run3*(-2498. + 292.3675511518382*nq + 3.45679012345679*nq2) +
      alphaS_run2*alpha_run*(-191.11111111111111 - 1.2887009409867445*nd - 
        6.932581541724756*nu) +
      alphaS_run*alpha_run2*(-101.92592592592592 - 4.621721027816504*nd + 
        1.1851851851851851*ne - 18.486884111266015*nu) + 
      alpha_run3*(-11.325102880658436 - 0.28637798688594324*nd +
        - 10.892699473080961*ne - 6.162294703755339*nu +
        0.5121170553269319*nd*nd + 3.072702331961591*nd*ne + 
        4.609053497942387*ne*ne + 4.096936442615455*nd*nu + 
        12.290809327846365*ne*nu + 8.19387288523091*nu*nu));

    beta_lnmd += THREELOOPFAC * (
      alphaS_run3*(-2498. + 292.3675511518382*nq + 3.45679012345679*nq2) +
      alphaS_run2*alpha_run*(-47.77777777777778 - 1.733145385431189*nd - 
        7.377025986169201*nu) + 
      alphaS_run*alpha_run2*(-6.37037037037037 - 1.155430256954126*nd + 
        0.2962962962962963*ne - 4.621721027816504*nu) + 
      alpha_run3*(-0.17695473251028807 - 0.09628585474617717*nd + 
        - 1.6393391080376*nu - 2.7972489423443148*ne + 
        0.12802926383173296*nd*nd + 0.7681755829903978*nd*ne + 
        1.1522633744855968*ne*ne + 1.0242341106538637*nd*nu + 
        3.072702331961591*ne*nu + 2.0484682213077274*nu*nu)); 

    beta_lnme += THREELOOPFAC * (
      alphaS_run*alpha_run2*(-11.287761201476023)*(nd + 4.*nu)  + 
      alpha_run3*(-129. - 0.273980100123002*nd - 12.383681601968032*nu +
        - 23.397462703321054*ne + 1.1522633744855966*nd*nd + 
        6.91358024691358*nd*ne + 10.37037037037037*ne*ne +
        9.218106995884773*nd*nu + 27.65432098765432*ne*nu + 
        18.436213991769545*nu*nu));
  }

  if (loopOrder > 3) {
    /* From 1606.08659 eq. (4) and 1701.01404 eqs (3.6) and (3.8). */     
    beta_alphaS += FOURLOOPFAC * alphaS_run5 * (-58485.92827238825 + 
      13892.579234007108*nq - 810.1780809197257*nq2 - 2.998628257887517*nq3);
   
    /* From 1402.6611 eq. (3.3) */
    beta_lnmq4 = FOURLOOPFAC * alphaS_run4 * (-50659.02809866388 + 
      9783.020502341325*nq - 141.39522633210606*nq2 - 2.9661298454314915*nq3);

    beta_lnmu += beta_lnmq4;
    beta_lnmd += beta_lnmq4;
  }

  if (loopOrder > 4) {
    /* From 1606.08659 eq. (5) and 1701.01404 eqs (3.7) and (3.9). */
    beta_alphaS += FIVELOOPFAC * alphaS_run6 * (-1074295.3481404716 + 
      372323.8990286561*nq - 35135.51530687367*nq2 + 
      462.55534530227294*nq3 + 3.6849488162478052*nq4);

    /* From 1402.6611 eq. (3.4)	*/
    beta_lnmq5 = FIVELOOPFAC * alphaS_run5 * (-1146279.7224661442 + 
      294269.88474770344*nq - 15323.91060923213*nq2 - 221.835929931524*nq3 +
      0.17481503204489446*nq4);

    beta_lnmu += beta_lnmq5;
    beta_lnmd += beta_lnmq5;
  }
    
  return;
}  
