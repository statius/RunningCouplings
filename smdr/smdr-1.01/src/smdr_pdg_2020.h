/* Experimental reference values from Review of Particle Properties 2020. */

#ifndef _PDG_2020_H_
#define _PDG_2020_H_

/* Fermi decay constant. Table 1.1 of 2020 RPP. */
SMDR_REAL SMDR_GFermi_EXPT =     0.000011663787;
SMDR_REAL SMDR_GFermi_EXPT_UNC = 0.000000000006;

/* Fine structure constant. Table 1.1 of 2020 RPP. */
SMDR_REAL SMDR_alpha_EXPT =     0.0072973525693; 
SMDR_REAL SMDR_alpha_EXPT_UNC = 0.0000000000011;

/* MSbar at Q=MZ, top decoupled. Equation (9.25) of 2020 RPP. */
SMDR_REAL SMDR_alphaS_MZ_EXPT =     0.1179;     
SMDR_REAL SMDR_alphaS_MZ_EXPT_UNC = 0.0010;

/* MSbar at Q=MZ, top decoupled. Table 10.2 of 2020 RPP. */
SMDR_REAL SMDR_s2W_MZ_EXPT =     0.23121;    
SMDR_REAL SMDR_s2W_MZ_EXPT_UNC = 0.00003;

/* MSbar at Q=MZ, top decoupled. From 1/(127.952 +- 0.009), 
   Section 10.2.2 of 2020 RPP, two lines before equation (10.9). */
SMDR_REAL SMDR_alpha_MZ_EXPT = 0.00781543; 
SMDR_REAL SMDR_alpha_MZ_EXPT_UNC = 0.00000055;

/* Hadronic contribution to \alpha. Section 10.2.2 of 2020 RPP, line
   after equation (10.9). */
SMDR_REAL SMDR_Delta_alpha_had_5_MZ_EXPT =     0.02766;
SMDR_REAL SMDR_Delta_alpha_had_5_MZ_EXPT_UNC = 0.00007;

/* Top-quark pole mass (from cross-section measurements). June 2020 pdgLive. */
SMDR_REAL SMDR_Mt_EXPT =   172.4;    
SMDR_REAL SMDR_Mt_EXPT_UNC = 0.7;

/* Higgs boson mass. June 2020 pdgLive. */
SMDR_REAL SMDR_Mh_EXPT =   125.10;
SMDR_REAL SMDR_Mh_EXPT_UNC = 0.14;

/* Z-boson Breit-Wigner mass. June 2020 pdgLive */
SMDR_REAL SMDR_MZ_EXPT =    91.1876; 
SMDR_REAL SMDR_MZ_EXPT_UNC = 0.0021;

/* W-boson Breit-Wigner mass. June 2020 pdgLive */
SMDR_REAL SMDR_MW_EXPT =    80.379;  
SMDR_REAL SMDR_MW_EXPT_UNC = 0.012;

/* Bottom-quark MSbar mass evaluated at itself. June 2020 pdgLive */
SMDR_REAL SMDR_mbmb_EXPT =        4.18;   
SMDR_REAL SMDR_mbmb_EXPT_UNC_hi = 0.03;
SMDR_REAL SMDR_mbmb_EXPT_UNC_lo = 0.02;

/* Charm-quark MSbar mass evaluated at itself. June 2020 pdgLive */
SMDR_REAL SMDR_mcmc_EXPT =        1.27;  
SMDR_REAL SMDR_mcmc_EXPT_UNC_hi = 0.02;
SMDR_REAL SMDR_mcmc_EXPT_UNC_lo = 0.02;

/* Strange-quark MSbar mass evaluated at 2 GeV. June 2020 pdgLive */
SMDR_REAL SMDR_ms_2GeV_EXPT =        0.093;   
SMDR_REAL SMDR_ms_2GeV_EXPT_UNC_hi = 0.011;
SMDR_REAL SMDR_ms_2GeV_EXPT_UNC_lo = 0.005;

/* Up-quark MSbar mass evaluated at 2 GeV. June 2020 pdgLive */
SMDR_REAL SMDR_mu_2GeV_EXPT =        0.00216;  
SMDR_REAL SMDR_mu_2GeV_EXPT_UNC_hi = 0.00049;
SMDR_REAL SMDR_mu_2GeV_EXPT_UNC_lo = 0.00026;

/* Down-quark MSbar mass evaluated at 2 GeV. June 2020 pdgLive */
SMDR_REAL SMDR_md_2GeV_EXPT =        0.00467;  
SMDR_REAL SMDR_md_2GeV_EXPT_UNC_hi = 0.00048;
SMDR_REAL SMDR_md_2GeV_EXPT_UNC_lo = 0.00017;

/* Tau lepton pole mass. June 2020 pdgLive */
SMDR_REAL SMDR_Mtau_EXPT =     1.77686; 
SMDR_REAL SMDR_Mtau_EXPT_UNC = 0.00012;

/* Muon pole mass. June 2020 pdgLive */
SMDR_REAL SMDR_Mmuon_EXPT =     0.1056583745; 
SMDR_REAL SMDR_Mmuon_EXPT_UNC = 0.0000000024;

/* Electron pole mass. June 2020 pdfLive */
SMDR_REAL SMDR_Melectron_EXPT =     0.0005109989461; 
SMDR_REAL SMDR_Melectron_EXPT_UNC = 0.0000000000031;

#endif
