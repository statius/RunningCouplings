// This program calculates and prints:

//   1) the pole masses of the Higgs boson, W and Z bosons, and the
//      top quark, using the state-of-the-art approximations for
//      each. This means full 2-loop plus leading 3-loop order for the
//      Higgs, full 2-loop for W and Z, and 4-loop pure QCD plus full
//      2-loop for the top quark.  Also given for the W and Z bosons
//      are the Breit-Wigner masses as usually quoted by
//      experimentalists:
//        MW_BreitWigner = sqrt(MW^2 + GammaW^2) and 
//        MZ_BreitWigner = sqrt(MZ^2 + GammaZ^2).

//   2) The MSbar quantities alpha_S(Q), alpha(Q), and sin^2_W(Q)
//      at Q = MZ, in the full Standard Model with nothing decoupled.

//   3) The MSbar quantities alpha_S(Q), alpha(Q), and sin^2_W(Q), 
//      at Q = MZ, in the PDG scheme where the top-quark is decoupled 
//      but the W boson is active.

//   4) The MSbar quantities alpha_S, alpha, m_b, m_c, m_s, m_d, m_u,
//      m_tau, m_muon, and m_electron, at the renormalization scale
//      Q=M_Z, in the 5-quark + 3-lepton QCD+QED effective theory.

//   5) The MSbar heavy quark masses m_b(m_b) in the 5-quark, 3-lepton
//      theory and m_c(m_c) in the 4-quark, 2-lepton theory.

//   6) The MSbar light quark masses m_s, m_d, m_u at Q = 2 GeV, in
//      the 4-quark, 3-lepton theory.

//   7) The charged lepton pole masses M_tau, M_muon, and M_electron.

//   8) The Sommerfeld fine structure constant alpha.

//   9) The Fermi constant GFermi.

// The inputs are the MSbar quantities in the non-decoupled Standard
// Model consisting of the gauge couplings g_3, g, g', the Higgs
// self-coupling lambda, the quark and lepton Yukawa couplings y_t,
// y_b, y_c, y_s, y_d, y_u, y_tau, y_mu, y_e, the Higgs vacuum
// expectation value v (defined to be the minimum of the effective
// potential in Landau gauge), and the non-perturbative light quark
// contribution to the fine-structure contant, Delta_hadronic^(5) alpha(MZ).

// Input values are read from the file ReferenceModel.dat which should
// be found in the current directory, unless a different file is
// specified using the -i option, see below. These inputs can then be
// modified if the -int (interactive mode) option is specified.

// Dimensionful quantities are in GeV. Timing information for the
// computation is also printed.

// This program takes four optional command line arguments:

//   -V                    Also finds Lambda and m2 from requiring the
//                         effective potential to be minimized with
//                         VEV v, and prints them. This doesn't affect
//                         the observables, and is much slower because
//                         many three-loop integrals have to be
//                         computed.

//   -i <input_filename>   Reads model data from <input_filename>; if
//                         not specified, data are read from the file
//                         ReferenceModel.dat.

//   -o <output_filename>  Writes complete model data to
//                         <output_filename>; if not specified, no
//                         output file is written and only the
//                         standard output mentioned above appears. If
//                         this flag is used, then -V is automatically
//                         also considered to be set. The output file
//                         is in a format that can be read into other
//                         programs (such as calc_fit or calc_RGrun)
//                         using SMDR_Read_Values ().

//   -int                  Runs the program in interactive mode. For
//                         each input parameter, the user is prompted
//                         to either accept the default value, or
//                         enter a new value.  The default values are
//                         either the ones from the input file
//                         specified by the -i flag, or else the ones
//                         found in ReferenceModel.dat.

// The program can be compiled from the current directory as, e.g.,

//    g++ -o calc_all calc_all.cpp -I.. -L.. -lsmdr -l3vil -ltsil -lm

// Run as, for example:
//    ./calc_all
// or
//    ./calc_all -i SampleModel.dat -V
// or
//    ./calc_all -i SampleModel.dat -o SampleOut.dat
// or
//    ./calc_all -i SampleModel.dat -int -o SampleOut.dat

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <ctime>

using namespace std;

#include "smdr.h"

int main (int argc, char *argv[])
{
  char inputFile[50], outputFile[50];
  int do_Veff, interactiveMode;
  /* char funcname[] = "calc_all"; */

  /* Input variables to be read in: */
  int nVars = 16;
  char *varList[] = {
    (char *)"SMDR_Q_in",
    (char *)"SMDR_v_in",
    (char *)"SMDR_lambda_in",
    (char *)"SMDR_g3_in",
    (char *)"SMDR_g_in",
    (char *)"SMDR_gp_in",
    (char *)"SMDR_yt_in",
    (char *)"SMDR_yb_in",
    (char *)"SMDR_yc_in",
    (char *)"SMDR_ys_in",
    (char *)"SMDR_yu_in",
    (char *)"SMDR_yd_in",
    (char *)"SMDR_ytau_in",
    (char *)"SMDR_ymu_in",
    (char *)"SMDR_ye_in",
    (char *)"SMDR_Delta_alpha_had_5_MZ_in"
  };

  /* Define arguments: */
  int nargs = 4;
  // char *arglist[] = {(char *)"-V",(char *)"-i",(char *)"-o",(char *)"-int"};
  char *arglist[] = {"-V","-i","-o","-int"};
  char *argtype[] = {(char *)"toggle",(char *)"string",(char *)"string",(char *)"toggle"};
  void *argvar[] = {&do_Veff, inputFile, outputFile, &interactiveMode};

  /* Set default values for optional args: */
  do_Veff = 0;
  strcpy (inputFile, "ReferenceModel.dat");
  strcpy (outputFile, "NULL");
  interactiveMode = 0;
  
  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);
  if (0 != strcmp ("NULL", outputFile)) do_Veff = 1; 

  SMDR_Read_Values (inputFile, nVars, varList);

  if (interactiveMode)
    SMDR_Set_Values_Interactively (nVars, varList);

  SMDR_Load_Inputs ();

  cout << "\nINPUT PARAMETERS read";
  if (0 == interactiveMode)
    cout << " from \"" << inputFile << "\"";
  cout << ":\n\n";

  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  SMDR_Display_Delta_alpha_had5 ();

  SMDR_Start_Timer();

  /* SMDR_Lambda and SMDR_m2 don't actually affect anything else, so
     only compute them if they are asked for, since the
     three-loop effective potential integrals are a speed bottleneck. 
     Otherwise just set them to 0.
  */
  if (1 == do_Veff) {
    SMDR_Lambda = 0;
    SMDR_Lambda_in = SMDR_Lambda = -SMDR_Eval_Veffmin (-1, 3.5);
    SMDR_m2_in = SMDR_m2;
  } else {
    SMDR_m2_in = SMDR_m2 = 0;
    SMDR_Lambda_in = SMDR_Lambda = 0;
  }

  /* Top-quark pole mass. */
  SMDR_Eval_Mt_pole (SMDR_Mt_EXPT, 1, 4, 2, &SMDR_Mt_pole, &SMDR_Gammat_pole);

  /* Higgs boson pole mass. */
  SMDR_Eval_Mh_pole (160., 2.5, &SMDR_Mh_pole, &SMDR_Gammah_pole);

  /* Z boson pole and Breit-Wigner masses */
  SMDR_Eval_MZ_pole (160., 2, &SMDR_MZ_pole, &SMDR_GammaZ_pole,
                     &SMDR_MZ_BreitWigner, &SMDR_GammaZ_BreitWigner);

  /* W boson pole and Breit-Wigner masses */
  SMDR_Eval_MW_pole (160., 2, &SMDR_MW_pole, &SMDR_GammaW_pole,
                     &SMDR_MW_BreitWigner, &SMDR_GammaW_BreitWigner);

  /* GFermi */
  SMDR_GFermi = SMDR_Eval_GFermi (SMDR_Mt_EXPT, 2);

  /* MSbar alpha(MZ), sin^2(thetaW), alphaS(MZ) in non-decoupled
     theory, Sommerfeld fine structure constant alpha, alpha(MZ) and
     sin^2(thetaW) in PDG MSbar scheme with only top decoupled.
  */
  SMDR_Eval_Gauge (SMDR_Mt_pole, SMDR_Mh_pole, SMDR_MW_BreitWigner);

  /* MSbar parameters in 5-quark, 3-lepton QCDQED theory at MZ */
  SMDR_Eval_QCDQED_at_MZ (SMDR_MZ_EXPT, SMDR_MZ_EXPT, 5);

  /* Running bottom mass evaluated at itself, 5-quark theory */
  SMDR_mbmb = SMDR_Eval_mbmb (SMDR_MZ_EXPT, 5);

  /* Running charm mass evaluated at itself, 4-quark theory */
  SMDR_mcmc = SMDR_Eval_mcmc (SMDR_Mtau_EXPT, SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 5);

  /* Running up, down, strange masses evaluated at Q = 2 GeV,  4-quark theory */
  SMDR_Eval_mquarks_2GeV (SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 5,
                          &SMDR_ms_2GeV, &SMDR_mu_2GeV, &SMDR_md_2GeV);

  /* Tau lepton pole mass. */
  SMDR_Mtau_pole = SMDR_Eval_Mtau_pole (SMDR_Mtau_EXPT,
                                        SMDR_mbmb_EXPT,
                                        SMDR_MZ_EXPT, 3);

  /* Muon pole mass. */
  SMDR_Mmuon_pole = SMDR_Eval_Mmuon_pole (SMDR_mcmc_EXPT,
                                          SMDR_mcmc_EXPT,
                                          SMDR_Mtau_EXPT,
                                          SMDR_mbmb_EXPT,
                                          SMDR_MZ_EXPT, 2);

  /* Electron pole mass. */
  SMDR_Melectron_pole = SMDR_Eval_Melectron_pole (SMDR_mcmc_EXPT,
                                                  SMDR_mcmc_EXPT,
                                                  SMDR_Mtau_EXPT,
                                                  SMDR_mbmb_EXPT,
                                                  SMDR_MZ_EXPT, 2);

  printf("\nOUTPUT QUANTITIES:\n");
 
  if (1 == do_Veff) {
    cout << endl;
    SMDR_Display_m2 ();
    SMDR_Display_Lambda ();
  }

  cout << endl;
  cout << "Mt = " << SMDR_Mt_pole << "; Gammat = " << SMDR_Gammat_pole 
       << ";   (* complex pole *)" << endl << endl;

  cout << "Mh = " << SMDR_Mh_pole << "; Gammah = " << SMDR_Gammah_pole 
       << ";   (* complex pole *)]" <<endl << endl;

  cout << "MZ = " << SMDR_MZ_pole << ";  GammaZ = " << SMDR_GammaZ_pole 
       << ";   (* complex pole *)" << endl;

  cout << "MZ = " << SMDR_MZ_BreitWigner << ";  GammaZ = " << SMDR_GammaZ_BreitWigner
       << ";   (* Breit-Wigner, compare to PDG *)" <<endl << endl;

  cout << "MW = " << SMDR_MW_pole << ";  GammaW = " << SMDR_GammaW_pole
       << ";   (* complex pole *)" << endl;

  cout << "MW = " << SMDR_MW_BreitWigner << ";  GammaW = " << SMDR_GammaW_BreitWigner
       << ";   (* Breit-Wigner, compare to PDG *)" <<endl << endl;

  cout << "MSbar quantities at Q = MZ, full Standard Model, nothing decoupled:" 
       << endl;
  cout << "  alphaS = " << SMDR_alphaS_MZ << "; "
       << "  alpha = 1/" << 1./SMDR_alpha_MZ << "; "
       << "  sin^2_thetaW = " << SMDR_s2W_MZ << ";" << endl << endl;

  cout << "MSbar quantities at Q = MZ, only top quark decoupled (PDG convention):"
       << endl;
  cout << "  alphaS = " << SMDR_alphaS_5_MZ << "; "
       << "  alpha = 1/" << 1./SMDR_alpha_MZ_PDG << "; "
       << "  sin^2_thetaW = " << SMDR_s2W_MZ_PDG << ";" << endl;

  cout << "\nMSbar quantities at Q = MZ, with top, Higgs, W, and Z all decoupled" 
       << endl;
  cout << "(5 quark + 3 lepton QCD+QED effective theory):" << endl;
  cout << "  alpha_S = " << SMDR_alphaS_5_MZ << ";   " 
       << "  alpha = 1/" << 1./SMDR_alpha_5_MZ << ";" << endl;

  cout << "  mb(MZ) = " << SMDR_mb_MZ << ";    "
       << "  mtau(MZ) = " << SMDR_mtau_MZ << ";" << endl;

  cout << "  mc(MZ) = " << SMDR_mc_MZ << ";    "
       << "  ms(MZ) = " << SMDR_ms_MZ << ";    "
       << "  mmu(MZ) = " << SMDR_mmuon_MZ << ";" << endl;

  cout << "  md(MZ) = " << SMDR_md_MZ << ";   "
       << "  mu(MZ) = " << SMDR_mu_MZ << ";   "
       << "  me(MZ) = " << SMDR_melectron_MZ << ";"
       <<endl << endl;

  cout << "MSbar bottom and charm masses:" << endl;
  cout << "  mb(mb) = " << SMDR_mbmb 
       << "; (MSbar mass in 5-quark + 3-lepton QCD+QED theory)" << endl;
  cout << "  mc(mc) = " << SMDR_mcmc
       << "; (MSbar mass in 4-quark + 2-lepton QCD+QED theory)" <<endl << endl;

  cout << "Light quark MSbar masses (at Q = 2 GeV, in 4-quark + 3-lepton QCD+QED theory):" << endl;
  cout << "  ms = " << SMDR_ms_2GeV << ";   "
       << "  mu = " << SMDR_mu_2GeV << ";   "
       << "  md = " << SMDR_md_2GeV << ";" <<endl << endl;

  cout << "Lepton pole masses:" << endl;
  cout << "  Mtau = "  << SMDR_Mtau_pole << ";    "
       << "  Mmuon = " << SMDR_Mmuon_pole << ";    "
       << "  Melectron = " << SMDR_Melectron_pole
       << ";" << endl << endl;

  cout << "Sommerfeld fine structure and Fermi constants:" << endl;
  cout << "  alpha = 1/" << 1./SMDR_alpha << ";   "
       << "  GFermi = " << 100000 * SMDR_GFermi << " 10^-5;" << endl;

  // GFermi (Degrassi, Gambino, Giardino version)
  // SMDR_GFermi_DGG = SMDR_Eval_GFermi_DGG (SMDR_Mt_pole, 
  //                                         SMDR_Mh_pole,  
  //                                         SMDR_MW_BreitWigner);
  // cout << "        alternate version GFermi_DGG = "
  //      << 100000 * SMDR_GFermi_DGG << " 10^-5;" << endl;

  SMDR_Timer();
  cout << "\nTotal calculation time: " << SMDR_Time_Total << " seconds" << endl;

  if (0 != strcmp (outputFile, "NULL"))
    SMDR_Write_Model_File (outputFile);

  return 0;
}
