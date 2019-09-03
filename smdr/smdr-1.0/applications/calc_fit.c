/* 
   This program calculates all of the MSbar quantities of the Standard Model, 
   using as inputs:

     1) The real part of the top-quark pole mass.

     2) The real part of the Higgs boson pole masses.

     3) The real part of the Z boson Breit-Wigner mass.

     4) The strong coupling constant alpha_S in the 5-quark theory at
        Q=MZ.

     5) The Sommerfeld fine structure constant alpha = 1/137.0...

     6) The non-perturbative light quark contribution to the
        fine-structure contant, Delta_hadronic^(5) alpha(MZ).

     7) The Fermi constant G_Fermi.

     8) The MSbar bottom quark mass evaluated at itself, mb(mb), in
        the 5-quark, 3-lepton QCD+QED effective theory.

     9) The MSbar charm quark mass evaluated at itself, mc(mc), in the
        4-quark, 2-lepton QCD+QED effective theory.

     10, 11, 12) The MSbar strange, down, and up masses, evaluated at
        Q=2 GeV, in the 4-quark, 3-lepton QCD+QED effective theory.

     13, 14, 15) The tau, muon, and electron pole masses.

   These input values are read from the file "ReferenceModel.dat",
   unless a different file is specified using the argument -i; see
   below. These inputs can then be modified if the -int (interactive
   mode) option is specified.

   Regardless of the choice of input file, the MSbar quantities in
   "ReferenceModel.dat" are always used as the initial guesses, but
   are otherwise not relevant to calc_fit.

   The program computes (by iteration) and then prints the MSbar
   quantities in the non-decoupled Standard Model consisting of: the
   gauge couplings g_3, g, g', the Higgs self-coupling lambda, the
   quark and lepton Yukawa couplings y_t, y_b, y_c, y_s, y_d, y_u,
   y_tau, y_mu, y_e, and the Higgs vacuum expectation value v (defined
   to be the minimum of the full effective potential in Landau gauge),
   as well as the Higgs squared mass parameter m2 (which is negative)
   and the value of the field-independent vacuum energy Lambda,
   determined by requiring the total loop-corrected vacuum energy to
   vanish.

   Dimensionful quantities are in GeV. Timing information for the
   computation is also printed.

   This program takes five optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes complete model data to
                          <output_filename>; if not specified, the
                          MSbar quantities appear on stdout but no
                          output file is written. The output file is
                          in a format that can be read into other
                          programs (such as calc_all or calc_RGrun)
                          using SMDR_Read_Values ().

   -e <error_tolerance>   Specifies error tolerance; if not specified,
                          a default value of 1.0e-7 is used.

   -Q <Q_target>          The renormalization scale at which the MSbar
                          parameters should be determined; if not
                          specified, Q_target will be the top-quark
                          pole mass.

   -int                   Runs the program in interactive mode. For
                          each input parameter, the user is prompted
                          to either accept the default value, or enter
                          a new value.  The default values are either
                          the ones from the input file specified by
                          the -i flag, or else the ones found in
                          ReferenceModel.dat.

   The executable for this program is automatically created within the
   smdr directory by make. You can also compile it separately as:

      gcc -o calc_fit calc_fit.c -L. -lsmdr -l3vil -ltsil -lm

   assuming the necessary archive files libsmdr.a and lib3vil.a and
   libtsil.a and the header file smdr.h are all in the current
   directory.

   Run as, for example:

     ./calc_fit

   for default choices with input ReferenceModel.dat, or

     ./calc_fit -i SampleModelOS.dat -o SampleOutput.dat -e 1.0e-9 -Q 160

   or

     ./calc_fit -int -o SampleOutput.dat -e 1.0e-9 -Q 160

   When running, it is important that the file ReferenceModel.dat (as
   well as the input model file, if it is different) can be found in
   the current directory.
*/

#include "smdr.h"

int main (int argc, char *argv[])
{
  SMDR_REAL Q_target;
  SMDR_REAL ERROR_TOLERANCE;
  char inputFile[50], outputFile[50];
  int interactiveMode;
  /* char funcname[] = "calc_fit"; */
  
  /* Input variables to be read in: */
  int nVars = 15;
  char *varList[] = {
    "SMDR_Mt_pole",
    "SMDR_Mh_pole",
    "SMDR_MZ_BreitWigner",
    "SMDR_alphaS_5_MZ",
    "SMDR_alpha",
    "SMDR_Delta_alpha_had_5_MZ_in",
    "SMDR_GFermi",
    "SMDR_mbmb",
    "SMDR_mcmc",
    "SMDR_ms_2GeV",
    "SMDR_md_2GeV",
    "SMDR_mu_2GeV",
    "SMDR_Mtau_pole",
    "SMDR_Mmuon_pole",
    "SMDR_Melectron_pole"
  };

  /* Define arguments: */
  int nargs = 5;
  char *arglist[] = {"-i","-o","-Q","-e","-int"};
  char *argtype[] = {"string","string","real","real","toggle"};
  void *argvar[] = {inputFile, outputFile, &Q_target, &ERROR_TOLERANCE, 
                    &interactiveMode};

  /* Set default values for optional args: */
  ERROR_TOLERANCE = 1.e-7;
  strcpy (inputFile, "ReferenceModel.dat");
  strcpy (outputFile, "NULL");
  Q_target = -1;
  interactiveMode = 0;

  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);

  /* Print version information: */
  SMDR_Display_Version ();

  SMDR_Read_Values (inputFile, nVars, varList);

  if (interactiveMode)
    SMDR_Set_Values_Interactively (nVars, varList);

  printf("\nOn-shell quantities read");
  if (0 == interactiveMode) {
    printf(" from \"%s\"", inputFile);
  } 
  printf(":\n\n");

  SMDR_Display_OS_Inputs ();
  printf("\n");

  SMDR_Start_Timer ();

  if (Q_target < 0) Q_target = SMDR_Mt_pole;

  SMDR_Fit_Inputs (Q_target,
                   SMDR_alphaS_5_MZ,
                   SMDR_alpha,
                   SMDR_GFermi,
                   SMDR_MZ_BreitWigner,
                   SMDR_Mh_pole,
                   SMDR_Mt_pole,
                   SMDR_mbmb,
                   SMDR_mcmc,
                   SMDR_ms_2GeV,
                   SMDR_md_2GeV,
                   SMDR_mu_2GeV,
                   SMDR_Mtau_pole,
                   SMDR_Mmuon_pole,
                   SMDR_Melectron_pole,
                   SMDR_Delta_alpha_had_5_MZ_in,
                   ERROR_TOLERANCE);

  printf("\nAfter fit by iteration:\n");

  printf("\n");
  printf("MW = %Lf;  (* pole mass *)\n", SMDR_MW_pole);
  printf("MW = %Lf;  (* Breit-Wigner mass, compare to PDG *)\n",
         SMDR_MW_BreitWigner);

  printf("\nMSbar parameters:\n");
  SMDR_Display_MSbar_Parameters ();
  SMDR_Display_v ();
  SMDR_Display_m2 ();
  SMDR_Display_Lambda ();
  SMDR_Display_Delta_alpha_had5 ();

  SMDR_Timer ();
  printf("\nTotal time: %.2f seconds\n\n", SMDR_Time_Total);

  /* Write output file, if one is specified: */
  if (0 != strcmp (outputFile, "NULL"))
    SMDR_Write_Model_File (outputFile);

  return 0;
}
