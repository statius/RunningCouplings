                     ********************************
		           Welcome to SMDR v1.0
                     ********************************

Copyright (C) 2019 S.P. Martin and D.G. Robertson

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.  See the file LICENSE.txt for further
details.

Contents of this file:

I.   Overview
II.  Building SMDR
III. Using SMDR
IV.  The SMDR API
V.   Sample Program
VI.  C++ Linkage

*********************************************************************  
I. Overview
********************************************************************* 

SMDR (/S/tandard /M/odel in /D/imensional /R/egularization) is a library 
of utilities for calculations involving the full Standard Model of 
particle physics. It provides a complete mapping between the MSbar 
lagrangian parameters and the observable quantities to which they most 
closely correspond. SMDR also contains utilities for minimizing the 
effective potential to determine the relation between the Higgs vacuum 
expectation value and the Lagrangian mass squared parameter, for 
renormalization group running of all parameters, and for matching onto 
low energy effective theories with various heavy particles decoupled. 
All computations are carried out at state-of-the-art precision. The code 
is written in C, and may be linked from C or C++.

SMDR requires:
* the TSIL library of utilities for the numerical calculation of
  dimensionally regularized 2-loop self-energy integrals [MR05]; and
* the 3VIL library for the evaluuation of 3-loop vacuum integrals
  [MR16].
The latest versions of these are always included with the SMDR
distribution, so typically users shouldn't have to worry about them.

Authors:
S.P. Martin [spmartin AT niu.edu] 
Department of Physics
Northern Illinois University
DeKalb, IL 60115
USA

D.G. Robertson [drobertson AT otterbein.edu]
Department of Physics
Otterbein University
Westerville, OH 43081
USA

To cite this program, please reference the paper (referred to in this
document as [MR19]):

  [MR19] "Standard Model parameters in the tadpole-free pure MSbar scheme" 
         by S.P. Martin and D.G. Robertson,
         [arXiv:1907.02500].

Other directly relevant papers are:

  [SM17] "Effective potential at three loops" Stephen P. Martin,
  	 Phys.Rev. D96 (2017) no.9, 096005 e-Print: arXiv:1709.02397
  	 [hep-ph]

  [MR16] "Evaluation of the general 3-loop vacuum Feynman integral",
  	 by S.P. Martin and D.G. Robertson Phys. Rev. D95 (2017) no.1,
  	 016008 [arXiv:1610.07720]

  [SM16] "Top-quark pole mass in the tadpole-free MS-bar scheme"
  	 Stephen P. Martin Phys.Rev. D93 (2016) no.9, 094017 e-Print:
  	 arXiv:1604.01134 [hep-ph]

[SPM15a] "Pole mass of the W boson at two-Loop order in the pure
         MS-bar scheme", by Stephen P. Martin, Phys.Rev. D91 (2015)
         no.11, 114003 e-Print: arXiv:1503.03782 [hep-ph]

[SPM15b] "Z-boson pole mass at two-loop order in the pure MS-bar
  	 scheme" Stephen P. Martin Phys.Rev. D92 (2015) no.1, 014026
  	 e-Print: arXiv:1505.04833 [hep-ph]

[SPM15c] "Four-Loop Standard Model effective potential at leading
  	 order in QCD" Stephen P. Martin Phys.Rev. D92 (2015) no.5,
  	 054029 e-Print: arXiv:1508.00912 [hep-ph]

  [MR14] "Higgs boson mass in the Standard Model at two-loop order and 
         beyond", by S.P. Martin and D.G. Robertson, Phys.Rev. D90
         (2014) no.7, 073010 [arXiv:1407.4336].

  [MR05] "TSIL: a program for the calculation of two-loop self-energy
         integrals", by S.P. Martin and D.G. Robertson,
         Comp. Phys. Comm. 174 (2006) 133, [hep-ph/0501132].

Also, the following papers contain results directly used in the
Standard Model effective potential minimization routines:

 [FJJ92] "The Standard model effective potential at two loops," by
         C. Ford, I. Jack and D.R.T. Jones, Nucl. Phys. B387, 373
         (1992) [Erratum-ibid. B504, 551 (1997)], [hep-ph/0111190]
 
 [SPM13] "Three-loop Standard Model effective potential at leading
         order in strong and top Yukawa couplings," by S.P. Martin,
         Phys. Rev. D 89, 013003 (2014), arXiv:1310.7553 [hep-ph].

 [SPM14] "Taming the Goldstone contributions to the effective
         potential," by S.P. Martin, arXiv:1406.2355 [hep-ph].

The multi-loop renormalization group running and decoupling routines
in SMDR use results listed in [MR19]. If your work relies on these,
you should cite the corresponding works.

SMDR is available from:

       http://www.niu.edu/spmartin/SMDR
       http://faculty.otterbein.edu/drobertson/SMDR

Version number: 1.0


********************************************************************* 
II. Building SMDR
********************************************************************* 

SMDR can be compiled on any system that supports the GNU Compiler
Collection (gcc), the Intel C compiler (icc), or a similar C compiler
with support for complex mathematics.

To compile SMDR, edit the Makefile in the main SMDR directory and
choose:

1. Compiler and optimization flags.

   Several sets are pre-defined in the Makefile; simply uncomment the
   appropriate one for your system, if present.  SMDR is currently
   known to compile with gcc (under Linux or Mac OS X).  Other C
   compilers should work provided that complex mathematics is
   supported, but in this case you will need to explicitly set the
   compiler name and optimization flags.

   If you succeed in building SMDR on a new platform, the authors
   would be grateful to know it, and to learn of any special measures
   that were needed to compile it.

2. Install directories, if desired.

   You can set INSTALL_INCS and INSTALL_LIBS to point to directories
   where you would like to place the header files and libraries,
   respectively, after compilation. If INSTALL_INCS and INSTALL_LIBS
   are not set, then the libraries and header files will be left after
   compilation in the main SMDR directory; they can then be moved by
   hand to any appropriate place.  Standard directories that are
   automatically searched by compilers and linkers typically include
   /usr/lib and /usr/include, but you will need root access to write
   to these directories.  If you specify other directories not on the
   standard search path, note that when compiling your own code it
   will be necessary to specify these directories using the options
   -I<dir> and -L<dir>.  See the compiler/linker man pages for
   complete details.

Users should not need to edit the Makefiles in any of the
subdirectories of the main SMDR installation (./src, ./applications,
and the TSIL and 3VIL directories).

Once these choices have been made, simply type

   make

to build SMDR, along with a suite of associated sample programs.  

Note that as a first step in the build process, make unpacks the
distributions of TSIL and 3VIL provided with SMDR, and compiles these.
It is strongly recommended that you use the very latest versions of
these packages, which will always be supplied automatically with SMDR.

After this compilation is complete, you can optionally type

   make install

to install the library and header files in the specified locations.
Congratulations, SMDR is ready for action!

As a quick introduction to the capabilities of SMDR, we recommend
starting by trying the command line executables calc_all and calc_fit
in interactive mode. These can be run (with self-explanatory prompted
input by the user) by issuing the commands

   ./calc_all -int

and

   ./calc_fit -int

In addition to calc_all and calc_fit, there are other command-line
executables to illustrate renormalization group running and effective
potential minimization:

       calc_RGrun
       calc_m2
       calc_vev

There are also executables that generate the data for all figures
appearing in [MR19]:

       fig_GFermi_vs_Q
       fig_MW_vs_Q
       fig_MZ_vs_Q
       fig_Mh_vs_Q
       fig_Mt_vs_Q
       fig_RGrun_QCDQED
       fig_RGrun_vs_Q
       fig_lambda_m2_vs_Q
       fig_lambda_vs_Q
       fig_m2_vs_Q
       fig_vev_vs_Q

All of these programs make use of a data file, ReferenceModel.dat, which
contains a benchmark set of MSbar and on-shell parameters. The
programs will require that this file be present in the directory from
which they are run. Users can edit this file, but be aware that it
will be overwritten with benchmark_data/ReferenceModel.dat each time
that 'make' is run. Users wishing to make their own persistent version
of this file should give it a different name, and then use the '-i'
option to the various applications to read it instead of
ReferenceModel.dat.

These programs are described in detail below.

Besides the sample programs and data files, the end product intended
for the user interested in writing their own programs to make use of
SMDR consists of the files:

       libsmdr.a     The static SMDR archive (will be placed in
		     INSTALL_LIBS upon "make install")

       smdr.h        The SMDR header file; must be included in any 
		     user code that uses SMDR (will be placed in
		     INSTALL_INCS upon "make install")

       libtsil.a     The static TSIL archive (will be placed in
		     INSTALL_LIBS upon "make install")

       tsil.h        The TSIL header file, will be included
		     automatically when using smdr.h (will be placed
		     in INSTALL_INCS upon "make install")

       lib3vil.a     The static 3VIL archive (will be placed in
		     INSTALL_LIBS upon "make install")

       3vil.h        The 3vil header file, will be included
		     automatically when including smdr.h (will be
		     placed in INSTALL_INCS upon "make install")


********************************************************************* 
III. Using SMDR
*********************************************************************

In SMDR, all dimensionful quantities are given in GeV.

To use SMDR functions in your code, you must:

1. Include the header file smdr.h in any source file that makes use
   of SMDR functions, e.g., by adding the line

	#include "smdr.h"

   This is appropriate if the file smdr.h is located in the directory
   where the code is being compiled; if it has been placed in a
   standard location such as /usr/include, then

	#include <smdr.h>

   would work.  If it is a nonstandard directory <inc_dir>, the
   compiler option

        -I<inc_dir>

   will generally be necessary.  See the compiler man pages on your
   system for further details.

   Note that smdr.h itself #includes tsil.h and 3vil.h, so there is no
   need to do this separately.

   #including smdr.h will define for you a set of global variables
   that represent the couplings of the Standard Model. These variable
   names are prefixed by SMDR_ so that they should not collide with
   other objects in users' code.  These are:

   	 SMDR_Q;      	    MSbar renormalization scale
	 SMDR_g3;     	    SU(3) gauge coupling
	 SMDR_gp;     	    U(1)_Y gauge coupling
	 SMDR_g;      	    SU(2) gauge coupling
	 SMDR_yt;     	    Top Yukawa coupling
	 SMDR_yb;     	    Bottom Yukawa coupling
	 SMDR_yc;     	    Charm Yukawa coupling
	 SMDR_ys;     	    Strange Yukawa coupling
	 SMDR_yu;     	    Up Yukawa coupling
	 SMDR_yd;     	    Down Yukawa coupling
	 SMDR_ytau;   	    Tau Yukawa coupling
	 SMDR_ymu;    	    Muon Yukawa coupling
	 SMDR_ye;     	    Electron Yukawa coupling
	 SMDR_lambda; 	    Higgs self coupling
	 SMDR_m2;     	    Higgs lagrangian mass^2 parameter
	 SMDR_v;      	    Higgs VEV
	 SMDR_Lambda; 	    Field-independent vacuum energy
	 SMDR_Delta_alpha_had_5_MZ;     Non-perturbative hadronic 
                                        contribution to alpha

   All of these are real parameters of type SMDR_REAL. Typical tasks
   will use most of these parameters as inputs. Users may set their
   values freely, and SMDR functions use them in their computations.

   Input values of Q and the MSbar lagrangian parameters are stored in
   global "input variables" SMDR_*_in, where * can be Q, g3, g, gp,
   yt, yb, yc, ys, yu, yd, ytau, ymu, ye, lambda, or v. These can be
   changed by the user, but typically they stay fixed while various
   calculations are performed using the "working variables"
   SMDR_*. After setting input variables, e.g., by reading input
   file(s), users can call

      SMDR_Load_Inputs ();

   to copy the input variables into working variables. 

   In addition, global parameters are defined that represent couplings
   in various effective theories with heavy particles decoupled. For
   example, the parameters in the 5-quark, 3-lepton effective theory
   obtained by integrating out the top quark and W, Z bosons, include:

   	  SMDR_alphaS_53       MSbar QCD coupling
   	  SMDR_alpha_53        MSbar EM coupling
          SMDR_mb_53           MSbar bottom quark mass

   etc. See smdr.h for a full listing.

   There are also global variables representing on-shell observable
   output quantities. These include, notably:

          SMDR_GFermi;      
          SMDR_alpha;          Sommerfeld fine structure constant 
          SMDR_alphaS_5_MZ;    QCD coupling, t,h,Z,W all decoupled 
          SMDR_Mt_pole;
          SMDR_Mh_pole;
          SMDR_MZ_pole;
          SMDR_MZ_BreitWigner;
          SMDR_MW_pole;
          SMDR_MW_BreitWigner;
          SMDR_mbmb;
          SMDR_mcmc;
          SMDR_ms_2GeV;
          SMDR_mu_2GeV;
          SMDR_md_2GeV;
          SMDR_Mtau_pole;
          SMDR_Mmuon_pole;
          SMDR_Melectron_pole;

   Again, see smdr.h for a full listing.

   Finally, there are constant experimental reference values taken from the
   Particle Data Group, and associated uncertainties.  These are also
   all of type SMDR_REAL:

          SMDR_GFermi_EXPT;
          SMDR_GFermi_EXPT_UNC;
          SMDR_alpha_EXPT;          /* Fine structure constant*/
          SMDR_alpha_EXPT_UNC;
          SMDR_s2W_MZ_EXPT;         /* MSbar at Q=MZ, top decoupled */
          SMDR_s2W_MZ_EXPT_UNC;
          SMDR_alphaS_MZ_EXPT;      /* MSbar at Q=MZ, top decoupled */
          SMDR_alphaS_MZ_EXPT_UNC;
          SMDR_alpha_MZ_EXPT;       /* MSbar at Q=MZ, top decoupled */
          SMDR_alpha_MZ_EXPT_UNC;
          SMDR_Delta_alpha_had_5_MZ_EXPT;
          SMDR_Delta_alpha_had_5_MZ_EXPT_UNC;
          SMDR_Mt_EXPT;             /* Pole mass, has renormalon ambiguity */
          SMDR_Mt_EXPT_UNC;
          SMDR_Mh_EXPT;
          SMDR_Mh_EXPT_UNC;
          SMDR_MZ_EXPT;             /* Experimental Breit-Wigner mass */
          SMDR_MZ_EXPT_UNC;
          SMDR_MW_EXPT;             /* Experimental Breit-Wigner mass */
          SMDR_MW_EXPT_UNC;
          SMDR_mbmb_EXPT;           /* MSbar mass evaluated at itself. */
          SMDR_mbmb_EXPT_UNC_hi;
          SMDR_mbmb_EXPT_UNC_lo;
          SMDR_mcmc_EXPT;           /* MSbar mass evaluated at itself. */
          SMDR_mcmc_EXPT_UNC_hi;
          SMDR_mcmc_EXPT_UNC_lo;
          SMDR_ms_2GeV_EXPT;        /* MSbar mass at Q = 2 GeV */
          SMDR_ms_2GeV_EXPT_UNC_hi;
          SMDR_ms_2GeV_EXPT_UNC_lo;
          SMDR_mu_2GeV_EXPT;        /* MSbar mass at Q = 2 GeV */
          SMDR_mu_2GeV_EXPT_UNC_hi;
          SMDR_mu_2GeV_EXPT_UNC_lo;
          SMDR_md_2GeV_EXPT;        /* MSbar mass at Q = 2 GeV */
          SMDR_md_2GeV_EXPT_UNC_hi;
          SMDR_md_2GeV_EXPT_UNC_lo;
          SMDR_Mtau_EXPT;           /* Pole mass */
          SMDR_Mtau_EXPT_UNC;
          SMDR_Mmuon_EXPT;          /* Pole mass */
          SMDR_Mmuon_EXPT_UNC;
          SMDR_Melectron_EXPT;      /* Pole mass */
          SMDR_Melectron_EXPT_UNC;

   At present these are set to values from the 2019 update to the
   2018 Review of Particle Proerties. Their values can be changed at run
   time, or updated in the file

          smdr_pdg.h.

2. Link to the SMDR, 3VIL, TSIL, and C math libraries at the end of
   the compilation process. This is accomplished via the (linker) flag

	-lsmdr -l3vil -ltsil -lm

   in that order.  (See examples in the SMDR Makefile.)  If the
   libraries are not in a standard location (including the case where
   they are in the current directory), you will generally need to add
   the flag

        -L<lib_dir>

   where <lib_dir> is the directory in which libtsil.a may be found.
   If this is the current directory, then

	 -L.

   may be used.  Again, consult the compiler man pages for complete
   details on making user libraries available to the linker.

Complete details regarding the SMDR functions are given in section IV
of this document.  In the rest of this section we will describe the
operation of the programs provided with the distribution. These can
also serve as potentially useful examples of using SMDR.

There are 5 programs that implement basic SMDR computations, and 11
others that generate the data plotted in the figures of [MR19].

----------------------------------
calc_all (Source file: calc_all.c)
----------------------------------

This program calculates and prints:

  1) The pole masses of the Higgs boson, W and Z bosons, and the top
     quark, using the state-of-the-art approximations for each. This
     means full 2-loop plus leading 3-loop order for the Higgs, full
     2-loop for W and Z, and 4-loop pure QCD plus full 2-loop for the
     top quark.  Also given for the W and Z bosons are the
     Breit-Wigner masses as usually quoted by experimentalists:
       MW_BreitWigner = sqrt(MW^2 + GammaW^2) and 
       MZ_BreitWigner = sqrt(MZ^2 + GammaZ^2).

  2) The MSbar quantities alpha_S(Q), alpha(Q), and sin^2_W(Q) at Q =
     MZ, in the full Standard Model with nothing decoupled.
  3) The MSbar quantities alpha_S(Q), alpha(Q), and sin^2_W(Q), at Q =
     MZ, in the PDG scheme where the top-quark is decoupled but the W
     boson is active.
  4) The MSbar quantities alpha_S, alpha, m_b, m_c, m_s, m_d, m_u,
     m_tau, m_muon, and m_electron, at the renormalization scale
     Q=M_Z, in the 5-quark + 3-lepton QCD+QED effective theory.
  5) The MSbar heavy quark masses m_b(m_b) in the 5-quark, 3-lepton
     theory and m_c(m_c) in the 4-quark, 2-lepton theory.
  6) The MSbar light quark masses m_s, m_d, m_u at Q = 2 GeV, in the
     4-quark, 3-lepton theory.
  7) The charged lepton pole masses M_tau, M_muon, and M_electron.
  8) The Sommerfeld fine structure constant alpha.
  9) The Fermi constant GFermi.

The inputs are the MSbar quantities in the non-decoupled Standard
Model consisting of the gauge couplings g_3, g, g', the Higgs
self-coupling lambda, the quark and lepton Yukawa couplings y_t, y_b,
y_c, y_s, y_d, y_u, y_tau, y_mu, y_e, the Higgs vacuum expectation
value v (defined to be the minimum of the effective potential in
Landau gauge), and the non-perturbative light quark contribution to
the fine-structure contant, Delta_hadronic^(5) alpha(MZ).

Input values are read from the file ReferenceModel.dat which should be
found in the current directory, unless a different file is specified
using the -i option; see below. These inputs can then be modified if
the -int (interactive mode) option is specified.

Dimensionful quantities are in GeV. Timing information for the
computation is also printed.

This program takes four optional command line arguments:

  -V                    Also finds Lambda and m2 from requiring the
                        effective potential to be minimized with
                        VEV v, and prints them. This doesn't affect
                        the observables, and is much slower because
                        many three-loop integrals have to be
                        computed.

  -i <input_filename>   Reads model data from <input_filename>; if
                        not specified, data are read from the file
                        ReferenceModel.dat.

  -o <output_filename>  Writes complete model data to
                        <output_filename>; if not specified, no
                        output file is written and only the
                        standard output mentioned above appears. If
                        this flag is used, then -V is automatically
                        also considered to be set. The output file
                        is in a format that can be read into other
                        programs (such as calc_fit or calc_RGrun)
                        using SMDR_Read_Values ().

  -int                  Runs the program in interactive mode. For
                        each input parameter, the user is prompted
                        to either accept the default value, or
                        enter a new value.  The default values are
                        either the ones from the input file
                        specified by the -i flag, or else the ones
                        found in ReferenceModel.dat.

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o calc_all calc_all.c -L. -lsmdr -l3vil -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and
lib3vil.a and the header files smdr.h and 3vil.h and tsil.h are all in
the current directory.

Run as, for example:

      ./calc_all

or

      ./calc_all -i SampleModel.dat -V

or

      ./calc_all -i SampleModel.dat -o SampleOut.dat

or

      ./calc_all -i SampleModel.dat -int -o SampleOut.dat


----------------------------------
calc_fit (Source file: calc_fit.c)
----------------------------------

This program calculates all of the MSbar quantities of the Standard
Model, using as inputs:

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

These input values are read from the file "ReferenceModel.dat", unless
a different file is specified using the argument -i; see below. These
inputs can then be modified if the -int (interactive mode) option is
specified.

Regardless of the choice of input file, the MSbar quantities in
"ReferenceModel.dat" are always used as the initial guesses, but are
otherwise not relevant to calc_fit.

The program computes (by iteration) and then prints the MSbar
quantities in the non-decoupled Standard Model consisting of: the
gauge couplings g_3, g, g', the Higgs self-coupling lambda, the quark
and lepton Yukawa couplings y_t, y_b, y_c, y_s, y_d, y_u, y_tau, y_mu,
y_e, and the Higgs vacuum expectation value v (defined to be the
minimum of the full effective potential in Landau gauge), as well as
the Higgs squared mass parameter m2 (which is negative) and the value
of the field-independent vacuum energy Lambda, determined by requiring
the total loop-corrected vacuum energy to vanish. It also computes and 
prints the resulting W boson mass.

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
libtsil.a and the header files smdr.h and 3vil.h and tsil.h are all in
the current directory.

Run as, for example:

     ./calc_fit

for default choices with input ReferenceModel.dat, or

     ./calc_fit -i SampleModelOS.dat -o SampleOutput.dat -e 1.0e-9 -Q 160

or

     ./calc_fit -int -o SampleOutput.dat -e 1.0e-9 -Q 160

When running, it is important that the file ReferenceModel.dat (or the 
input model file, if it is different) can be found in the current directory.


----------------------------------
calc_vev (Source file: calc_vev.c)
----------------------------------

This program calculates the Higgs VEV, obtained by minimizing the
Landau gauge effective potential with resummed Goldston boson
contributions, using, by default, the full 3-loop results of
arXiv:1709.02397 plus the 4-loop QCD contribution from 1508.00912.  It
takes one optional command-line inputs:

   -l <loopOrder>

The allowed values for loopOrder are:

   0   [tree-level]
   1   [1-loop]
   2   [2-loop, from 1406.2355 eq. (4.20), or SMDeltas.anc of
       1709.02397]
   2.5 [adds the leading 3-loop contributions from 1406.2355
       eq. (4.21)]
   3   [full 3-loop, from ancillary file SMDeltas.anc of 1709.02397]
   3.5.[full 3-loop plus the 4-loop contribution at leading order in
        QCD, from 1508.00912 equation (5.5)]

If -l <loopOrder> is omitted, it is taken to be 3.5 by default.  

The program then prompts the user for the quantities:
     Q, m2, lambda, g3, g, gp, yt, yb, ytau 
with default values taken from the file ReferenceModel.dat, which
should be in the current directory. It then calculates and prints the
full set of MSbar quantities including the value of the VEV v.

Dimensionful quantities are in appropriate powers of GeV.  Timing
information for the computation is also printed.

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o calc_vev calc_vev.c -L. -lsmdr -ltsil -l3vil -lm

assuming the necessary archive files libsmdr.a and lib3vil.a and
libtsil.a and the header files smdr.h and 3vil.h and tsil.h are all in
the current directory.

Run as, for example:

   ./calc_vev 

or 

   ./calc_vev -l 2 

for the result including only the 2-loop contributions.


--------------------------------
calc_m2 (Source file: calc_m2.c)
--------------------------------

This program calculates the negative Higgs lagrangian squared mass
parameter m2, obtained by minimizing the Landau gauge effective
potential with resummed Goldston boson contributions, using, by
default, the full 3-loop results of arXiv:1709.02397 plus the leading
order in QCD at 4-loop order from 1508.00912. It takes one optional
command-line input:

   -l <loopOrder>

The values for loopOrder are:
 
   0   [tree-level]
   1   [1-loop]
   2   [2-loop, from 1406.2355 eq. (4.20), or SMDeltas.anc of
        1709.02397]
   2.5 [adds the leading 3-loop contributions from 1406.2355
        eq. (4.21)]
   3   [full 3-loop, from ancillary file SMDeltas.anc of 1709.02397]
   3.5 [full 3-loop plus the 4-loop contribution at leading order in
        QCD, from 1508.00912 equation (5.5)]
   
If -l <loopOrder> loopOrder is omitted, it is taken to be 3.5 by
default.

The program then prompts the user for the quantities:
     Q, v, lambda, g3, g, gp, yt, yb, ytau

with default values taken from the file ReferenceModel.dat, which
should be in the current directory. It then calculates and prints the
full set of MSbar quantities including the value of m2.

Dimensionful quantities are in appropriate powers of GeV.  Timing
information for the computation is also printed.

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

     gcc -o calc_m2 calc_m2.c -L. -lsmdr -ltsil -l3vil -lm

assuming the necessary archive files libsmdr.a and lib3vil.a and
libtsil.a and the header files smdr.h and 3vil.h and tsil.h are all in
the current directory.

Run as, for example:

   ./calc_m2 

for the result with all known contributions, or

   ./calc_m2 -l 2

for the result including only the 2-loop contributions.


--------------------------------------
calc_RGrun (Source file: calc_RGrun.c)
--------------------------------------

This program implements renormalization group running of the Standard
Model MSbar parameters from one specified scale to another.  It takes
one required argument Qfinal and up to three optional arguments:

  Q_final               This required command line argument is the 
                        destination RG scale for the running.

  -i <input_filename>   Reads model data from <input_filename>; if not
                        specified, data are read from the file
                        ReferenceModel.dat, which should be in the
                        current directory.

  -int                  Runs the program in interactive mode. For each
                        input parameter, the user is prompted to
                        either accept the default value, or enter a
                        new value.  The default values are either the
                        ones from the input file specified by the -i
                        flag, or else the ones found in
                        ReferenceModel.dat.

  -l <loopOrder>        Runs the RG equations at the specified loop
                        order.  Allowed values for loopOrder are
                        1,2,3,4,5.  If not specified, by default
                        loopOrder = 5.

As output, the program prints the values of the MSbar quantities at
the scales Q = Q_in and Q_final. General and timing information for
the computation are also printed.

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

      gcc -o calc_RGrun calc_RGrun.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and lib3vil.a and
libtsil.a and the header files smdr.h and 3vil.h and tsil.h are all in
the current directory.

Run as, for example:

      ./calc_RGrun 2.4e18

for full 5-loop running and input data in ReferenceModel.dat, or

      ./calc_RGrun 2.4e18 -l 2 -i myInFile -int

for 2-loop running and model data in myInFile, or

      ./calc_RGrun 2.4e18 -i myInFile -int

for prompts that allow interactive modification of the initial
parameter values.


----------------------------------------
fig_Mh_vs_Q (Source file: fig_Mh_vs_Q.c)
----------------------------------------

Calculates the Higgs pole mass as a function of the MSbar
renormalization scale Q at which it is computed, in various
approximations, using the results in 1407.4336. The input parameters
are obtained from the file "ReferenceModel.dat" unless a different
file is specified using the "-i" option; see below.

This program produces by default an output file: "FIG_Mh_vs_Q.dat"
with data in 7 columns, respectively
 
   1  Q  (MSbar renormalization scale)
   2  Mh (loopOrder = 0, tree-level)
   3  Mh (loopOrder = 1, 1-loop) 
   4  Mh (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   5  Mh (loopOrder = 2, full 2-loop) 
   6  Mh (loopOrder = 2.3, includes 3-loop QCD) 
   7  Mh (loopOrder = 2.5, includes 3-loop QCD and yt corrections) 

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 4.1a graphs columns 3, 4, 5, and 7 as a function of column 1.

The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_Mh_vs_Q.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_Mh_vs_Q fig_Mh_vs_Q.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and the
header files smdr.h and tsil.h are all in the current directory.

Run as:

   ./fig_Mh_vs_Q

Running time is typically of order 4 minutes, depending on your
hardware, so you might want to run it in the background.


------------------------------------------------------
fig_lambda_m2_vs_Q (Source file: fig_lambda_m2_vs_Q.c)
------------------------------------------------------

Finds the Standard Model MSbar Higgs parameters lambda and m^2 at
their best fit to current experimental observables, and also the
envelope of values that lambda and m^2 take as the observables are
varied within 1 sigma and 2 sigma of their central values, as a
function of the renormalization scale from Q = MZ to Q = 10^19 GeV.
The current experimental values are all taken from the file
"smdr_pdg.h".

Output data are given by default in the file
"FIG_lambda_m2_highQ.dat", in 12 columns:

   1  Q
   2  Log[10,Q]
   3  lambda (lower bound for observables within 2 sigma)  
   4  lambda (lower bound for observables within 1 sigma)  
   5  lambda (central value)  
   6  lambda (upper bound for observables within 1 sigma)  
   7  lambda (upper bound for observables within 2 sigma)  
   8  sqrt(-m2) (lower bound for observables within 2 sigma)  
   9  sqrt(-m2) (lower bound for observables within 1 sigma)  
   10 sqrt(-m2) (central value)  
   11 sqrt(-m2) (upper bound for observables within 1 sigma)  
   12 sqrt(-m2) (upper bound for observables within 2 sigma)  

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 4.3a graphs columns 3, 4, 5, 6, and 7 as a function of column 2.
Figure 4.3b graphs columns 8, 9, 10, 11, and 12 as a function of
column 2.

The program takes an optional command line argument:

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_lambda_m2_highQ.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_lambda_m2_vs_Q fig_lambda_m2_vs_Q.c -L. -lsmdr -ltsil -l3vil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and
lib3vil.a and the header files smdr.h and 3vil.h and tsil.h are all in
the current directory.

Run as:

   ./fig_lambda_m2_vs_Q

Running in the background is highly recommended, because the running
time can be of order an hour, depending on your hardware.


----------------------------------------
fig_MW_vs_Q (Source file: fig_MW_vs_Q.c)
----------------------------------------

Calculates the physical mass of the W boson as a function of the MSbar
renormalization scale Q at which it is computed, in various
approximations, using the results in 1503.03782. The input parameters
are obtained from the file "ReferenceModel.dat" unless a different
file is specified using the "-i" option; see below.

This program produces by default an output file "FIG_MW_vs_Q.dat",
with data in 14 columns:

   1  Q  (MSbar renormalization scale)
   2  MW_pole            (loopOrder = 0, tree-level)
   3  MW_pole            (loopOrder = 1, 1-loop)
   4  MW_BreitWigner     (loopOrder = 1, 1-loop)
   5  GammaW_pole        (loopOrder = 1, 1-loop)
   6  GammaW_BreitWigner (loopOrder = 1, 1-loop)
   7  MW_pole            (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   8  MW_BreitWigner     (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   9  GammaW_pole        (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   10 GammaW_BreitWigner (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   11 MW_pole            (loopOrder = 2, full 2-loop)
   12 MW_BreitWigner     (loopOrder = 2, full 2-loop)
   13 GammaW_pole        (loopOrder = 2, full 2-loop)
   14 GammaW_BreitWigner (loopOrder = 2, full 2-loop)

where the complex pole squared mass is:

   s_pole = MW_pole^2 - i GammaW MW.

Note that the Breit-Wigner mass reported by experimentalists is
related by
 
   MW_BreitWigner^2 = MW_pole^2 + GammaW^2

so that MW_BreitWigner is approximately 27 MeV larger than MW.

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 4.1d graphs columns 4, 8, and 12 as a function of column 1.

The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_MW_vs_Q.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_MW_vs_Q fig_MW_vs_Q.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and the
header files smdr.h and tsil.h are all in the current directory.

Run as:

   ./fig_MW_vs_Q

The running time is typically of order 20 seconds, depending on your
hardware.


----------------------------------------
fig_MZ_vs_Q (Source file: fig_MZ_vs_Q.c)
----------------------------------------

Calculates the physical mass of the Z boson as a function of the MSbar
renormalization scale Q at which it is computed, in various
approximations, using the results in 1505.04833. The input parameters
are obtained from the file "ReferenceModel.dat" unless a different
file is specified using the "-i" option; see below.

This program produces by default the output file "FIG_MZ_vs_Q.dat",
with data in 14 columns:

   1  Q                  (MSbar renormalization scale)
   2  MZ_pole            (loopOrder = 0, tree-level)
   3  MZ_pole            (loopOrder = 1, 1-loop)
   4  MZ_BreitWigner     (loopOrder = 1, 1-loop)
   5  GammaZ_pole        (loopOrder = 1, 1-loop)
   6  GammaZ_BreitWigner (loopOrder = 1, 1-loop)
   7  MZ_pole            (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   8  MZ_BreitWigner     (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   9  GammaZ_pole        (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   10 GammaZ_BreitWigner (loopOrder = 1.5, 1-loop plus 2-loop QCD)
   11 MZ_pole            (loopOrder = 2, full 2-loop)
   12 MZ_BreitWigner     (loopOrder = 2, full 2-loop)
   13 GammaZ_pole        (loopOrder = 2, full 2-loop)
   14 GammaZ_BreitWigner (loopOrder = 2, full 2-loop)
 
where the complex pole squared mass is denoted:

   s_pole = MZ_pole^2 - i GammaZ_pole MZ_pole.

Note that the Breit-Wigner mass reported by experimentalists is
related by

   MZ_BreitWigner^2 = MZ_pole^2 + GammaZ_pole^2

so that MZ_BreitWigner is approximately 34 MeV larger than MZ.

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 4.1c graphs columns 4, 8, and 12 as a function of column 1.

The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_MZ_vs_Q.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_MZ_vs_Q fig_MZ_vs_Q.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and the
header files smdr.h and tsil.h are all in the current directory.

Run as:

   ./fig_MZ_vs_Q 

The running time is typically of order 15 seconds, depending on your
hardware.


----------------------------------------
fig_Mt_vs_Q (Source file: fig_Mt_vs_Q.c)
----------------------------------------

Calculates the top-quark pole mass as a function of the MSbar
renormalization scale Q at which it is computed, in various
approximations, using the results of 1604.01134 and references
therein. The input parameters are obtained from the file
"ReferenceModel.dat" unless a different file is specified using the
"-i" option; see below. This program produces output data in 12
columns:

   1    Q 
   2    Mt tree
   3    Mt 1-loop QCD
   4    Mt 2-loop QCD
   5    Mt 3-loop QCD
   6    Mt 4-loop QCD
   7    Mt     4-loop QCD + 1-loop nonQCD
   8    Gammat 4-loop QCD + 1-loop nonQCD
   9    Mt     4-loop QCD + 1-loop nonQCD + 2-loop mixed
   10   Gammat 4-loop QCD + 1-loop nonQCD + 2-loop mixed
   11   Mt     4-loop QCD + 1-loop nonQCD + 2-loop full  
   12   Gammat 4-loop QCD + 1-loop nonQCD + 2-loop full

Two output files are produced: "FIG_Mt_vs_Q_0.dat" and
"FIG_Mt_vs_Q_1.dat", using expansions around the tree-level and real
pole masses, respectively.

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 4.1b graphs columns 6, 7, 9, and 11 as a
function of column 1 from FIG_Mt_vs_Q_1.dat.

The program takes an optional command line argument:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          "ReferenceModel.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_Mt_vs_Q fig_Mt_vs_Q.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and the
header files smdr.h and tsil.h are all in the current directory.

Run as:

   ./fig_Mt_vs_Q
   
The running time is typically of order 6 minutes, depending on your
hardware, so you might want to run it in the background.


------------------------------------------------
fig_GFermi_vs_Q (Source file: fig_GFermi_vs_Q.c)
------------------------------------------------

Calculates GFermi as a function of the MSbar renormalization scale Q
at which it is computed, in various approximations. The input
parameters are obtained from the file "ReferenceModel.dat". This
program produces by default an output file "FIG_GFermi_vs_Q.dat", with
data in 11 columns, respectively

   1  Q  (MSbar renormalization scale) 
   2  GFermi 10^5 (loopOrder = 0, tree-level) 
   3  GFermi 10^5 (loopOrder = 1, 1-loop) 
   4  GFermi 10^5 (loopOrder = 1.5, 1-loop plus 2-loop LO in QCD) 
   5  GFermi 10^5 (loopOrder = 2, full 2-loop) 
   6, 7, 8, 9 = same, but divided by GFermi(experiment)

In the paper: arXiv:1907.02500 
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson, 
Figure 3.3 graphs columns 3, 4, and 5 as a function of column 1.

This program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_GFermi_vs_Q.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_GFermi_vs_Q fig_GFermi_vs_Q.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and the
header files smdr.h and tsil.h are all in the current directory.

Run as:

   ./fig_GFermi_vs_Q 


------------------------------------------------
fig_lambda_vs_Q (Source file: fig_lambda_vs_Q.c)
------------------------------------------------

Calculates the Higgs self coupling lambda as a function of the MSbar
renormalization scale Q. The inputs are the experimental Higgs pole
mass SMDR_Mh_EXPT from smdr_pdg.h and the other MSbar input parameters
obtained from the file "ReferenceModel.dat" unless a different file is
specified using the "-i" option; see below. This program produces (by
default) an output file "FIG_lambda_vs_Q.dat" containing data in 13
columns:

   1  Q (MSbar renormalization scale)
   2  lambda (tree-level)
   3  lambda (1-loop)
   4  lambda (1-loop plus 2-loop QCD)
   5  lambda (2-loop)
   6  lambda (2-loop plus leading 3-loop QCD)
   7  lambda (2-loop plus leading 3-loop QCD and yt)
   8  lambda/lambdarun (tree-level)
   9  lambda/lambdarun (1-loop)
   10 lambda/lambdarun (1-loop plus 2-loop QCD)
   11 lambda/lambdarun (2-loop)
   12 lambda/lambdarun (2-loop plus leading 3-loop QCD)
   13 lambda/lambdarun (2-loop plus leading 3-loop QCD and yt)
 
where lambda is the value obtained at Q from the input Higgs pole mass
with all other pertinent MSbar parameters having been run from their
benchmark values, and lambdarun is instead obtained by directly RG
running its benchmark input value to the scale Q.

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme 
by Stephen P. Martin and David G. Robertson,
Figure 4.2a graphs columns 3, 4, 5, and 7 as a function of column 1.
Figure 4.2b graphs columns 9, 10, 11, and 13 as a function of column 1.

The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          "ReferenceModel.dat".

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_lambda_vs_Q.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_lambda_vs_Q fig_lambda_vs_Q.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and the
header file smdr.h and tsil.h are all in the current directory.
 
Run as:

   ./fig_lambda_vs_Q

The running time is typically of order 3 minutes, depending on your
hardware.


----------------------------------------------
fig_RGrun_vs_Q (Source file: fig_RGrun_vs_Q.c)
----------------------------------------------

Calculates the Standard Model MSbar parameters as a function of the
renormalization scale from Q = 100 GeV to Q = 10^19 GeV, given their
input values at a benchmark point defined in the file
"ReferenceModel.dat" unless a different file is specified using the
"-i" option; see below.  The full state-of-the art beta functions are
used.

The data is output by default in the file "FIG_RGrun_vs_Q.dat" in 24
columns, respectively:

   1  Q
   2  Log[Q]
   3  g3 
   4  g 
   5  gp 
   6  alpha3 
   7  alpha2
   8  alpha1 
   9  1/alpha3 
   10 1/alpha2
   11 1/alpha1 
   12  yt 
   13  yb 
   14  yc 
   15  ys 
   16  yd 
   17  yu 
   18  ytau 
   19  ymu 
   20  ye 
   21  lambda 
   22  sqrt(-m2) 
   23  v
   24  Lambda

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 2.1a graphs columns 9, 10, and 11 as a function of column 2.
Figure 2.1b graphs columns 12, 13, 14, 15, 16, 17, 18, 19, and 20 
as a function of column 2.

The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

   -o <output_filename>   Writes output data to <output_filename>; if
                          not specified, the results appear in
                          "FIG_RGrun_vs_Q.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_RGrun_vs_Q fig_RGrun_vs_Q.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and the
header files smdr.h and tsil.h are all in the current directory.

Run as:

   ./fig_RGrun_vs_Q


--------------------------------------------------
fig_RGrun_QCDQED (Source file: fig_RGrun_QCDQED.c)
--------------------------------------------------

Prints the running MSbar parameters in the QCD+QED effective theory as
a function of the renormalization scale from Q = MZ down to Q = 1 GeV,
decoupling the bottom, tau, and charm quarks at appropriate
thresholds.  The input parameters come from the file
"ReferenceModel.dat" unless a different file is specified using the
"-i" option; see below. The state-of-the-art beta functions and
decoupling relations are used.

Data are output in the following files: 

"FIG_RGrun_53.dat" for the 5-quark, 3-lepton theory, in 11 columns,
     1  Q (from SMDR_MZ_EXPT down to SMDR_mbmb_EXPT)
     2  alphaS 
     3  alpha 
     4  mb 
     5  mc 
     6  ms 
     7  md
     8  mu
     9  mtau 
     10  mmuon 
     11  melectron 

"FIG_RGrun_43.dat" for the 4-quark, 3-lepton theory, in 10 columns,
     1  Q (from SMDR_mbmb_EXPT down to SMDR_Mtau_EXPT)
     2  alphaS 
     3  alpha 
     4  mc 
     5  ms 
     6  md
     7  mu
     8  mtau 
     9  mmuon 
     10  melectron 

"FIG_RGrun_42.dat" for the 4-quark, 2-lepton theory, in 9 columns,
     1  Q (from SMDR_Mtau_EXPT down to SMDR_mcmc_EXPT)
     2  alphaS 
     3  alpha 
     4  mc 
     5  ms 
     6  md
     7  mu
     8  mmuon 
     9  melectron 

"FIG_RGrun_32.dat" for the 3-quark, 2-lepton theory, in 8 columns,
     1  Q (from SMDR_mcmc_EXPT down to 1 GeV)
     2  alphaS 
     3  alpha 
     4  ms 
     5  md
     6  mu
     7  mmuon 
     8  melectron 

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 2.2a graphs columns 2 and 3 as a function of column 1 from
     FIG_RGrun_53.dat FIG_RGrun_43.dat FIG_RGrun_42.dat FIG_RGrun_32.dat
Figure 2.2b graphs 
     columns 4, 5, 6, 7, 8, 9, 10, and 11 as a function of column 1 
       from FIG_RGrun_53.dat
     columns 4, 5, 6, 7, 8, 9, and 10 as a function of column 1 
       from FIG_RGrun_43.dat
     columns 4, 5, 6, 7, 8, and 9 as a function of column 1 
       from FIG_RGrun_42.dat
     columns 4, 5, 6, 7, and 8 as a function of column 1 
       from FIG_RGrun_32.dat

The program takes an optional command line argument:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          ReferenceModel.dat.

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_RGrun_QCDQED fig_RGrun_QCDQED.c -L. -lsmdr -ltsil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and the
header files smdr.h and tsil.h are all in the current directory.

Run as: 

   ./fig_RGrun_QCDQED


----------------------------------------
fig_m2_vs_Q (Source file: fig_m2_vs_Q.c)
----------------------------------------

Calculates the Higgs lagrangian mass^2 parameter m2 as a function of
the MSbar renormalization scale Q. Here m2 is determined by
minimization of the effective potential, in various approximations
from:
      hep-ph/0111190 (full 1-loop and 2-loop).
      1310.7553 (3-loop at leading order in g3 and yt), 
      1709.02397 (full 3-loop), and 
      1508.00912 (4-loop at leading order in g3). 

The input parameters are obtained from the file "ReferenceModel.dat"
unless a different file is specified using the "-i" option; see below.
This program produces by default an output file "FIG_m2_vs_Q.dat",
with data in 13 columns:

   1  Q (MSbar renormalization scale)
   2  sqrt(-m2(0))   (loopOrder = 0, tree-level),
   3  sqrt(-m2(1))   (loopOrder = 1, 1-loop),
   4  sqrt(-m2(2))   (loopOrder = 2, 2-loop),
   5  sqrt(-m2(2.5)) (loopOrder = 2.5, includes leading 3-loop),
   6  sqrt(-m2(3))   (loopOrder = 3, full 3-loop),
   7  sqrt(-m2(3.5)) (loopOrder = 3, full 3-loop plus QCD 4-loop),
   8   m2(0)/m2run,
   9   m2(1)/m2run,
   10  m2(2)/m2run,
   11  m2(2.5)/m2run,
   12  m2(3)/m2run,
   13  m2(3.5)/m2run

where m2(n) is the value of m2 obtained by requiring the effective
potential to be minimized at the scale Q at loopOrder=n, and m2run is
obtained by directly running the input value of m2 to the scale Q
using 3-loop RG equations.

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 3.1a graphs columns 3, 4, 5, 6, and 7 as a function of column 1.
Figure 3.1b graphs columns 10, 11, 12, and 13 as a function of column 1.

The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          "ReferenceModel.dat".

   -o <output_filename>   Writes data to <output_filename>; if not
                          specified, the results appear in
                          "FIG_m2_vs_Q.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_m2_vs_Q fig_m2_vs_Q.c -L. -lsmdr -ltsil -l3vil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and
lib3vil.a and the header file smdr.h and tsil.h and 3vil.h are all in
the current directory.

Run as:

   ./fig_m2_vs_Q

Running time is typically of order 15 minutes, depending on your
hardware.


------------------------------------------
fig_vev_vs_Q (Source file: fig_vev_vs_Q.c)
------------------------------------------

Calculates the Higgs VEV as a function of the MSbar renormalization
scale Q.  Here v is determined by minimization of the effective
potential, in various approximations from:
      hep-ph/0111190 (full 1-loop and 2-loop).
      1310.7553 (3-loop at leading order in g3 and yt),
      1709.02397 (full 3-loop), and
      1508.00912 (4-loop at leading order in g3).

The input parameters are obtained from the file "ReferenceModel.dat"
unless a different file is specified using the "-i" option; see
below. This program produces by default an output file
"FIG_vev_vs_Q.dat", with data in 13 columns:
   
   1  Q (MSbar renormalization scale) 
   2  v(0)   (loopOrder = 0,   tree-level), 
   3  v(1)   (loopOrder = 1,   1-loop), 
   4  v(2)   (loopOrder = 2,   2-loop), 
   5  v(2.5) (loopOrder = 2.5, includes leading 3-loop), 
   6  v(3)   (loopOrder = 3,   full 3-loop), 
   7  v(3.5) (loopOrder = 3.5, full 3-loop + QCD 4-loop), 
   8  v(0)/vrun, 
   9  v(1)/vrun, 
   10  v(2)/vrun, 
   11  v(2.5)/vrun, 
   12  v(3)/vrun
   13  v(3.5)/vrun

where v(n) is the value of v obtained by minimizing the effective
potential at the scale Q in the given approximation, and vrun is
obtained by directly running the input value of v to the scale Q using
3-loop RG equations.

In the paper: arXiv:1907.02500
Standard Model parameters in the tadpole-free pure MSbar scheme
by Stephen P. Martin and David G. Robertson,
Figure 3.2a graphs columns 3, 4, 5, 6, and 7 as a function of column 1.
Figure 3.2b graphs columns 10, 11, 12, and 13 as a function of column 1.

The program takes two optional command line arguments:

   -i <input_filename>    Reads model data from <input_filename>; if 
                          not specified, data are read from the file
                          "ReferenceModel.dat".

   -o <output_filename>   Writes data to <output_filename>; if not
                          specified, the results appear in
                          "FIG_vev_vs_Q.dat".

The executable for this program is automatically created within the
smdr directory by make. You can also compile it separately as:

   gcc -o fig_vev_vs_Q fig_vev_vs_Q.c -L. -lsmdr -ltsil -l3vil -lm

assuming the necessary archive files libsmdr.a and libtsil.a and
lib3vil.a and the header file smdr.h and tsil.h and 3vil.h are all in
the current directory.

Run as:

   ./fig_vev_vs_Q

Running in the background is highly recommended, as the running time
is typically of order 2 hours, depending on your hardware.


*********************************************************************
IV. The SMDR Application Programmer Interface
*********************************************************************

In this section we define the signatures of all SMDR functions in the
user API, and describe their operation.  SMDR function names are all
preceded by SMDR_ to avoid collisions with user defined functions and
variables.

---------------------------------------------------------------------
             ==== RENORMALIZATION GROUP RUNNING ====
---------------------------------------------------------------------

For any of these functions, loopOrder should be 1, 2, 3, 4, 5.  In the 
case of loopOrder = 5, the known partial 5-loop order results are 
included.

---------------------------------------------------------------------
int SMDR_RGeval_SM (SMDR_REAL Q_final, int loopOrder);

Copies the input MSbar parameter global variables:

   Q_in, g3_in, g_in, ..., m2_in, v_in, Lambda_in,
   Delta_alpha_had_5_MZ_in

to the corresponding working global variables:

   Q, g3, g, ..., m2, v, Lambda, Delta_alpha_had_5_MZ.

Then, runs all of these variables (except the last, which is defined
at MZ) from the scale Q = Q_in to the scale Q_final. Uses beta
functions with loop order loopOrder. The results are stored in the
global variables:

   Q = Q_final, g3, g, gp, yt, yb, ..., m2, v, Lambda.

---------------------------------------------------------------------
int SMDR_RGrun_SM (SMDR_REAL Q_final, int loopOrder);

Runs the full Standard Model MSbar parameters in the non-decoupled
theory:

   Q, g3, g, gp, yt, yb, ..., m2, v, Lambda

from the current scale Q to the scale Q_final. Uses beta functions
with loop order loopOrder. The results are stored in the global
variables:

   Q = Q_final, g3, gp, g, yt, yb, ..., m2, v, Lambda.

---------------------------------------------------------------------
int SMDR_RGeval_QCDQED_53 (SMDR_REAL Qfinal,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder);

Evaluates, using RG running and decoupling, the 5-quark,
3-charged-lepton QCD+QED effective theory MSbar running parameter
global variables:

   alphaS, alpha, mb, mc, ms, md, mu, mtau, mmuon, melectron 

at the scale Qfinal, taking the decoupling scale for top, Higgs, Z, W
to be Q_thZW_dec. The inputs are the global variables corresponding to
the non-decoupled Standard Model input MSbar parameters:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., lambda_in, v_in

The results are stored as the global variables:

   Q_53 = Qfinal, alphaS_53, alpha_53, m<fermion>_53

with <fermion> = b, c, s, d, u, tau, muon, electron.  Typically,
Q_final is expected to be roughly between 4 GeV and MZ.

---------------------------------------------------------------------
int SMDR_RGeval_QCDQED_43 (SMDR_REAL Qfinal,
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder);

Evaluates, using RG running and decoupling, the 4-quark,
3-charged-lepton QCD+QED effective theory MSbar running parameters:

   alphaS, alpha, mc, ms, md, mu, mtau, mmuon, melectron

at the scale Qfinal, taking decoupling scales of Q_thZW_dec for top,
Higgs, Z, W, and taking Q_b_dec for the bottom quark decoupling scale.
The inputs are the global variables corresponding to the non-decoupled
Standard Model input MSbar parameters:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., lambda_in, v_in.

The results are stored as the global variables:

   Q_43 = Qfinal, alphaS_43, alpha_43, m<fermion>_43,

with <fermion> = c, s, d, u, tau, muon, electron.  Typically, Q_final
is expected to be roughly between 1.5 GeV and 5 GeV.

---------------------------------------------------------------------
int SMDR_RGeval_QCDQED_42 (SMDR_REAL Qfinal,
                           SMDR_REAL Q_tau_dec,
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder);

Evaluates, using RG running and decoupling, the 4-quark,
2-charged-lepton QCD+QED effective theory MSbar running parameters:

   alphaS, alpha, mc, ms, md, mu, mmuon, melectron 

at the scale Qfinal, taking decoupling scales of Q_thZW_dec for top,
Higgs, Z, W, and taking Q_b_dec for the bottom quark decoupling scale,
and taking Q_tau_dec for the tau lepton decoupling scale. The inputs
are the global variables corresponding to the non-decoupled Standard
Model input MSbar parameters:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., lambda_in, v_in.

The results are stored as the global variables:

   Q_42 = Qfinal, alphaS_42, alpha_42, m<fermion>_42

with <fermion> = c, s, d, u, muon, electron.  Typically, Q_final is
expected to be roughly between 1 GeV and 2 GeV.

---------------------------------------------------------------------
int SMDR_RGeval_QCDQED_32 (SMDR_REAL Qfinal,
                           SMDR_REAL Q_c_dec,
                           SMDR_REAL Q_tau_dec,
                           SMDR_REAL Q_b_dec,
                           SMDR_REAL Q_thZW_dec,
                           int loopOrder);

Evaluates, using RG running and decoupling, the 3-quark,
2-charged-lepton QCD+QED effective theory MSbar running parameters:

   alphaS, alpha, ms, mu, md, mmuon, melectron

at the scale Qfinal, taking decoupling scales of Q_thZW_dec for top,
Higgs, Z, W, and taking Q_b_dec for the bottom quark decoupling scale,
and Q_tau_dec for the tau lepton decoupling scale, and Q_c_dec for the
charm quark decoupling scale. The inputs are the non-decoupled
Standard Model input MSbar parameters:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., lambda_in, v_in.

The results are stored as the global variables:

   Q_32 = Qfinal, alphaS_32, alpha_32, m<fermion>_32, 

with <fermion> = s, d, u, muon, electron. Typically, Q_final is
expected to be less than about 1.5 GeV.

---------------------------------------------------------------------
int SMDR_RGrun_QCDQED_53 (SMDR_REAL Qfinal, int loopOrder);

Runs the MSbar parameters of the 5-quark, 3-charged lepton QCD+QED
effective theory:

   Q_53, alphaS_53, alpha_53, m<fermion>_53, 

with fermion = b, c, s, d, u, tau, muon, electron, from the current
scale Q_53 to the scale Q_final. Uses beta functions with loop order
loopOrder. The results are stored in the same global variables:
Q_53 = Qfinal, alphaS_53, alpha_53, m<fermion>_53.

---------------------------------------------------------------------
int SMDR_RGrun_QCDQED_43 (SMDR_REAL Qfinal, int loopOrder);

Runs the MSbar parameters of the 4-quark, 3-charged lepton QCD+QED
effective theory:

   Q_43, alphaS_43, alpha_43, m<fermion>_43, 

with <fermion> = c, s, d, u, tau, muon, electron, from the current
scale Q_43 to the scale Q_final. Uses beta functions with loop order
loopOrder. The results are stored in the same global variables:
Q_43 = Qfinal, alphaS_43, alpha_43, m<fermion>_43.

---------------------------------------------------------------------
int SMDR_RGrun_QCDQED_42 (SMDR_REAL Qfinal, int loopOrder);

Runs the MSbar parameters of the 4-quark, 2-charged lepton QCD+QED
effective theory:

   Q_42, alphaS_42, alpha_42, m<fermion>_42, 

with <fermion> = c, s, d, u, muon, electron from the current scale
Q_42 to the scale Q_final. Uses beta functions with loop order
loopOrder. The results are stored in the same global variables:
Q_42 = Qfinal, alphaS_42, alpha_42, m<fermion>_42.

---------------------------------------------------------------------
int SMDR_RGrun_QCDQED_32 (SMDR_REAL Qfinal, int loopOrder);

Runs the MSbar parameters of the 3-quark, 2-charged lepton QCD+QED
effective theory:

   Q_32, alphaS_32, alpha_32, m<fermion>_32,

with <fermion> = s, d, u, muon, electron, from the current scale Q_32
to the scale Q_final. Uses beta functions with loop order
loopOrder. The results are stored in the same global variables:

   Q_32 = Qfinal, alphaS_32, alpha_32, m<fermion>_32.

---------------------------------------------------------------------
int SMDR_RGrun_QCDQED (SMDR_REAL Q_init, SMDR_REAL Q_final,
                       int loopOrder, int nu, int nd, int ne,
                       SMDR_REAL alphaS_init, SMDR_REAL alpha_init,
                       SMDR_REAL *alphaS_final, SMDR_REAL *alpha_final,
                       SMDR_REAL *cu, SMDR_REAL *cd, SMDR_REAL *ce);

This function is a tool used by the preceding functions. 

Runs the MSbar parameters of the QCD+QED theory with nu up-type
quarks, nd down-type quarks, and ne charged leptons, from the scale
Q_init to the scale Q_final. The values of the QCD and QED couplings
at Q_init are alphaS_init and alpha_init. The output results at
Q_final are: alphaS_final, alpha_final, and the ratios
cu=mu_final/mu_init, cd=md_final/md_init, ce=me_final/me_init.


---------------------------------------------------------------------
            ==== DECOUPLING AND MATCHING ROUTINES ====
---------------------------------------------------------------------

void SMDR_Decouple_thZW (int loopOrder);

Simultaneously decouples the top quark and the Higgs, Z, and W bosons,
at the current MSbar renormalization scale Q. The inputs are the
global variables corresponding to the full Standard Model:

   Q, g3, g, gp, yt, yb, yc, ys, yd, yu, ytau, ymuon, yelectron,
   lambda, v

The results are stored in the global variables: 

   alphaS_53, alpha_53, and m<fermion>_53, 

where <fermion> = b, c, s, u, d, tau, muon, electron.  These are the
MSbar QCD and EM couplings and the running masses in the QCD+QED
effective theory that contains only the 5 quarks b,c,s,u,d and the
charged leptons tau,mu,e, and the gluons and the photons.  At present,
the argument loopOrder is ignored, and all known effects are included.

---------------------------------------------------------------------
void SMDR_Eval_Gauge (SMDR_REAL Mtpole,
     		      SMDR_REAL Mhpole,
		      SMDR_REAL MWpole);

Evaluates:
  1) alphaS(MZ) and alphaS(MZ) and s^2(thetaW) in the MSbar scheme of
     the non-decoupled Standard Model theory, putting the results in
     the global variables: SMDR_alphaS_MZ, SMDR_alpha_MZ, SMDR_s2W_MZ.
  2) The Sommerfeld fine structure constant, using the results of
     1411.7040 Degrassi, Gambino, and Giardino. The result is stored
     in the global variable SMDR_alpha.
  3) alpha(MZ) and s^2(thetaW) in PDG MSbar scheme with only top
     decoupled but W not decoupled. The results are stored in the
     global variables: SMDR_alpha_MZ_PDG, SMDR_s2W_MZ_PDG.

Uses as inputs the global variable Standard Model input values

   Q_in, g3_in, g_in, gp_in, yt_in, ...

as well as the arguments Mtpole, Mhpole, MWpole.

---------------------------------------------------------------------

For the following four functions, loopOrder = 0, 1, 2, 3, 4.

EM contributions beyond 2 loops are not included.

---------------------------------------------------------------------
void SMDR_Decouple_bottom (int loopOrder);

Decouples the bottom quark at the scale given by the global variable
Q_53 = the MSbar renormalization scale of the 5-quark,
3-charged-lepton QCD+QED effective field theory. The inputs are the
global variables:

   Q_53, alphaS_53, alpha_53, m<Fermion>_53

with <Fermion> = b, c, s, d, u, tau, muon, electron.  The output MSbar
parameters of the 4-quark, 3-charged-lepton QCD+QED effective theory
are put in the global variables:

   Q_43 = Q_53, alphaS_43, alpha_43, m<fermion>_43,

with <fermion> = c, s, d, u, tau, muon, electron,

---------------------------------------------------------------------
void SMDR_Decouple_tau (int loopOrder);

Decouples the tau lepton at the scale given by the global variable
Q_43 = the MSbar renormalization scale of the 4-quark,
3-charged-lepton QCD+QED effective field theory. The input are the
global variables:

   Q_43, alphaS_43, alpha_43, m<Fermion>_43

with <Fermion> = c, s, d, u, tau, muon, electron.  The output MSbar
parameters of the 4-quark, 2-charged-lepton QCD+QED effective theory
are put in the global variables:

   Q_42 = Q_43, alphaS_42, alpha_42, m<fermion>_42,

with <fermion> = c, s, d, u, muon, electron,

---------------------------------------------------------------------
void SMDR_Decouple_charm (int loopOrder);

Decouples the charm quark at the scale given by the global variable
Q_42 = the MSbar renormalization scale of the 4-quark,
2-charged-lepton QCD+QED effective field theory. The inputs are the
global variables:

   Q_42, alphaS_42, alpha_42, m<Fermion>_42

with <Fermion> = c, s, d, u, muon, electron.  The output MSbar
parameters of the 3-quark, 2-charged-lepton QCD+QED effective theory
are put in the global variables:

   Q_32 = Q_42, alphaS_32, alpha_32, m<fermion>_32,

with <fermion> = s, d, u, muon, electron,

---------------------------------------------------------------------
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
                       SMDR_REAL *ze);

This function is a tool used by the preceding functions.

The input arguments are:

  Fermion_type = "u", "d", or "e" = type of heavy fermion being
                 decoupled
  Q_match      = MSbar decoupling scale.
  m_Fermion    = MSbar mass of heavy fermion being decoupled, in full
               	 theory, at Q_match
  alphaS_hi    = strong coupling in full theory, at scale Q_match.
  alpha_hi     = EM coupling in full theory, at scale Q_match
  nqlight      = number of light quarks, not including the fermion
                 being decoupled. Note that the number of light
                 leptons doesn't matter in the approximation being
                 used. (It would matter if 3-loop EM contributions
                 were included.)

The outputs are:

  alphaS_lo = strong coupling in decoupled theory, at scale Q_match.
  alpha_lo  = EM coupling in decoupled theory, at scale Q_match.
  zu        = ratio of (MSbar light up-type mass in decoupled theory)/
              (MSbar light up-type mass in full theory).
  zd        = ratio of (MSbar light down-type mass in decoupled
              theory)/ (MSbar light down-type mass in full theory).
  ze        = ratio of (MSbar light charged lepton mass in decoupled
              theory)/ (MSbar light charged lepton mass in full
              theory).

Small mass corrections are not included here, but are added in the
functions SMDR_Decouple_b and SMDR_Decouple_c.


---------------------------------------------------------------------
                  ==== EFFECTIVE POTENTIAL ====
---------------------------------------------------------------------

For these functions, loopOrder may take the following values:

   0     tree-level
   1     1-loop
   2     2-loop
   2.5   2-loop plus 3-loop contributions in the large
         g3, yt limit
   3     full 3-loop
   3.5   full 3-loop plus 4-loop at leading order in QCD

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_Veffmin (SMDR_REAL Q_eval, float loopOrder);

Returns the minimum value of the effective potential (in Landau
gauge, in the MSbar scheme, with Goldstone boson resummation).

If the argument Q_eval is positive, then the inputs are obtained 
by first RG running the MSbar parameter global variables:

   Q_in, v_in, lambda_in, g3_in, g_in, gp_in, yt_in, yb_in, ytau_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead
taken to be the current values of the MSbar parameter global
working variables:

   Q, v, lambda, g3, g, gp, yt, yb, ytau

The function adjusts the global working variable m2 to assure that
the resummed effective potential is minimized, and returns the
corresponding minimum value of the effective potential.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_m2 (SMDR_REAL Q_eval, float loopOrder);

Returns the value of m2 that minimizes the effective potential (in
Landau gauge, in the MSbar scheme, with Goldstone boson resummation).

If the argument Q_eval is positive, then the inputs are obtained
by first RG running the MSbar parameter global variables:

   Q_in, v_in, lambda_in, g3_in, g_in, gp_in, yt_in, yb_in, ytau_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead
taken to be the current values of the MSbar parameter global
working variables:

   Q, v, lambda, g3, g, gp, yt, yb, ytau

which are assumed to be fixed. This function does not automatically
change the global variable m2.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_vev (SMDR_REAL Q_eval, float loopOrder);

Returns the value of the VEV v at the minimum of the effective
potential (in Landau gauge, in the MSbar scheme, with Goldstone
boson resummation, at loop order loopOrder)

If the argument Q_eval is positive, then the inputs are obtained
by first RG running the MSbar parameter global variables:

   Q_in, m2_in, lambda_in, g3_in, g_in, gp_in, yt_in, yb_in, ytau_in,

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead
taken to be the current values of the MSbar parameter global
working variables:

   Q, m2, lambda, g3, g, gp, yt, yb, ytau,

which are assumed to be fixed. This function does not automatically
change the global variable v.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_vevDelta (SMDR_REAL Q_eval, float loopOrder);

Returns the value of the quantity

   Delta = -G = (vtree^2 - v^2)*lambda

defined in section V.3 of 1709.02397.

If the argument Q_eval is positive, then the inputs are obtained
by first RG running the MSbar parameter global variables:

   Q_in, v_in, lambda_in, g3_in, g_in, gp_in, yt_in, yb_in, ytau_in,

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead
taken to be the current values of the MSbar parameter global working variables:

   Q, v, lambda, g3, g, gp, yt, yb, ytau,

which are assumed to be fixed.


---------------------------------------------------------------------
                     ==== HIGGS SECTOR ====
---------------------------------------------------------------------

For these functions, loopOrder may take the following values:

   0    tree level
   1    1-loop
   1.5  1-loop plus 2-loop QCD corrections
   2    full 2-loop
   2.3  full 2-loop plus leading 3-loop QCD
   2.5  full 2-loop plus leading 3-loop QCD and nonQCD terms.

The computations use equations (2.43), (2.46)-(2.48), and
(3.2)-(3.4) of [MR14].

---------------------------------------------------------------------
void SMDR_Eval_Mh_pole (SMDR_REAL Q_eval,
                        float loopOrder, 
                        SMDR_REAL *Mhpoleresult,
                        SMDR_REAL *Gammahpoleresult);

Computes the Higgs pole mass. The results M_h and Gamma_h for the
complex pole squared mass

   M_h^2 - i Gamma_h M_h 

are returned as Mhpoleresult and Gammahpoleresult.

If the argument Q_eval is positive, then the inputs are obtained by
first RG running the MSbar parameter global variables:

   Q_in, g3_in, g_in, gp_in, lambda_in, yt_in, ..., v_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead given
directly by the current values of the MSbar parameter global
variables:

   Q, g3, g, gp, lambda, yt, ..., v.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_lambda (SMDR_REAL mhpole,
                            SMDR_REAL Q_eval, 
                            float loopOrder);

Returns the Higgs self-coupling lambda, given the real part of the
Higgs pole mass and the other MSbar Standard Model parameters,
determined as follows:

If the argument Q_eval is positive, then the inputs are obtained by
first RG running the MSbar parameter global variables:

   Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead taken to be 
the current values of the MSbar parameter global variables:

   Q, g3, g, gp, yt, ..., v.


---------------------------------------------------------------------
                      ==== W, Z BOSONS ====
---------------------------------------------------------------------

For these functions, loopOrder may take the following values:

   0    tree level
   1    1-loop
   1.5  1-loop plus 2-loop QCD corrections
   2    full 2-loop

---------------------------------------------------------------------
void SMDR_Eval_MW_pole (SMDR_REAL Q_eval,
                        float loopOrder,
                        SMDR_REAL *MWpoleresult,
                        SMDR_REAL *GammaWpoleresult,
                        SMDR_REAL *MWBreitWignerresult,
                        SMDR_REAL *GammaWBreitWignerresult);

Computes the complex pole and Breit-Wigner squared masses of the W
boson, using the calculation from 1503.03782.

If the argument Q_eval is positive, then the inputs are obtained by
first RG running the MSbar parameter global variables:

   Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead given
by the current values of the MSbar parameter global variables:

   Q, g3, g, gp, yt, ..., v.

The results for the complex pole squared mass 

   M_W^2 - i Gamma_W M_W

are returned in MWpoleresult, GammaWpoleresult.

The results for the complex Breit-Wigner squared mass are also
returned as MWBreitWignerresult, GammaWBreitWignerresult.

---------------------------------------------------------------------
void SMDR_Eval_MZ_pole (SMDR_REAL Q_eval,
                        float loopOrder,
                        SMDR_REAL *MZpoleresult,
                        SMDR_REAL *GammaZpoleresult,
                        SMDR_REAL *MZBreitWignerresult,
                        SMDR_REAL *GammaZBreitWignerresult);

Computes the complex pole and Breit-Wigner squared masses of the Z
boson, using the calculation from 1505.04833.

If the argument Q_eval is positive, then the inputs are obtained by
first RG running the MSbar parameter global variables:

   Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead given
by the current values of the MSbar parameter global variables:

   Q, g3, g, gp, yt, ..., v.

The results for the complex pole squared mass

   M_Z^2 - i Gamma_Z M_Z

are returned as MZpoleresult, GammaZpoleresult.

The results for the complex Breit-Wigner squared mass are also
returned as MZBreitWignerresult, GammaZBreitWignerresult.


---------------------------------------------------------------------
                       ==== GFERMI ====
---------------------------------------------------------------------

SMDR_REAL SMDR_Eval_GFermi (SMDR_REAL Q_eval, float loopOrder);

Returns GFermi at up to two loops. Based on Deltartilde given in
section 3 of arXiv:1907.02500 and the ancillary file Deltartilde.txt.

See also Deltarbar in
   Kniehl and Veretin 1401.1844 eqs (37)-(40) and
   Kniehl, Pikelner, Veretin 1503.02138 eqs (60)-(62). Appendix A.2
   Kniehl, Pikelner, Veretin 1601.08143, file
   mr-1.3.2/mr/dr/drbar20.cpp in the computer code "mr".
 
The argument loopOrder has allowed values:
    0    tree-level
    1    1-loop
    1.3  1-loop plus 2-loop leading order in QCD
    1.5  1-loop plus 2-loop leading order in QCD and yt
    2     Full 2-loop result

If the argument Q_eval is positive, then the inputs are obtained by
first RG running the MSbar parameters:

    Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead given
by the current values of the MSbar parameters:

    Q, g3, g, gp, yt, ..., v.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_GFermi_DGG (SMDR_REAL Mtpole,
                                SMDR_REAL Mhpole,
                                SMDR_REAL MWpole);

Returns the result for the Fermi decay constant G_Fermi as found by
1411.7040 Degrassi, Gambino, and Giardino. Uses as inputs the global
variable Standard Model input values

   Q_in, g3_in, g_in, gp_in, yt_in, ...

as well as the arguments Mtpole, Mhpole, MWpole.  This is an
alternative to SMDR_Eval_GFermi().


---------------------------------------------------------------------
                      ==== TOP QUARK ====
---------------------------------------------------------------------

For these functions,

   method = 0 (expand around tree-level) or
            1 (expand around pole. Recommended; more accurate but
               slower).

   QCDLoopOrder = 0, 1, 2, 3, or 4.

   otherLoopOrder = 0 (tree-level) or
                    1 (1-loop) or
                    1.5 (2-loop mixed QCD/EW) or
                    2 (2-loop full).

The sources for the relevant results are:
   4-loop pure QCD from 1502.01030, updated and perfected in
   1606.06754.
   2-loop non-pure-QCD from 1604.01134.

---------------------------------------------------------------------
void SMDR_Eval_Mt_pole (SMDR_REAL Q_eval,
                        int method, 
                        int QCDloopOrder, 
                        float otherLoopOrder,
                        SMDR_REAL *Mtpoleresult, 
                        SMDR_REAL *Gammatpoleresult);

Computes the complex pole mass of the top quark. The results for the
complex pole squared mass M_t^2 - i Gamma_t M_t are returned as
Mtpoleresult, Gammatpoleresult.

If the argument Q_eval is positive, then the inputs are obtained by
first RG running the MSbar parameter global variables:

   Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead given
by the current values of the MSbar parameter global variables:

   Q, g3, g, gp, yt, ..., v.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_yt (SMDR_REAL MTpoletarget,
                        SMDR_REAL Q_eval,
                        int method,
                        int QCDLoopOrder,
                        float otherLoopOrder);

Returns the top-quark Yukawa coupling, given the real part of the pole
mass of the top quark, specified as MTpoletarget. The other arguments
are as described for the previous function SMDR_Eval_Mt_pole().

If the argument Q_eval is positive, then the inputs are obtained
by first RG running the MSbar parameter global variables:

   Q_in, lambda_in, g3_in, g_in, gp_in, yb_in, ..., v_in

to the scale Q_eval, using SMDR_RGeval_SM().

If the argument Q_eval is negative, then the inputs are instead taken 
to be the current values of the MSbar parameter global variables:

   Q, lambda, g3, g, gp, yb, ..., v.


---------------------------------------------------------------------
                ==== LIGHT QUARKS AND LEPTONS ====
---------------------------------------------------------------------

void SMDR_Eval_Light_Masses ();

Has the same effect as calling all of the following individual
functions:

   SMDR_Eval_mbmb(); 
   SMDR_Eval_mcmc(); 
   SMDR_Eval_mquarks_2GeV();
   SMDR_Eval_Mtau_pole(); 
   SMDR_Eval_Mmuon_pole(); 
   SMDR_Eval_Melectron_pole();

However, it is slightly more efficient than calling them individually, 
because it doesn't unnecessarily redo the RG running each time.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_mbmb (SMDR_REAL Q_dec_thZW, int loopOrder);

Returns the bottom-quark MSbar running mass evaluated at itself,
mb(mb), in the 5-quark, 3-charged-lepton QCD+QED effective theory.
Uses beta functions with loop order loopOrder = 1, 2, 3, 4, 5.  The
inputs are the global variables corresponding to the MSbar parameters
of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, ..., v_in.

Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_mcmc (SMDR_REAL Q_dec_tau,
                          SMDR_REAL Q_dec_bottom,
                          SMDR_REAL Q_dec_thZW, int loopOrder);

Returns the charm-quark MSbar running mass evaluated at itself,
mc(mc), in the 4-quark, 2-charged-lepton QCD+QED effective theory.
Uses beta functions with loop order loopOrder = 1, 2, 3, 4, 5.  The
inputs are the global variables corresponding to the MSbar parameters
of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
decoupling of the bottom and tau occur at Q_dec_bottom, Q_dec_tau.

---------------------------------------------------------------------
void SMDR_Eval_mquarks_2GeV (SMDR_REAL Q_dec_bottom,
                             SMDR_REAL Q_dec_thZW, 
                             int loopOrder,
                             SMDR_REAL *ms,
                             SMDR_REAL *mu,
                             SMDR_REAL *md);

Evaluates the strange, down, and up-quark MSbar running masss
evaluated at Q=2 GeV, in the 4-quark, 3-charged-lepton QCD+QED
effective theory.  Uses beta functions with loop order loopOrder = 1,
2, 3, 4, 5.  The inputs are the global variables corresponding to the
MSbar parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
decoupling of the bottom at Q_dec_bottom.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_Mb_pole (SMDR_REAL Q_eval, 
                             SMDR_REAL Q_dec_thZW, 
                             int loopOrder);

Returns the pole mass of the bottom quark, computed in the 5-quark,
3-charged-lepton QCD+QED effective field theory at the scale Q_eval,
in the approximation of QCD loopOrder = 0, 1, 2, 3, or 4.  (The RG
running is done including all known effects up to 5-loop order.)  The
starting inputs are the global variables corresponding to the MSbar
parameters of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW.

NOTE: This just doesn't converge fast enough to be useful!

The Review of Particle Properties uses the 2-loop approximation, which
is horrible. See 1606.06754 for discussion of the atrocious
convergence properties of the expansion. Therefore, the bottom (and
charm) pole masses are kind of useless as observables and are
deprecated.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_Mtau_pole (SMDR_REAL Q_eval,
                               SMDR_REAL Q_dec_bottom,
                               SMDR_REAL Q_dec_thZW,
                               int loopOrder);

Returns the pole mass of the tau lepton, computed in the 4-quark,
3-charged-lepton QCD+QED effective field theory at the scale Q_eval,
in the approximation of loopOrder = 0, 1, 2, 3.  (The RG running is
done including all known effects up to 5-loop order.)  The starting
inputs are the global variables corresponding to the MSbar parameters
of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
decoupling of the bottom quark occurs at the scale Q_dec_bottom.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_Mmuon_pole (SMDR_REAL Q_eval,
                                SMDR_REAL Q_dec_charm,
                                SMDR_REAL Q_dec_tau,
                                SMDR_REAL Q_dec_bottom,
                                SMDR_REAL Q_dec_thZW,
                                int loopOrder);

Returns the pole mass of the muon, computed in the 3-quark,
2-charged-lepton QCD+QED effective field theory at the scale Q_eval,
in the approximation of loopOrder = 0, 1, or 2.  (The RG running is
done including all known effects up to 5-loop order.)  The starting
inputs are the global variables corresponding to the MSbar parameters
of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
decoupling of the bottom, tau, and charm occur at the scales
Q_dec_bottom, Q_dec_tau, and Q_dec_charm, respectively.

---------------------------------------------------------------------
SMDR_REAL SMDR_Eval_Melectron_pole (SMDR_REAL Q_eval,
                                    SMDR_REAL Q_dec_charm,
                                    SMDR_REAL Q_dec_tau,
                                    SMDR_REAL Q_dec_bottom,
                                    SMDR_REAL Q_dec_thZW,
                                    int loopOrder);

Returns the pole mass of the electron, computed in the 3-quark,
2-charged-lepton QCD+QED effective field theory at the scale Q_eval,
in the approximation of loopOrder = 0, 1, or 2.  (The RG running is
done including all known effects up to 5-loop order.)  The starting
inputs are the global variables corresponding to the MSbar parameters
of the full Standard Model:

   Q_in, g3_in, g_in, gp_in, yt_in, yb_in, yc_in, ..., v_in.

Decoupling of t, h, Z, W occurs at the scale Q_dec_thZW, and
decoupling of the bottom, tau, and charm occur at the scales
Q_dec_bottom, Q_dec_tau, and Q_dec_charm, respectively.


---------------------------------------------------------------------
                      ==== UTILITIES ====
--------------------------------------------------------------------- 

int SMDR_Fit_Inputs (SMDR_REAL Q_target,
                     SMDR_REAL alphaS_target,
                     SMDR_REAL alpha_target,
                     SMDR_REAL GFermi_target,
                     SMDR_REAL MZ_target,
                     SMDR_REAL Mh_target,
                     SMDR_REAL Mt_target,
                     SMDR_REAL mbmb_target,
                     SMDR_REAL mcmc_target,
                     SMDR_REAL ms_2GeV_target,
                     SMDR_REAL md_2GeV_target,
                     SMDR_REAL mu_2GeV_target,
                     SMDR_REAL Mtau_target,    
                     SMDR_REAL Mmuon_target,
                     SMDR_REAL Melectron_target,
                     SMDR_REAL Delta_alpha_had,
                     SMDR_REAL error_target);

Finds the (non-decoupled) Standard Model MSbar parameters at the scale
Q_target. The results are stored in the global variables:

   Q_in = Q_target, g3_in, g_in, gp_in, v_in, lambda_in, y<fermion>_in

where fermion = t, b, c, s, d, u, tau, muon, electron.

The inputs are the arguments:

   alphaS_target (in the 5-quark, 3-charged-lepton QCD+QED theory at
                  Q=MZ),
   alpha_target (in the 5-quark, 3-charged-lepton theory at Q=MZ),
   MZ_target (the Z-boson Breit-Wigner mass),
   GFermi_target (the Fermi constant),
   Mh_target (the Higgs boson pole mass),
   Mt_target (the top-quark pole mass),
   mbmb_target (the MSbar bottom mass evaluated at itself in the
                5-quark, 3-charged-lepton QCD+QED theory),
   mcmc_target (the MSbar charm mass evaluated at itself in the
                4-quark, 2-charged-lepton QCD+QED theory),
   ms_2GeV_target (the MSbar strange quark mass at Q=2 GeV, in the
                   3-quark, 2-charged-lepton QCD+QED theory),
   md_2GeV_target (the MSbar down quark mass at Q=2 GeV, in the
                   3-quark, 2-charged-lepton QCD+QED theory),
   mu_2GeV_target (the MSbar up quark mass at Q=2 GeV, in the 3-quark,
                   2-charged-lepton QCD+QED theory),
   Mtau_target (the tau lepton pole mass),
   Mmuon_target (the muon pole mass),
   Melectron_target (the electron pole mass),
   Delta_alpha_had (the non-perturbative hadronic contribution to the
                    fine-structure constant at MZ)

The calculation proceeds by iteration, until the largest of the
fractional errors in the target quantities is less than
error_target. About 5 or fewer iterations are typically expected for
error_target = 10^-8.

---------------------------------------------------------------------
void SMDR_Load_Inputs (void);

Loads the Standard Model MSbar parameters from the "static" global input 
variables X_in to the corresponding global working variables X, where

   X = Q, g3, g, gp, yt, yb, yc, ys, yu, yd, ytau, ymu, ye, lambda,
       m2, v, Lambda, and Delta_alpha_had_5_MZ.

This function should be called before calling other functions that need 
the working variables X. It should also be called whenever a preceding 
function may have altered one or more of the working variables X, if it 
is necessary to reset them to the original input variables X_in. In 
practice, many other functions start by calling SMDR_Load_Inputs() 
anyway, in which case it is not necessary for the user application to do 
so. This includes all functions whose names contain "RGeval" or "Eval". 
The exception to this rule is functions with "Eval" in the special case 
that they are called with the first argument Q_eval negative; this 
specifies that the working variables at the current value of Q will be 
used. In those cases, one should have either previously called 
SMDR_Load_Inputs(), or else SMDR_RGeval_SM() to set the working 
variables at the desired current scale Q.

Another common situation where the function SMDR_Load_Inputs() is 
necessary in user applications is when one wants to write or display the 
*input* values of the MSbar variables but using one of the functions 
that are designed to write or display the working variables, namely: 
SMDR_Write_MSbar_Parameters() and SMDR_Write_v() and SMDR_Write_m2() and 
SMDR_Write_Lambda() or the corresponding functions 
SMDR_Display_MSbar_Parameters() and SMDR_Display_v() and 
SMDR_Display_m2() and SMDR_Display_Lambda().

---------------------------------------------------------------------
void SMDR_Save_Inputs (void);

Loads the Standard Model MSbar parameters from the current global 
working variables X to the corresponding global input variables X_in, where

   X = Q, g3, g, gp, yt, yb, yc, ys, yu, yd, ytau, ymu, ye, lambda,
       m2, v, Lambda, and Delta_alpha_had_5_MZ.

This function is the inverse of SMDR_Load_Inputs().

---------------------------------------------------------------------
void SMDR_Update (void);

Updates various useful parameter combinations.

In general, this function is called automatically when needed by the
API functions, so users need not worry about it or use it under normal
circumstances.

---------------------------------------------------------------------
void SMDR_Check_Ranges (void);
void SMDR_Check_Q_Range (SMDR_REAL Q_range_lo, SMDR_REAL Q_range_hi);
void SMDR_Check_VEV_Range (SMDR_REAL vev_range_lo, SMDR_REAL vev_range_hi);
void SMDR_Check_m2_Range (SMDR_REAL m2_range_lo, SMDR_REAL m2_range_hi);
void SMDR_Check_lambda_Range (SMDR_REAL k_range_lo, SMDR_REAL k_range_hi);
void SMDR_Check_yt_Range (SMDR_REAL yt_range_lo, SMDR_REAL yt_range_hi);
void SMDR_Check_yb_Range (SMDR_REAL yb_range_hi);
void SMDR_Check_ytau_Range (SMDR_REAL ytau_range_hi);
void SMDR_Check_g3_Range (SMDR_REAL g3_range_lo, SMDR_REAL g3_range_hi);
void SMDR_Check_g_Range (SMDR_REAL g_range_lo, SMDR_REAL g_range_hi);
void SMDR_Check_gp_Range (SMDR_REAL gp_range_lo, SMDR_REAL gp_range_hi);
void SMDR_Check_Mhpole_Range (SMDR_REAL Mhpole, 
                              SMDR_REAL Mh_range_lo,
                              SMDR_REAL Mh_range_hi);
void SMDR_Check_Mtpole_Range (SMDR_REAL Mtpole, 
                              SMDR_REAL Mt_range_lo,
                              SMDR_REAL Mt_range_hi);

These functions performs sanity checks for the Standard Model MSbar
parameters g3, g, gp, lambda, yb, ytau, v, Q.  If out of range, a
warning message is printed (but execution is not halted).

---------------------------------------------------------------------
void SMDR_Eval_thZW_pole (void);

Computes the complex pole masses of the top quark, Higgs boson, Z
boson, and W bozon, and stores the results in the global variables:

   Mt_pole, Gammat_pole
   Mh_pole, Gammah_pole
   MZ_pole, GammaZ_pole, MZ_BreitWigner, GammaZ_BreitWigner
   MW_pole, GammaW_pole, MW_BreitWigner, GammaW_BreitWigner

Works by calling the functions SMDR_Eval_Mt_pole(),
SMDR_Eval_Mh_pole(), SMDR_Eval_MZ_pole(), and SMDR_Eval_MW_pole(). In
each case, the maximum available loopOrder is used.

The input parameters are the Standard Model MSbar global variables
inputs:

   Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in.

---------------------------------------------------------------------
void SMDR_Eval_QCDQED_at_MZ (SMDR_REAL Q_MZ, 
                             SMDR_REAL Q_dec_thZW,
                             int loopOrder);

Computes the MSbar couplings and masses in the 5-quark,
3-charged-lepton QCD+QED effective theory, at the scale Q_MZ, taking
t, h, Z, W to be decoupled at Q_dec_thZW. The RG running is evaluated
at loopOrder = 1, 2, 3, 4, or 5.

The input parameters are the Standard Model MSbar global variables
inputs:

   Q_in, g3_in, g_in, gp_in, yt_in, ..., v_in.

The results are put in the global variables:

   alphaS_5_MZ, alpha_5_MZ, m<fermion>_53, with

fermion = b, c, s, d, u, tau, muon, electron.

---------------------------------------------------------------------
void SMDR_Make_Accelcoeffs ();

Evaluates coefficients that allow for accelerated convergence of the
iteration process in SMDR_Fit_Inputs(), and stores them in the file
src/accelcoeffs.h. This function is called by the command line program

   ./make_coeffs

It should only be necessary to run this when a new version of SMDR is
being prepared involving a different version of the reference model in
ReferenceModel.dat.

---------------------------------------------------------------------
void SMDR_Error (char *funcName, char *errMsg, int status);

Prints a message of the form:

       ERROR (funcName): errMsg

and exits with status code status. A stack trace is also printed for
diagnostic purposes. SMDR error messages appear on stdout by default,
so they are seen even if warnings have been disabled or redirected.

---------------------------------------------------------------------
void SMDR_Warn (char *, char *);

Prints a message of the form:

       WARN (funcName): errMsg

and execution continues. SMDR warnings appear on stderr so they can be
redirected by the shell or otherwise suppressed, e.g., by

       fclose (stderr);

---------------------------------------------------------------------
SMDR_REAL SMDR_SGNSQRT (SMDR_REAL x);

Returns sgn(x) sqrt(|x|).

---------------------------------------------------------------------
void   SMDR_Start_Timer (void);
double SMDR_Timer (void);

Functions related to timing. See sample programs for usage.

---------------------------------------------------------------------
int SMDR_Process_Arguments (int argc,
                            char *argv[],
                            int nargs, 
                            char *arglist[],
                            char *argtype[],
                            void *argvar[]);

Simple utility to process command line arguments.  The arguments and
their types are defined in the calling program, and passed in via the
variables:

nargs
-----

The total number of defined arguments, equal to the sum of the number
of required arguments and the number of possible optional arguments.

arglist 
-------
A list of strings (of dimension nargs) specifying the
arguments. *Required* arguments must come first, and are given in a
specified order that cannot be modified. A required argument is
indicated with "req". After all required args are given, the arglist
entries should define a set of flags that begin with '-'. Each of
these specifies an *optional* argument. As an example, if there are
three required arguments and two optional ones, then arglist might
look like this:

   char arglist[] = {"req","req","req","-a","-i"};

In this case the two optional argw would be specified on the command
line as

   [-a <arg1> [-i <arg2>]]

Optional arguments can be given in any order on the commad line, but
the flags specified here must match up with the correspondng elements
in the following arrays that define them.
   
argtype
-------
An array of strings (of dimension nargs) defining the types of the
arguments, in order. Allowed type strings are currently "real", "int",
"string", and "toggle". As an example,

   char argtype[] = {"string","int","real","real","string"};

Type "toggle" is used for an argument that acts as a toggle, i.e.,
without an associated value.

argvar
------
This is an array of pointers-to-void (of dimension nargs), to data
objects in the calling program to which the read values should be
assigned. As an example,

   void argvar[] = {inFile, &i, &x, &Q, outFile};

would assign the first argument in the variable inFile, the second in
the variable i, and so on. Note that no ampersand is required for
"string"-type variables.

For arguments of type "toggle", the assigned variable should normally
be of type int.  The effect of including the toggle will be to set
this int to 1 (YES).

As another example, for three optional arguments corresponding to an
input filename, and output filename, and a real error tolerance, one
could specify:

   int nargs = 3;
   char *arglist[] = {"-e","-i","-o"};
   char *argtype[] = {"real","string","string"};
   void *argvar[] = {&ERROR_TOLERANCE, inputFile, outputFile};

Users should be sure to specify default values for all optional
arguments in the calling program.

Some simple tests for errors are included, but these are probably not
bulletproof. Caveat emptor!

---------------------------------------------------------------------
                     ==== INPUT/OUTPUT ====
--------------------------------------------------------------------- 

int SMDR_Read_MSbar_Inputs (char *fileName);

Reads the values of the following set of MSbar inputs from a file:
    SMDR_Q_in
    SMDR_g3_in
    SMDR_gp_in
    SMDR_g_in
    SMDR_yt_in
    SMDR_yb_in
    SMDR_yc_in
    SMDR_ys_in
    SMDR_yu_in
    SMDR_yd_in
    SMDR_ytau_in
    SMDR_ymu_in
    SMDR_ye_in
    SMDR_lambda_in

Note that the quantities:

  SMDR_v_in, SMDR_m2_in, SMDR_Lambda_in, and SMDR_Delta_alpha_had_5_MZ_in

are not read in by this function. This is because they are not always
needed by all applications. They should be read in separately using
SMDR_Read_Value() or SMDR_Read_Values(), described below.

The file format should be as in ReferenceModel.dat, although the lines
can be in any order, the trailing semicolons can be included or not,
and any line starting with a '#' is treated as a comment. Blank lines
are ignored.

If any of the listed variables are not found, they are listed and a
fatal error results. 

---------------------------------------------------------------------
int SMDR_Read_OS_Inputs (char *fileName);

Reads the values of the following on-shell observables from a file:
    SMDR_Delta_alpha_had_5_MZ_in
    SMDR_Mt_pole
    SMDR_Mh_pole
    SMDR_MZ_BreitWigner
    SMDR_GFermi
    SMDR_alpha
    SMDR_alphaS_5_MZ
    SMDR_mbmb
    SMDR_mcmc
    SMDR_ms_2GeV
    SMDR_md_2GeV
    SMDR_mu_2GeV
    SMDR_Mtau_pole
    SMDR_Mmuon_pole
    SMDR_Melectron_pole

The file format should be as in ReferenceModel.dat, although the lines
can be in any order, the trailing semicolons can be included or not,
and any line starting with a '#' is treated as a comment. Blank lines
are ignored.

If any of the listed variables are not found, they are listed and a
fatal error results.

--------------------------------------------------------------------- 
int SMDR_Get_Value (FILE *file, char *varName);

Searches the file stream file for the variable varName, and sets the
corresponding SMDR variable to the value found. If the value is not
found, a warning is issued.

The list of variables that may be searched can be found in

   src/smdr_ioparameters.h

--------------------------------------------------------------------- 
int SMDR_Read_Value (char *fileName, char *varName);

Reads from a file a single input value with the name specified by the
character string varName. Here, varName can be "SMDR_Q_in",
"SMDR_g3_in", etc.

If the specified variable is not found in fileName, a fatal error
results.

The list of variables that may be searched can be found in

   src/smdr_ioparameters.h

--------------------------------------------------------------------- 
int SMDR_Read_Values (char *fileName, int n, char *listOfNames[]);

Reads from a file a set of n input values with the name specified by
the array of strings listOfNames. Here, elements of the array
listOfNames can be "SMDR_Q_in", "SMDR_g3_in", etc.

If any of the specified variables is not found in fileName, a fatal error
results.

The list of variables that may be searched can be found in

   src/smdr_ioparameters.h

--------------------------------------------------------------------- 
int SMDR_Set_Values_Interactively (int, char *listOfNames[]);

A simple function that takes an integer number of variables and a list 
of their names, and steps through the list, printing the current value 
of each to stdout and prompting the user for a new value.

The list of variables that may be set interactively can be found in

   src/smdr_ioparameters.h

--------------------------------------------------------------------- 
int SMDR_Read_Model_File (char *fileName);

Reads a complete input file specifying a model. Format should be as in
ReferenceModel.dat, although the lines can be in any order, the
trailing semicolons can be included or not, and any line starting with
a '#' is treated as a comment.

If any of the MSbar lagrangian or on-shell parameters are not found in
the file, the missing variables are listed and a fatal error results.
See ReferenceModel.dat for a complete listing of what constitutes a
model file.

---------------------------------------------------------------------
int SMDR_Write_Model_File (char *fileName);

Writes a model file in the same format as ReferenceModel.dat, using the 
current values of the MSbar inputs and the observable outputs. The 
resulting file can be used as an alternative input file for any 
application that could have used ReferenceModel.dat.

--------------------------------------------------------------------- 
void SMDR_Write_Column_Data (FILE *,
                             int n,
                             char *strings[],
                             char *prepend);

This utility prints a pre-defined array of strings to a file stream,
with a prepended character. See the files applications/fig_*.c for
examples of its use.

---------------------------------------------------------------------
---------------------------------------------------------------------
For the remaining seven function pairs, the first form (*_Display_*)
prints to stdout.  The second (*_Write_*) prints to a general file
stream, and prepends a string to each line.  Thus, e.g.,

   SMDR_Display_MSbar_Parameters () 

is equivalent to

   SMDR_Write_MSbar_Parameters (stdout, "");

To write values to a file pointer fp with a prepending '#', use e.g.

   SMDR_Write_MSbar_Parameters (fp, "# ");

---------------------------------------------------------------------
---------------------------------------------------------------------
void SMDR_Display_MSbar_Parameters (void);
void SMDR_Write_MSbar_Parameters (FILE *file, char *prepend);

Print the current scale Q and Standard Model MSbar working parameters 
  SMDR_g3, SMDR_g,  SMDR_gp, 
  SMDR_yt, SMDR_yb, SMDR_ytau, 
  SMDR_yc, SMDR_ys, SMDR_ymu
  SMDR_yu, SMDR_yd, SMDR_ye
  SMDR_lambda

Note that SMDR_v, SMDR_m2, SMDR_Lambda, and SMDR_Delta_alpha_had_5_MZ
are not included here, and have their own separate Write and Display
functions, described below.

---------------------------------------------------------------------
void SMDR_Display_OS_Inputs (void);
void SMDR_Write_OS_Inputs (FILE *fileName, char *prepend);

Print the Standard Model on-shell input parameters.

---------------------------------------------------------------------
void SMDR_Display_Delta_alpha_had5 (void);
void SMDR_Write_Delta_alpha_had5 (FILE *, char *);

Print the value of SMDR_Display_Delta_alpha_had5.

--------------------------------------------------------------------- 
void SMDR_Display_Lambda (void);
void SMDR_Write_Lambda (FILE *, char *);

Print the working value of SMDR_Lambda.

---------------------------------------------------------------------  
void SMDR_Display_m2 (void);
void SMDR_Write_m2 (FILE *, char *);

Print the working value of SMDR_m2.

---------------------------------------------------------------------  
void SMDR_Display_v (void);
void SMDR_Write_v (FILE *, char *);

Print the working value of SMDR_v.

---------------------------------------------------------------------  
void SMDR_Display_Version ();
void SMDR_Write_Version (FILE *, char *);

Print SMDR name and current version.


*********************************************************************  
V. Sample Program
********************************************************************* 

Below is a code fragment that illustrates the use of SMDR. It is based
on the code applications/calc_all.c, provided with the distribution.

#include "smdr.h" /* Required for all applications of SMDR. Headers
                     for TSIL and 3VIL are also included
                     automatically. */

int main (int argc, char *argv[])
{
  char inputFile[50], outputFile[50];
  int do_Veff, interactiveMode;

/* Here is the list of input variables to be read inm consisting of
   the basic MSbar parameters plus SMDR_v_in and
   SMDR_Delta_alpha_had_5_MZ_in: */

  int nVars = 16;
  char *varList[] = {
    "SMDR_Q_in",
    "SMDR_v_in",
    "SMDR_lambda_in",
    "SMDR_g3_in",
    "SMDR_g_in",
    "SMDR_gp_in",
    "SMDR_yt_in",
    "SMDR_yb_in",
    "SMDR_yc_in",
    "SMDR_ys_in",
    "SMDR_yu_in",
    "SMDR_yd_in",
    "SMDR_ytau_in",
    "SMDR_ymu_in",
    "SMDR_ye_in",
    "SMDR_Delta_alpha_had_5_MZ_in"
  };

  /* Define possible command line arguments. In this case there are
     four, all optional. The first and last are of type "toggle", so
     that if specified the corresponding variable will be set equal to
     1. There are also two "string" type arguments, corresponding to
     posssible input and output files. The final array holds pointers
     to the variables that will be set. */
  int nargs = 4;
  char *arglist[] = {"-V","-i","-o","-int"};
  char *argtype[] = {"toggle","string","string","toggle"};
  void *argvar[] = {&do_Veff, inputFile, outputFile, &interactiveMode};

  /* Users should be sure to set default values for any optional
     arguments: */
  do_Veff = 0;
  strcpy (inputFile, "ReferenceModel.dat");
  strcpy (outputFile, "NULL");
  interactiveMode = 0;

  /* This processes the arguments and sets the corresponding values: */  
  SMDR_Process_Arguments (argc, argv, nargs, arglist, argtype, argvar);

  /* If an outputFile has been specified, we will automatically set
     -V: */
  if (0 != strcmp ("NULL", outputFile)) do_Veff = 1; 

  /* Print version information: */
  SMDR_Display_Version ();

  /* Read in values from input file. If -i was not set, default is
     "ReferenceModel.dat", in the current directory. 

     Note that in place of this one could equivalently eliminate the
     declarations of  nVars and varList above, and instead call:

        SMDR_Read_MSbar_Inputs (inputFile);
        SMDR_Read_Value (inputFile, "SMDR_v_in");
        SMDR_Read_Value (inputFile, "SMDR_Delta_alpha_had_5_MZ_in");
  */
  SMDR_Read_Values (inputFile, nVars, varList);

  /* If -int was specified we are in interactive mode; offer the user
     a chance to change any values by hand: 
  */
  if (interactiveMode)
    SMDR_Set_Values_Interactively (nVars, varList);

  /* (Definition of varLits allows all the parameters that were read
     to also be modified interactively.) */

  /* Copies input values into working variables; necessary before 
     displaying them. */
  SMDR_Load_Inputs ();

  printf("\nINPUT PARAMETERS read");
  if (0 == interactiveMode) {
    printf(" from \"%s\"", inputFile);
  } 
  printf(":\n\n");

  /* Print input values to stdout; */
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

  /* Evaluate top-quark pole mass. */
  SMDR_Eval_Mt_pole (SMDR_Mt_EXPT, 1, 4, 2, &SMDR_Mt_pole,
                                            &SMDR_Gammat_pole);

  /* Evaluate Higgs boson pole mass. */
  SMDR_Eval_Mh_pole (160., 2.5, &SMDR_Mh_pole, &SMDR_Gammah_pole);

  /* Evaluate Z boson pole and Breit-Wigner masses */
  SMDR_Eval_MZ_pole (160., 2, &SMDR_MZ_pole, &SMDR_GammaZ_pole,
                     &SMDR_MZ_BreitWigner, &SMDR_GammaZ_BreitWigner);

  /* Evaluate GFermi */
  SMDR_GFermi = SMDR_Eval_GFermi (SMDR_Mt_EXPT, 2);

  /* Evaluate MSbar alpha(MZ), sin^2(thetaW), alphaS(MZ) in
     non-decoupled theory, Sommerfeld fine structure constant alpha,
     alpha(MZ) and sin^2(thetaW) in PDG MSbar scheme with only top
     decoupled.  */
  SMDR_Eval_Gauge (SMDR_Mt_pole, SMDR_Mh_pole, SMDR_MW_BreitWigner);

  /* Evaluate MSbar parameters in 5-quark, 3-lepton QCDQED theory at
  MZ */
  SMDR_Eval_QCDQED_at_MZ (SMDR_MZ_EXPT, SMDR_MZ_EXPT, 5);

  printf("\nOUTPUT QUANTITIES:\n");
 
  if (1 == do_Veff) {
    printf("\n");
    SMDR_Display_m2();
    SMDR_Display_Lambda();
  }

  printf("\n");
  printf("Mt = %Lf; Gammat = %Lf;   (* complex pole *)\n\n",
          SMDR_Mt_pole, SMDR_Gammat_pole);

  printf("Mh = %Lf; Gammah = %Lf;   (* complex pole *)]\n\n",
          SMDR_Mh_pole, SMDR_Gammah_pole);

  ... etc. ...

  SMDR_Timer();
  printf("\nTotal calculation time: %.2f seconds\n", SMDR_Time_Total);

  /* If an outputFile was specified, write a full model file: */
  if (0 != strcmp (outputFile, "NULL"))
    SMDR_Write_Model_File (outputFile);

  return 0;
}

*********************************************************************  
VI. C++ Linkage
********************************************************************* 

The SMDR library functions are all callable "as is" from C++ code. The
header file smdr.h must of course still be included in any C++ source
code that makes use of SMDR.

Some of the utilities that take or use string arguments may produce
(harmless) compiler warnings if string literals are used.  This is
because the C++ type for a string literal is slightly different (const
char[]) than in C (char *).  These warnings can be avoided by
explicitly casting the C++ strings to type (char *). The sample code
C++/calc_all.cpp shows this approach in action.  Note, however, that
the source file applications/calc_all.c will compile just fine under
C++ (albeit with warnings).  Indeed, calc_all.cpp is identical to
calc_all.c as far as SMDR usage is concerned; only the headers and IO
statements have been re-written in a C++ idiom.


********************************************************************* 
End of README.txt
*********************************************************************
