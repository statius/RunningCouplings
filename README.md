# RunningCouplings

Wolfram Language package calculating the couplings and mixings (including uncertainty) of the standard model (SM) and minimal supersymmetric standard model (MSSM) as functions of the renormalization scale

*RunningCouplings* uses [*SMDR*](https://www.niu.edu/spmartin/SMDR/) (Standard Model in Dimensional Reduction) to extract the Standard Model masses and couplings.

**Code:** [github.com/statius/runningcouplings](https://github.com/statius/runningcouplings)

## Installation

- Download the latest release and unpack it or clone the repository somewhere on the Wolfram Language `$Path` (e.g. the `Applications` folder in `$UserBaseDirectory` for *Mathematica*).

- Load *RunningCouplings* as

  ```mathematica
  Needs @ "RunningCouplings`"
  ```

  On the first initialization, *RunningCouplings* will attempt to build the included copy of *SMDR* from source. This will take some time. If the build fails, options may be specified using

  ```mathematica
  BuildSMDR[
    "CCompiler" -> "<command>",
    "CComplierFlags" -> "<flags>"
    "ArchiveCommand" -> "<command>"
  ]
  ```

  Alternatively, the path to a directory containing precompiled SMDR executables can be (persistantly) specified as

  ```mathematica
  ConfigureRunningCouplings["SMDRExecutableDirectory" -> "<path>"]
  ```

  The automatic *SMDR* compilation procedure has only been tested on macOS.

## Renormalization Procedure

- SM MS-bar couplings are extracted from experimental data at the top quark pole mass scale using the *SMDR* function `calc_fit`. 

- The couplings are then run to the Z boson pole mass scale using the *SMDR* program `calc_run`.

- PMNS matrix parameters and the measured neutrino mass differences are used to construct either a neutrino Yukawa coupling matrix or a neutrino Weinberg operator matrix.

  *Note:* The gauge-fixing dependence of the Higgs VEV (*SMDR* employs the Landau gauge) is neglected when the neutrino Yukawa or Weinberg operator matrix is calculated, as it is minimal at this scale.

- The parameters may be run to a specified renormalization scale using the Wolfram Language function `NDSolve` to solve the SM beta functions at up to two loops. Parameter uncertainties are calculated by performing the renormalization with each parameter individually varied across its uncertain range and summing the resulting deviations in quadrature.

- For renormalization in the MSSM, the SM parameters are first run to a specified scale of supersymmetry, where they are matched to the MSSM at tree level, and converted (at one loop) into the DR-bar renormalization scheme. The parameters may then be run to a specified renormalization scale using `NDSolve` to solve the MSSM beta functions at up to two loops. Parameter uncertainties are calculated as in the SM procedure. In addition, uncertainty due to neglected theshold corrections is estimated by varying the scale of supersymmetry over the interval [mSUSY/2, 2 mSUSY] and calculating the resulting deviations, which are added linearly to the parameter uncertainties.

## Usage

For usage examples, see *usage-examples.nb*.

## Project Information

### Licensing

This project is released under the GPL license.

*SMDR* is released under the GPL license

### Contributions

This package is maintained by Andrew Miller (and primarily created for my own needs). Pull requests and suggestions are always welcomed.

*SMDR* is maintained by S. P. Martin and D. G. Robertson. Please see [arxiv.org/abs/1907.02500](http://arxiv.org/abs/1907.02500) for more information.