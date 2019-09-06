(* ::Package:: *)

(* Wolfram language version check *)
If[
  $VersionNumber < 12,
  Print["RunningCouplings requires Wolfram Language 12 or later."];
  Abort[]
]

(* unprotect package symbols in case RunningCouplings is double-loaded *)
Unprotect["RunningCouplings`*", "RunningCouplings`Developer`*", "RunningCouplings`SMDR`*", "RunningCouplings`Config`*"];

(* load the package *)
Get["RunningCouplings`RunningCouplings`"]

(* protect all package symbols *)
SetAttributes[
  Flatten[
    Names /@ {"RunningCouplings`*", "RunningCouplings`Developer`*", "RunningCouplings`SMDR`*", "RunningCouplings`Config`*"}
  ] // Evaluate,
  {Protected, ReadProtected}
]
