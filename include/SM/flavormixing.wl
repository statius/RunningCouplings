(* CKM and PMNS experimental reference parameters *)


(* CKM Wolfenstein parameters with 3-sigma uncertainty *)
(* CKMfitter group (ckmfitter.in2p3.fr) global fit (summer 2019) *)
$ckm["WolfensteinLambda", "CentralValue"] = 0.22484`5;
$ckm["WolfensteinLambda", "Uncertainty"] = {-0.00018, 0.00097};
$ckm["WolfensteinA", "CentralValue"] = 0.823`3;
$ckm["WolfensteinA", "Uncertainty"] = {-0.029, 0.015};
$ckm["WolfensteinRhoBar", "CentralValue"] = 0.157`3;
$ckm["WolfensteinRhoBar", "Uncertainty"] = {-0.019, 0.036};
$ckm["WolfensteinEtaBar", "CentralValue"] = 0.350`3;
$ckm["WolfensteinEtaBar", "Uncertainty"] = {-0.023, 0.027};

(* CKM Jarlskog invariant with 3-sigma uncertainty *)
(* CKMfitter group (ckmfitter.in2p3.fr) global fit (summer 2019) *)
$ckm["JarlskogInvariant", "CentralValue"] = 3.06`3*^-5;
$ckm["JarlskogInvariant", "Uncertainty"] = {-0.18*^-5, 0.27*^-5};

(* PMNS standard parameters with 3-sigma uncertainty *)
(* NuFIT group (nu-fit.org) global fit (July 2020) v5.0 *)
$pmns["NormalHierarchy", "SineSquaredTheta12", "CentralValue"] = 0.304`3;
$pmns["NormalHierarchy", "SineSquaredTheta12", "Uncertainty"] = {0.269, 0.343} - 0.304;
$pmns["NormalHierarchy", "Theta12", "CentralValue"] = 33.44`4 Degree;
$pmns["NormalHierarchy", "Theta12", "Uncertainty"] = ({31.27, 35.86} - 33.44) Degree;
$pmns["NormalHierarchy", "SineSquaredTheta23", "CentralValue"] = 0.570`3;
$pmns["NormalHierarchy", "SineSquaredTheta23", "Uncertainty"] = {0.407, 0.618} - 0.570;
$pmns["NormalHierarchy", "Theta23", "CentralValue"] = 49.0`3 Degree;
$pmns["NormalHierarchy", "Theta23", "Uncertainty"] = ({39.6, 51.8} - 49.0) Degree;
$pmns["NormalHierarchy", "SineSquaredTheta13", "CentralValue"] = 0.02221`4;
$pmns["NormalHierarchy", "SineSquaredTheta13", "Uncertainty"] = {-0.02034, 0.02430} - 0.02221;
$pmns["NormalHierarchy", "Theta13", "CentralValue"] = 8.57`3 Degree;
$pmns["NormalHierarchy", "Theta13", "Uncertainty"] = ({8.20, 8.97} - 8.57) Degree;
$pmns["NormalHierarchy", "Delta", "CentralValue"] = 195`3 Degree;
$pmns["NormalHierarchy", "Delta", "Uncertainty"] = ({107, 403} - 195) Degree;

$pmns["InvertedHierarchy", "SineSquaredTheta12", "CentralValue"] = 0.304`3;
$pmns["InvertedHierarchy", "SineSquaredTheta12", "Uncertainty"] = {0.269, 0.343} - 0.304;
$pmns["InvertedHierarchy", "Theta12", "CentralValue"] = 33.44`4 Degree;
$pmns["InvertedHierarchy", "Theta12", "Uncertainty"] = ({31.27, 35.87} - 33.44) Degree;
$pmns["InvertedHierarchy", "SineSquaredTheta23", "CentralValue"] = 0.575`3;
$pmns["InvertedHierarchy", "SineSquaredTheta23", "Uncertainty"] = {0.411, 0.621} - 0.575;
$pmns["InvertedHierarchy", "Theta23", "CentralValue"] = 49.3`3 Degree;
$pmns["InvertedHierarchy", "Theta23", "Uncertainty"] = ({39.9, 52.0} - 49.3) Degree;
$pmns["InvertedHierarchy", "SineSquaredTheta13", "CentralValue"] = 0.02240`4;
$pmns["InvertedHierarchy", "SineSquaredTheta13", "Uncertainty"] = {0.02053, 0.02436} - 0.02240;
$pmns["InvertedHierarchy", "Theta13", "CentralValue"] = 8.61`3 Degree;
$pmns["InvertedHierarchy", "Theta13", "Uncertainty"] = ({8.24, 8.98} - 8.61) Degree;
$pmns["InvertedHierarchy", "Delta", "CentralValue"] = 286`3 Degree;
$pmns["InvertedHierarchy", "Delta", "Uncertainty"] = ({192, 360} - 286) Degree;
