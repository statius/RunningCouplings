(* CKM and PMNS experimental reference parameters *)


(* CKM Wolfenstein parameters with 3-sigam uncertainty *)
(* CKMfitter group (ckmfitter.in2p3.fr) global fit (summer 2018) *)
$ckm["WolfensteinLambda", "CentralValue"] = 0.22475`5;
$ckm["WolfensteinLambda", "Uncertainty"] = {-0.00018, 0.00106};
$ckm["WolfensteinA", "CentralValue"] = 0.840`3;
$ckm["WolfensteinA", "Uncertainty"] = {-0.043, 0.016};
$ckm["WolfensteinRhoBar", "CentralValue"] = 0.158`3;
$ckm["WolfensteinRhoBar", "Uncertainty"] = {-0.020, 0.036};
$ckm["WolfensteinEtaBar", "CentralValue"] = 0.349`3;
$ckm["WolfensteinEtaBar", "Uncertainty"] = {-0.025, 0.029};

(* CKM Jarlskog invariant with 3-sigma uncertainty *)
(* CKMfitter group (ckmfitter.in2p3.fr) global fit (summer 2018) *)
$ckm["JarlskogInvariant", "CentralValue"] = 3.17`3*^-5;
$ckm["JarlskogInvariant", "Uncertainty"] = {-0.34*^-5, 0.29*^-5};

(* PMNS standard parameters with 3-sigma uncertainty *)
(* NuFIT group (nu-fit.org) global fit (July 2019) v4.1 *)
$pmns["NormalHierarchy", "SineSquaredTheta12", "CentralValue"] = 0.310`3;
$pmns["NormalHierarchy", "SineSquaredTheta12", "Uncertainty"] = {-0.035, 0.040};
$pmns["NormalHierarchy", "Theta12", "CentralValue"] = 33.82`4 Degree;
$pmns["NormalHierarchy", "Theta12", "Uncertainty"] = {-2.21, 2.45} Degree;
$pmns["NormalHierarchy", "SineSquaredTheta23", "CentralValue"] = 0.563`3;
$pmns["NormalHierarchy", "SineSquaredTheta23", "Uncertainty"] = {-0.130, 0.046};
$pmns["NormalHierarchy", "Theta23", "CentralValue"] = 48.6`3 Degree;
$pmns["NormalHierarchy", "Theta23", "Uncertainty"] = {-7.5, 2.7} Degree;
$pmns["NormalHierarchy", "SineSquaredTheta13", "CentralValue"] = 0.02237`4;
$pmns["NormalHierarchy", "SineSquaredTheta13", "Uncertainty"] = {-0.00193, 0.00198};
$pmns["NormalHierarchy", "Theta13", "CentralValue"] = 8.60`3 Degree;
$pmns["NormalHierarchy", "Theta13", "Uncertainty"] = {-0.38, 0.38} Degree;
$pmns["NormalHierarchy", "Delta", "CentralValue"] = 221`3 Degree;
$pmns["NormalHierarchy", "Delta", "Uncertainty"] = {-77, 136} Degree;

$pmns["InvertedHierarchy", "SineSquaredTheta12", "CentralValue"] = 0.310`3;
$pmns["InvertedHierarchy", "SineSquaredTheta12", "Uncertainty"] = {-0.035, 0.040};
$pmns["InvertedHierarchy", "Theta12", "CentralValue"] = 33.82`4 Degree;
$pmns["InvertedHierarchy", "Theta12", "Uncertainty"] = {-2.21, 2.45} Degree;
$pmns["InvertedHierarchy", "SineSquaredTheta23", "CentralValue"] = 0.565`3;
$pmns["InvertedHierarchy", "SineSquaredTheta23", "Uncertainty"] = {-0.129, 0.045};
$pmns["InvertedHierarchy", "Theta23", "CentralValue"] = 48.8`3 Degree;
$pmns["InvertedHierarchy", "Theta23", "Uncertainty"] = {-7.4, 2.5} Degree;
$pmns["InvertedHierarchy", "SineSquaredTheta13", "CentralValue"] = 0.02259`4;
$pmns["InvertedHierarchy", "SineSquaredTheta13", "Uncertainty"] = {-0.00195, 0.00198};
$pmns["InvertedHierarchy", "Theta13", "CentralValue"] = 8.64`3 Degree;
$pmns["InvertedHierarchy", "Theta13", "Uncertainty"] = {-0.38, 0.38} Degree;
$pmns["InvertedHierarchy", "Delta", "CentralValue"] = 282`3 Degree;
$pmns["InvertedHierarchy", "Delta", "Uncertainty"] = {-77, 66} Degree;
