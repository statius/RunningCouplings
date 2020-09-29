(* neutrino mass experimental reference parameters *)


(* neutrino mass differences (in eV^2) with 3-sigma uncertainty *)
(* NuFIT group (nu-fit.org) global fit (July 2020) v5.0 *)
$neutrino["NormalHierarchy", "DeltaMassSquared21", "CentralValue"] = 7.42`3*^-5;
$neutrino["NormalHierarchy", "DeltaMassSquared21", "Uncertainty"] = {6.82*^-5, 8.04*^-5} - 7.42*^-5;
$neutrino["NormalHierarchy", "DeltaMassSquared31", "CentralValue"] = 2.514`4*^-3;
$neutrino["NormalHierarchy", "DeltaMassSquared31", "Uncertainty"] = {2.431*^-3, 2.598*^-3} - 2.514*^-3;

$neutrino["InvertedHierarchy", "DeltaMassSquared21", "CentralValue"] = 7.42`3*^-5;
$neutrino["InvertedHierarchy", "DeltaMassSquared21", "Uncertainty"] = {6.82*^-5, 8.04*^-5} - 7.42*^-5;
$neutrino["InvertedHierarchy", "DeltaMassSquared32", "CentralValue"] = -2.497`4*^-3;
$neutrino["InvertedHierarchy", "DeltaMassSquared32", "Uncertainty"] = {-2.583*^-3, -2.412*^-3} + 2.497*^-3;
