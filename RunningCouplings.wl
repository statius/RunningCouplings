(* ::Package:: *)

(* ::Title:: *)
(*RunningCouplings package*)


(* ::Subsubsection:: *)
(*metadata*)


(* : Title : RunningCouplings *)
(* : Author : Andrew Miller < amiller@physics.umn.edu > *)
(* : Context : RunningCouplings` *)
(* : Version : 0.1.0 *)
(* : Date : 2019-08-01 *)
(* : Mathematica Version : 12.0 *)
(* : Copyright : (c) 2019 Andrew Miller *)


(* ::Title::GrayLevel[0]:: *)
(*context RunningCouplings`*)


(* ::Subsubsection::GrayLevel[0]::Closed:: *)
(*begin package context*)


BeginPackage["RunningCouplings`"];


(* ::Section:: *)
(*usage messages*)


(* ::Subsection::Closed:: *)
(*configuration*)


ConfigureRunningCouplings::usage = "ConfigureRunningCouplings [\!\(\*SubscriptBox[\(name\), \(1\)]\) \[Rule] \!\(\*SubscriptBox[\(value\), \(1\)]\), \!\(\*SubscriptBox[\(name\), \(2\)]\) \[Rule] \!\(\*SubscriptBox[\(value\), \(2\)]\), \[Ellipsis]] sets the configuration options for RunningCouplings and stores them persistantly.
ConfigureRunningCouplings[] returns the current configuration.";


(* ::Subsection::Closed:: *)
(*utilities*)


TakagiFactorization::usage = "TakagiFactorization[m] gives the Takagi factorization for a numerical symmetric matrix m as a list of matrices {u, d, e}, where d is a diagonal matrix, m can be written as u.d.Transpose[u], and e is a matrix of residuals."


ComplexAround::usage = "ComplexAround[x, \[Delta]] attempts to extend the functionality of Around to numeric complex numbers. It has the same behaviors as Around."


(* ::Subsection::Closed:: *)
(*SMDR interface*)


BuildSMDR::usage = "BuildSMDR[opts] attempts to build the included SMDR package.";


SMDRCalcFit::usage = "SMDRCalcFit[opts] runs the SMDR function calc_fit.
SMDRCalcFit[args, opts] runs calc_fit for the association of input values args.";


SMDRCalcRGrun::usage = "SMDRCalcRGrun[q, opts] runs the SMDR function calc_RGrun to the RG scale q.
SMDRCalcFit[q, args, opts] runs calc_RGrun for the association of input values args.";


(* ::Subsection::Closed:: *)
(*mixing*)


MixingMatrix::usage = "MixingMatrix[\"name\", params, opts] calculates the mixing matrix \"name\" from the list of parameters params.
MixingMatrix[\"name\", params] calculates the matrix from the association of parameters params.
MixingMatrix[\"name\", {\!\(\*SubscriptBox[\(m\), \(1\)]\), \!\(\*SubscriptBox[\(m\), \(2\)]\)}, opts] calculates matrix from the Yukawa or mass matrices \!\(\*SubscriptBox[\(m\), \(1\)]\) and \!\(\*SubscriptBox[\(m\), \(2\)]\).";


MixingParameters::usage = "MixingParameters[\"name\", m, opts] calculates the parameters of mixing matrix \"name\" from the matrix m.
MixingParameters[\"name\", {\!\(\*SubscriptBox[\(m\), \(1\)]\), \!\(\*SubscriptBox[\(m\), \(2\)]\)}, opts] calculates the parameters from the Yukawa or mass matrices \!\(\*SubscriptBox[\(m\), \(1\)]\) and \!\(\*SubscriptBox[\(m\), \(2\)]\).
MixingParameters[\"name\", params, \"Parametrization\" \[Rule] \"para\"] calculates the parameters in the parametrization \"para\" from the association of parameters params.";


(* ::Subsection::Closed:: *)
(*renormalization*)


RenormalizeModel::usage = "RenormalizeModel[{\"model\", \"method\"}, q, params] runs the method \"method\" to renormalize the parameters params of the model \"model\" to the scale q.";


(* ::Subsection::Closed:: *)
(*RunningCouplings*)


RunningCouplings::usage = "RunningCouplings[\"model\", q, opts] calculates the running couplings in model \"model\" at the scale q using options opts.";


(* ::Subsection::Closed:: *)
(*running quantities*)


RunningGaugeCoupling::usage = "RunningGaugeCoupling[{\"model\", i}, q, opts] gives the \!\(\*SuperscriptBox[\(i\), \(th\)]\) gauge coupling at the scale q in the model \"model\".";


RunningYukawaCoupling::usage = "RunningYukawaCoupling[{\"model\", \"family\"}, q, opts] gives the Yukawa couplings for the family \"family\" at the scale q in the model \"model\".
YukawaCoupling[{\"model\", \"family\", i}, q, opts] gives the Yukawa coupling for the \!\(\*SuperscriptBox[\(i\), \(th\)]\) generation.
YukawaCoupling[{\"model\", \"name\"}, q, opts] gives the Yukawa coupling the fermion \"name\".";


RunningWeinbergOperator::usage = "RunningWeinbergOperator[\"model\", q, opts] gives the neutrino Weinberg operators at the scale q in the model \"model\".
RunningWeinbergOperator[{\"model\", i}, q, opts] gives the operator for the \!\(\*SuperscriptBox[\(i\), \(th\)]\) generation.";


RunningMixingParameter::usage = "RunningMixingParameters[{\"model\", \"matrix\"}, q, opts] gives the parameters of the mixing matrix \"matrix\" at the scale q in the model \"model\".
RunningMixingParameters[{\"model\", \"matrix\", \"param\"}, q, opts] gives the parameter \"param\".";


(* ::Section:: *)
(*package information*)


(* ::Subsubsection::Closed:: *)
(*general*)


`Developer`$Version = "0.1.0 (01 August 2019)";


`Developer`$VersionNumber = StringReplace[`Developer`$Version, "." ~~ Except["."] .. ~~ EndOfString :> ""];


`Developer`$ReleaseNumber = StringSplit[`Developer`$Version, {" ", "."}][[3]];


`Developer`$CreationDate := DateObject @ Last @ StringSplit[`Developer`$Version, {"(", ")"}]


`Developer`$PackageDirectory = DirectoryName @ $InputFileName;


(* ::Subsubsection::Closed:: *)
(*SMDR*)


`SMDR`$BaseDirectory = FileNameJoin @ {`Developer`$PackageDirectory, "smdr"};


`SMDR`$Version = StringDelete[Last @ Sort[StringCases[FileNames[All, `SMDR`$BaseDirectory], "smdr-" ~~ __] // Flatten], "smdr-"];


`SMDR`$SourceDirectory = FileNameJoin @ {`SMDR`$BaseDirectory, "smdr-" <> `SMDR`$Version};


`SMDR`$ExecutableDirectory = FileNameJoin @ {`SMDR`$SourceDirectory, "bin"};


`SMDR`$IncludeDirectory = FileNameJoin @ {`SMDR`$SourceDirectory, "include"};


`SMDR`$LibraryDirectory = FileNameJoin @ {`SMDR`$SourceDirectory, "lib"};


(* ::Chapter:: *)
(*context RunningCouplings`Private`*)


(* ::Subsubsection::GrayLevel[0]::Closed:: *)
(*begin private context*)


Begin @ "`Private`";


(* ::Subsubsection::GrayLevel[0]::Closed:: *)
(*imports*)


Needs @ "DifferentialEquations`InterpolatingFunctionAnatomy`"


Needs @ "JLink`"


(* ::Section:: *)
(*configuration settings*)


(* ::Subsection::Closed:: *)
(*configuration file*)


RunningCouplings`Config`$ConfigDirectory = If[
  SameQ[ParentDirectory @ RunningCouplings`Developer`$PackageDirectory, FileNameJoin @ {$UserBaseDirectory, "Applications"}], 
  FileNameJoin @ {$UserBaseDirectory, "ApplicationData", FileNameTake @ RunningCouplings`Developer`$PackageDirectory},
  RunningCouplings`Developer`$PackageDirectory
 ];


RunningCouplings`Config`$ConfigFile = FileNameJoin @ {RunningCouplings`Config`$ConfigDirectory, "config.m"};


(* ::Subsection::Closed:: *)
(*configuration options*)


$configDefault = Association[
  "SMDRExecutableDirectory" :> RunningCouplings`SMDR`$ExecutableDirectory,
  "NumericTolerance" -> 10^(-24),
  "DoubletSingularVectors" -> Right
];


$configAlternatives = Association[
  "SMDRExecutableDirectory" -> _String | _File,
  "NumericTolerance" -> _Integer,
  "DoubletSingularVectors" -> Left | "Left" | Right | "Right"
];


(* ::Subsection:: *)
(*configuration API*)


(* ::Subsubsection::Closed:: *)
(*ConfigureRunningCouplings*)


ConfigureRunningCouplings::optx = "Unknown configuration option ``.";
ConfigureRunningCouplings::valng = "Value of configuration option `1` is not one of `2`.";


SyntaxInformation @ ConfigureRunningCouplings = {"ArgumentsPattern" -> {OptionsPattern[]}, "OptionNames" -> Keys @ $configDefault};


ConfigureRunningCouplings[rules__Rule] := (
  iConfigureRunningCouplings[rules];
  Export[RunningCouplings`Config`$ConfigFile, RunningCouplings`Config`$Config, "Package"];
  Normal @ RunningCouplings`Config`$Config
)


ConfigureRunningCouplings[] := Normal @ RunningCouplings`Config`$Config


(* ::Subsubsection::Closed:: *)
(*iConfigureRunningCouplings*)


iConfigureRunningCouplings[] := Module[
  {goodopts, badopts},

  goodopts = KeyTake[RunningCouplings`Config`$Config, Keys @ $configDefault];
  badopts = KeyDrop[RunningCouplings`Config`$Config, Keys @ $configDefault];

  Scan[Message[ConfigureRunningCouplings::optx, First @ #] &, Normal @ badopts];
  KeyValueMap[
    If[!MatchQ[#2, $configAlternatives @ #1],
      Message[ConfigureRunningCouplings::valng, #1 -> #2, TextString[List @@ $configAlternatives @ #1]];
      KeyDropFrom[goodopts, #1]
    ] &,
    goodopts
   ];
   
  RunningCouplings`Config`$Config = Join[$configDefault, goodopts]
  
]


iConfigureRunningCouplings[rules__Rule] := Module[
  {goodopts, badopts},

  goodopts = KeyTake[<|rules|>, Keys @ $configDefault];
  badopts = KeyDrop[<|rules|>, Keys @ $configDefault];

  Scan[Message[ConfigureRunningCouplings::optx, First @ #] &, Normal @ badopts];
  KeyValueMap[
    If[!MatchQ[#2, $configAlternatives @ #1],
      Message[ConfigureRunningCouplings::valng, #1 -> #2, TextString[List @@ $configAlternatives @ #1]];
      KeyDropFrom[goodopts, #1]
    ] &,
    goodopts
   ];
   
  AppendTo[RunningCouplings`Config`$Config, goodopts]

]


(* ::Subsection::Closed:: *)
(*initial setting*)


If[
  FileExistsQ @ RunningCouplings`Config`$ConfigFile, 
  RunningCouplings`Config`$Config = Normal @ Import[RunningCouplings`Config`$ConfigFile, "Package"], 
  RunningCouplings`Config`$Config = $configDefault
];


(* ::Subsection::Closed:: *)
(*aliases*)


$config := RunningCouplings`Config`$Config


$smdrExecutableDirectory := $config @ "SMDRExecutableDirectory"


$context = Context[];


(* ::Section:: *)
(*utilities*)


(* ::Subsubsection::Closed:: *)
(*relativeFilePath*)


InstallJava[];


relativeFilePath[base_String, target_String] := JavaBlock[

  With[
    {
     basepath = JavaNew["java.io.File", base],
     targetpath = JavaNew["java.io.File", target]
    },
    
    basepath @ toPath[] @ relativize[targetpath @ toPath[]] @ toString[]
  
  ]
]


(* ::Subsubsection::Closed:: *)
(*fortranNumberString*)


fortranNumberString = NumberString ~~ ((("e" | "E") ~~ ("+" | "-" | "") ~~ DigitCharacter ..) | "") 


(* ::Subsubsection::Closed:: *)
(*toPreciseNumber*)


toPreciseNumber[s_String] := SetPrecision[
  Interpreter["Number"][s], 
  StringLength @ StringDelete[s, {"-", StartOfString ~~ "0." ~~ ("0" ...), ".", "e" ~~ __ ~~ EndOfString}, IgnoreCase -> True]
]


toPreciseNumber[n : (_Real | _Integer), ref_Real] := SetPrecision[
  n, 
  Length[RealDigits[ref][[1]]] + RealDigits[n][[2]] - RealDigits[ref][[2]]
]


toPreciseNumber[n_, ref_] := n


(* ::Subsubsection::Closed:: *)
(*toRoundedNumber*)


SetAttributes[toRoundedNumber, Listable];


toRoundedNumber[i_Integer] := i


toRoundedNumber[n_Real] := With[
  {digits = Length @ #1 - #2 & @@ (RealDigits @ n)},
  
  Round[n, 10^(-digits)]

]  


toRoundedNumber[s_String] := With[
  {digits = Total[StringCases[s, {"." ~~ val : DigitCharacter .. :> - StringLength[val], "e" ~~ val : ( ("+" | "-") ~~ DigitCharacter .. ) ~~ EndOfString :> ToExpression @ val}] // Flatten]},

  Round[ToExpression @ s, 10^digits]

]


(* ::Subsubsection::Closed:: *)
(*roundTo*)


roundTo[n_? NumericQ, ref_? NumericQ] := Round[n, (Length @ #1 - #2) & @@ (RealDigits @ ref)]


(* ::Subsubsection::Closed:: *)
(*numericMatrixQ*)


numericMatrixQ = MatrixQ[#, NumericQ] &


(* ::Subsubsection::Closed:: *)
(*numericSymmetricMatrixQ*)


Options @ numericSymmetricMatrixQ = Options @ SymmetricMatrixQ;


numericSymmetricMatrixQ[m_, opts : OptionsPattern[]] := And[MatrixQ[m, NumericQ], SymmetricMatrixQ[m, opts]]


(* ::Subsubsection::Closed:: *)
(*orderColumns*)


orderColumns[m_? MatrixQ, Reverse] := Transpose[Reverse @ Transpose[m]]


orderColumns[m_? MatrixQ, order_List] /; ContainsExactly[order, Range[Length @ m]] := Transpose[Transpose[m][[order]]]


(* ::Subsubsection::Closed:: *)
(*flattenAssociation*)


flattenAssociation[a_Association] := Association[
  Normal[a] //. (k_ -> v_Association) :> Normal[KeyMap[Append[{k}, #] &, v]]
]


(* ::Subsubsection::Closed:: *)
(*nestedMerge*)


nestedMerge[a : {__Association}] := Merge[a, nestedMerge]


nestedMerge[a_] := Last @ a


(* ::Subsubsection::Closed:: *)
(*nestAssociation*)


nestAssociation[a_Association] := nestedMerge[
  KeyValueMap[
    Association[
      ReplaceAll[
        #1, 
        {
         {k1_, k2_} :> (k1 -> <|k2 -> #2|>), 
         {k1_, k2_, k3__} :> (k1 -> <|{k2, k3} -> #2|>), 
         {k_} :> (k -> #2), 
         _ :> (#1 -> #2)
        }
      ]
    ] &, 
    a
  ]
]


(* ::Subsubsection::Closed:: *)
(*aroundMod*)


aroundMod[Around[v_, u_], n_, d_:0] := Around[Mod[v, n, d], u]


aroundMod[m_, n_, d_:0] := Mod[m, n, d]


(* ::Subsubsection::Closed:: *)
(*TakagiFactorization*)


Options @ TakagiFactorization = Options @ SingularValueDecomposition;


SyntaxInformation @ TakagiFactorization = {"ArgumentsPattern" -> {_, OptionsPattern[]}, "OptionNames" -> Keys @ Options @ TakagiFactorization};


TakagiFactorization::symx = "Warning: matrix `1` is not symmetric within numeric precision.";


(* algorithm following https://mathoverflow.net/questions/125960/diagonalizing-a-complex-symmetric-matrix *)

TakagiFactorization[
  m_? numericMatrixQ, 
  opts : OptionsPattern[]
] := With[
  {svd = SingularValueDecomposition[N @ m, opts]},
  
  If[! SymmetricMatrixQ @ m, Message[TakagiFactorization::symx, m]];
  
  With[
    {d = ConjugateTranspose[svd[[1]]].N[m].Conjugate[svd[[1]]]},
    
    With[
      {v = svd[[1]].DiagonalMatrix[Exp[I Diagonal[Arg @ d] / 2]]},
      
      {v, svd[[2]], Abs[m - v.svd[[2]].Transpose[v]]}
    
    ]
  ]
]


(* ::Subsubsection:: *)
(*ComplexAround*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


SyntaxInformation @ ComplexAround = SyntaxInformation @ Around;


(* ::Subsubsubsection::Closed:: *)
(*single argument*)


ComplexAround[v_] := v


(* ::Subsubsubsection::Closed:: *)
(*real case*)


ComplexAround[
  v_? NumericQ /; Im[v] == 0,
  u : (_? NumericQ | {_? NumericQ, _? NumericQ}) /; AllTrue[{u} // Flatten, Im[#] == 0 &]
] := Around[Re @ v, Re @ u]


(* ::Subsubsubsection::Closed:: *)
(*properties*)


ComplexAround[
  v_? NumericQ,
  (_? NumericQ | {_? NumericQ, _? NumericQ})
]["Value"] := v


ComplexAround[
  _? NumericQ,
  u : (_? NumericQ | {_? NumericQ, _? NumericQ})
]["Uncertainty"] := u


(* ::Subsubsubsection::Closed:: *)
(*complex number behaviors*)


ComplexAround[
  v_? NumericQ,
  u_? NumericQ
] /; Im[u] < 0 := ComplexAround[v, Conjugate @ u]


ComplexAround[
  v_? NumericQ,
  u : {_? NumericQ, _? NumericQ}
] /; AnyTrue[Im[u], Negative] := ComplexAround[v, u /. c_Complex /; Im[c] < 0 :> Conjugate @ c]


ComplexAround /: func_[ComplexAround[v_? NumericQ, u : (_? NumericQ | {_? NumericQ, _? NumericQ})]] /; MatchQ[func, (Re | Im)] := Around[func @ v, func @ u]


ComplexAround /: Conjugate[ComplexAround[v_? NumericQ, u : (_? NumericQ | {_? NumericQ, _? NumericQ})]] := ComplexAround[Conjugate @ v, u]


ComplexAround /: Abs[c : ComplexAround[_? NumericQ, (_? NumericQ | {_? NumericQ, _? NumericQ})]] := Sqrt[Re[c]^2 + Im[c]^2]


ComplexAround /: Arg[c : ComplexAround[_? NumericQ, (_? NumericQ | {_? NumericQ, _? NumericQ})]] := With[
  {
   atan = Piecewise[
     {
      {ArcTan[#1 / #2], #2["Value"] > 0},
      {ArcTan[#1 / #2] + \[Pi], And[#2["Value"] < 0, #1["Value"] >= 0]},
      {ArcTan[#1 / #2] - \[Pi], And[#2["Value"] < 0, #1["Value"] < 0]},
      {\[Pi] / 2, And[#2["Value"] == 0, #1["Value"] > 0]},
      {-\[Pi] / 2, And[#2["Value"] == 0, #1["Value"] < 0]}
     },
     0
   ] &
  },
  
  atan[Im @ c, Re @ c]

]


(* ::Subsubsubsection::Closed:: *)
(*arithmetic operations*)


ComplexAround /: Plus[
  x : (_? NumericQ | _Around | _ComplexAround),
  y : ComplexAround[_, _]
] := With[
  {
   re = Re @ x + Re @ y,
   im = Im @ x + Im @ y
  },
  
  ComplexAround[
    Replace[re, a_Around :> a @ "Value"] + I Replace[im, a_Around :> a @ "Value"],
    Replace[re, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}] + I Replace[im, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}]
  ]
  
]


ComplexAround /: Times[
  x : (_? NumericQ | _Around | _ComplexAround),
  y : ComplexAround[_, _]
] := With[
  {
   re = (Re[x] Re[y] - Im[x] Im[y]),
   im = (Im[x] Re[y] + Re[x] Im[y])
  },
  
  ComplexAround[
    Replace[re, a_Around :> a @ "Value"] + I Replace[im, a_Around :> a @ "Value"],
    Replace[re, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}] + I Replace[im, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}]
  ]
  
]


ComplexAround /: Power[x_ComplexAround, y_] := With[
  {
   re = (Re[x]^2 + Im[x]^2)^(y / 2) Cos[y Arg[x]],
   im = (Re[x]^2 + Im[x]^2)^(y / 2) Sin[y Arg[x]]
  },
  
  ComplexAround[
    Replace[re, a_Around :> a @ "Value"] + I Replace[im, a_Around :> a @ "Value"],
    Replace[re, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}] + I Replace[im, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}]
  ]
  
]


(* ::Subsubsubsection::Closed:: *)
(*boxes*)


ComplexAround /: MakeBoxes[
  c : ComplexAround[_? NumericQ, (_? NumericQ | {_? NumericQ, _? NumericQ})], 
  form : (StandardForm | TraditionalForm)
] := With[
  {
   re = Replace[ToBoxes[Re @ c, form], InterpretationBox[t_TemplaceBox, _] :> t],
   im = Replace[ToBoxes[Im @ c, form], InterpretationBox[t_TemplaceBox, _] :> t]
  },
  
  InterpretationBox[RowBox @ {re, "+", im, "\[ImaginaryI]"}, c]
]


(* ::Subsubsection::Closed:: *)
(*NAroundReplace*)


Options @ NAroundReplace = {"UncertaintyFunction" -> Automatic};


SyntaxInformation @ NAroundReplace = {"ArgumentsPattern" -> {_, {__}, OptionsPattern[]}, "OptionNames" -> Keys @ Options @ NAroundReplace};


NAroundReplace[expr_, rules : {__Rule}, opts : OptionsPattern[]] := With[
  {
   values = rules /. a_Around :> a @ "Value",
   keys = Cases[rules, ((k_ -> _Around) | (k_ :> _Around)) :> k],
   ufunc = OptionValue @ "UncertaintyFunction" /. Automatic -> (Sqrt[{Total[Select[#, Negative]^2], Total[Select[#, Positive]^2]}] &)
  },
  
  With[
    {
     centralvalue = expr /. values,
     uncertainty = Flatten[
       Map[
         {
          expr /. <|values, # -> Min[Lookup[rules, #]["Interval"]]|>, 
          expr /. <|values, # -> Max[Lookup[rules, #]["Interval"]]|>
         } &, 
         keys
       ],
       1
     ]
    },
    
    If[
      Length @ uncertainty > 0,
      
      If[
        NumericQ @ centralvalue,
        
        ComplexAround[
          centralvalue, 
            (ufunc[Re @ uncertainty - Re @ centralvalue] /. {x_, y_} /; x == y :> x) 
          + (ufunc[Im @ uncertainty - Im @ centralvalue] /. {x_, y_} /; x == y :> x) I
        ],
        
        ReplacePart[
          centralvalue,
          Map[
            With[
              {v = Extract[centralvalue, #]},
      
              # -> ComplexAround[
                v, 
                  (ufunc[Re @ uncertainty[[All, Sequence @@ #]] - Re @ v] /. {x_, y_} /; x == y :> x) 
                + (ufunc[Im @ uncertainty[[All, Sequence @@ #]] - Im @ v] /. {x_, y_} /; x == y :> x) I
              ]
        
            ] &,
            Flatten[Position[centralvalue, #] & /@ DeleteDuplicates[Cases[centralvalue, _? NumericQ, \[Infinity]]], 1]
          ]
        ]
      ],
      
      centralvalue
    ]
  
  ]
]


(* ::Subsubsection::Closed:: *)
(*NAroundApply*)


Options @ NAroundApply = {"UncertaintyFunction" -> Automatic};


SyntaxInformation @ NAroundApply = {"ArgumentsPattern" -> {_, {__}, OptionsPattern[]}, "OptionNames" -> Keys @ Options @ NAroundApply};


NAroundApply[func_, args_List, opts : OptionsPattern[]] := With[
  {
   rules = AssociationThread[Range[Length @ args] -> args],
   keys = Position[args, _Around, 1] // Flatten,
   values = (# /. a_Around :> a @ "Value") &,
   ufunc = OptionValue @ "UncertaintyFunction" /. Automatic -> (Sqrt[{Total[Select[#, Negative]^2], Total[Select[#, Positive]^2]}] &)
  },
  
  With[
    {
     centralvalue = func @@ (values @ args),
     uncertainty = Flatten[
       Map[
         {
          func @@ Values[<|values @ rules, # -> Min[Lookup[rules, #]["Interval"]]|>], 
          func @@ Values[<|values @ rules, # -> Max[Lookup[rules, #]["Interval"]]|>]
         } &, 
         keys
       ],
       1
     ]
    },
    
    If[
      Length @ uncertainty > 0,
      
      If[
        NumericQ @ centralvalue,
        
        ComplexAround[
          centralvalue, 
            (ufunc[Re @ uncertainty - Re @ centralvalue] /. {x_, y_} /; x == y :> x) 
          + (ufunc[Im @ uncertainty - Im @ centralvalue] /. {x_, y_} /; x == y :> x) I
        ],
      
        ReplacePart[
          centralvalue,
          Map[
            With[
              {v = Extract[centralvalue, #]},
      
              # -> ComplexAround[
                v, 
                  (ufunc[Re @ uncertainty[[All, Sequence @@ #]] - Re @ v] /. {x_, y_} /; x == y :> x) 
                + (ufunc[Im @ uncertainty[[All, Sequence @@ #]] - Im @ v] /. {x_, y_} /; x == y :> x) I
              ]
        
            ] &,
            Flatten[Position[centralvalue, #] & /@ DeleteDuplicates[Cases[centralvalue, _? NumericQ, \[Infinity]]], 1]
          ]
        ]
      ],
      
      centralvalue
    ]
  
  ]
]
  


(* ::Subsubsection::Closed:: *)
(*ComplexAroundReplace*)


SyntaxInformation @ ComplexAroundReplace = {"ArgumentsPattern" -> {_, _}};


ComplexAroundReplace[
  expr_,
  rules : {__Rule}
] := With[
  {
   re = (Re[expr] // ComplexExpand) /. rules,
   im = (Im[expr] // ComplexExpand) /. rules
  },
  
  If[
    And[MatchQ[re, (_Around | _? NumericQ)], MatchQ[im, (_Around | _? NumericQ)]],
    
    ComplexAround[
      Replace[re, a_Around :> a @ "Value"] + I Replace[im, a_Around :> a @ "Value"],
      Replace[re, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}] + I Replace[im, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}]
    ],
    
    ReplacePart[
      re,
      Map[
        With[
          {
           r = Extract[re, #],
           i = Extract[im, #]
          },
      
          # -> ComplexAround[
            Replace[r, a_Around :> a @ "Value"] + I Replace[i, a_Around :> a @ "Value"],
            Replace[r, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}] + I Replace[i, {a_Around :> a @ "Uncertainty", _? NumericQ :> 0.0}]
          ]
          
        ] &,
        Union[
          Flatten[Position[re, #] & /@ DeleteDuplicates[Cases[re, _Around, \[Infinity]]], 1],
          Flatten[Position[im, #] & /@ DeleteDuplicates[Cases[im, _Around, \[Infinity]]], 1]
        ]
      ]
    ]
    
  ]
]


(* ::Section:: *)
(*SMDR interface*)


(* ::Subsection:: *)
(*compile SMDR*)


(* ::Subsubsection::Closed:: *)
(*getCCompiler*)


getCCompiler[] := StringCases[
  FindList[FileNameJoin @ {RunningCouplings`SMDR`$SourceDirectory, "Makefile"}, "CC"],
  StartOfString ~~ (WhitespaceCharacter ...) ~~ "CC" ~~ (WhitespaceCharacter ...) ~~ "=" ~~ c__ :> StringTrim @ c
] // Flatten // First


(* ::Subsubsection::Closed:: *)
(*getCCompilerFlags*)


getCCompilerFlags[compiler_] := Module[{makefile, flags},

  makefile = OpenRead[FileNameJoin @ {RunningCouplings`SMDR`$SourceDirectory, "Makefile"}];
  
  Find[makefile, compiler];
  
  flags = ReadLine @ makefile;
  flags = If[
    StringQ @ flags,
    StringCases[flags, StartOfString ~~ (WhitespaceCharacter ...) ~~ "SMDR_OPT" ~~ (WhitespaceCharacter ...) ~~ "=" ~~ f : (Except["#"] ...) :> StringTrim @ f] // First,
    ""
  ];
  
  Close @ makefile;
  
  If[
    StringContainsQ[flags, Whitespace],
    "\"" <> flags <> "\"",
    flags
  ]
  
]


(* ::Subsubsection::Closed:: *)
(*BuildSMDR*)


Options @ BuildSMDR = {"CCompiler" -> Automatic, "CCompilerFlags" -> Automatic, "ArchiveCommand" -> Automatic};


SyntaxInformation @ BuildSMDR = {"ArgumentsPattern" -> {OptionsPattern[]}, "OptionNames" -> Keys[Options @ BuildSMDR]};


BuildSMDR::compx = "SMDR compilation failed. Please run BuildSMDR with adjusted options or specify the configuration option \"SMDRExecutableDirectory\"."


mem : BuildSMDR[opts : OptionsPattern[]] := Monitor[
      
  Module[
    {make},
    
    With[
      {ccompiler = OptionValue @ "CCompiler" /. Automatic :> getCCompiler[]},
    
      With[
        {
         ccompilerflags = OptionValue @ "CCompilerFlags" /. Automatic :> getCCompilerFlags @ ccompiler,
         ar = OptionValue @ "ArchiveCommand"
        },
   
        make = RunProcess[
          {
           "make",
           "-C", 
           RunningCouplings`SMDR`$SourceDirectory, 
           "CC=" <> ccompiler, 
           "SMDR_OPT=" <> ccompilerflags,
           If[StringQ @ ar, "AR=" <> ar, Nothing],
           "INSTALL_INCS=" <> RunningCouplings`SMDR`$IncludeDirectory,
           "INSTALL_LIBS=" <> RunningCouplings`SMDR`$LibraryDirectory
          }
        ];
    
        If[
          make @ "ExitCode" != 0,
      
          RunProcess @ {"make", "-C", RunningCouplings`SMDR`$SourceDirectory, "clean"};
          Failure[
            "CompilationError",
            Association[
              "Message" -> BuildSMDR::compx,
              "CCompiler" -> ccompiler,
              "CCompilerFlags" -> ccompilerflags,
              "ArchiveCommand" -> ar,
              "MakeCommand" -> StringRiffle[
                {
                 "make", "-C", 
                 RunningCouplings`SMDR`$SourceDirectory, 
                 "CC=" <> ccompiler, 
                 "SMDR_OPT=" <> ccompilerflags,
                 If[StringQ @ ar, "AR=" <> ar, Nothing],
                 "INSTALL_INCS=" <> RunningCouplings`SMDR`$IncludeDirectory,
                 "INSTALL_LIBS=" <> RunningCouplings`SMDR`$LibraryDirectory
                }
              ],
              "StandardError" -> "\n" <> make @ "StandardError"
            ]
          ] // Throw
        ];
    
        If[! DirectoryQ @ RunningCouplings`SMDR`$ExecutableDirectory, CreateDirectory @ RunningCouplings`SMDR`$ExecutableDirectory];
        If[! DirectoryQ @ RunningCouplings`SMDR`$IncludeDirectory, CreateDirectory @ RunningCouplings`SMDR`$IncludeDirectory];
        If[! DirectoryQ @ RunningCouplings`SMDR`$LibraryDirectory, CreateDirectory @ RunningCouplings`SMDR`$LibraryDirectory];
    
        RunProcess[{"make", "-C", RunningCouplings`SMDR`$SourceDirectory, "install"}];
    
        RenameFile[#, FileNameJoin @ {RunningCouplings`SMDR`$ExecutableDirectory, FileNameTake @ #}] & /@ Select[
          FileNames[{"calc_*", "fig_*"}, RunningCouplings`SMDR`$SourceDirectory],
          ! DirectoryQ @ # &
        ];
    
        mem = Success[
          "CompilationComplete", 
          Association[
            "Message" -> "SMDR successfully built.",
            "CCompiler" -> ccompiler,
            "CCompilerFlags" -> ccompilerflags,
            "ArchiveCommand" -> ar
          ]
        ]
    
      ]
    ]
    
  ] // Catch,
  
  Row @ {"Building SMDR ", ProgressIndicator[Appearance -> "Ellipsis"]}
  
]


(* ::Subsection:: *)
(*initialization*)


(* ::Subsubsection::Closed:: *)
(*check for SMDR*)


smdrAvailableQ := FileExistsQ[FileNameJoin @ {$smdrExecutableDirectory, "calc_fit"}];


(* ::Subsubsection::Closed:: *)
(*initial compilation*)


If[! smdrAvailableQ, Echo[BuildSMDR[], "RunningCouplings initialization: \n"]];


(* ::Subsection:: *)
(*calc_fit wrapper*)


(* ::Subsubsection::Closed:: *)
(*SMDR input central values and uncertainties*)


mem : smdrInputParameters[] := mem = With[
  {
   file = Import[FileNameJoin @ {RunningCouplings`SMDR`$SourceDirectory, "src", "smdr_pdg.h"}, "Text"]
  },

  AssociationThread[
  
    {
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
    },
      
    Map[
      With[
        {
         centralvalue = StringCases[
           file, 
            # <> "_EXPT" ~~ (WhitespaceCharacter ...) ~~ "=" ~~ Shortest[val__] ~~ ";" :> toPreciseNumber[StringTrim @ val]
          ] // First
        },
        
        Around[
          centralvalue,
          Flatten[
            {
             StringCases[
               file,
               # <> "_EXPT_UNC" ~~ (WhitespaceCharacter ...) ~~ "=" ~~ Shortest[val__] ~~ ";" :> roundTo[ToExpression[StringTrim @ val], centralvalue]
             ],
             StringCases[
               file, 
               # <> "_EXPT_UNC_lo" ~~ (WhitespaceCharacter ...) ~~ "=" ~~ Shortest[val__] ~~ ";" :> roundTo[ToExpression[StringTrim @ val], centralvalue]
             ],
             StringCases[
               file, 
               # <> "_EXPT_UNC_hi" ~~ (WhitespaceCharacter ...) ~~ "=" ~~ Shortest[val__] ~~ ";" :> roundTo[ToExpression[StringTrim @ val], centralvalue]
             ]
            } 
          ] /. {x_} :> x
        ] 
      ] &,
      {
       "SMDR_Mt",
       "SMDR_Mh",
       "SMDR_MZ",
       "SMDR_alphaS_MZ",
       "SMDR_alpha",
       "SMDR_Delta_alpha_had_5_MZ",
       "SMDR_GFermi",
       "SMDR_mbmb",
       "SMDR_mcmc",
       "SMDR_ms_2GeV",
       "SMDR_md_2GeV",
       "SMDR_mu_2GeV",
       "SMDR_Mtau",
       "SMDR_Mmuon",
       "SMDR_Melectron"
      }
    ]
    
  ]

]


mem : smdrInputParameters["Values"] := mem = #["Value"] & /@ smdrInputParameters[]


mem : smdrInputParameters["Ranges"] := mem = MinMax[# @ "Interval"] & /@ smdrInputParameters[]


(* ::Subsubsection:: *)
(*SMDRCalcFit*)


(* ::Subsubsubsection::Closed:: *)
(*iSMDRCalcFit*)


Options @ iSMDRCalcFit = {"InputFile" -> Automatic, "OutputFile" -> Automatic, "ErrorTolerance" -> Automatic, "RenormalizationScale" -> Automatic};


iSMDRCalcFit[opts : OptionsPattern[]] := With[
  {
   input = OptionValue @ "InputFile",
   output = OptionValue @ "OutputFile",
   error = OptionValue @ "ErrorTolerance",
   q = OptionValue @ "RenormalizationScale",
   dir = Directory[]
  },
  
  Module[{calcfit},
  
    SetDirectory[DirectoryName @ input];

    calcfit = RunProcess[
      {
       FileNameJoin @ {$smdrExecutableDirectory, "calc_fit"},
       If[StringQ @ error, {"-e", ToString[error, FortranForm]}, Nothing],
       If[NumericQ @ q, {"-Q", ToString[q, FortranForm]}, Nothing],
       If[StringQ @ input, {"-i", FileNameJoin @ {".", FileNameTake @ input}}, Nothing],
       If[StringQ @ output, {"-o", FileNameJoin @ {".", relativeFilePath[DirectoryName @ input, output]}}, Nothing]
      } // Flatten
    ];
    
    SetDirectory @ dir;
    
    If[
      calcfit @ "ExitCode" != 0, 
      Failure[
        "ExternalProgramError", 
        Association[
          "ExternalCommand" -> StringRiffle[
            {
             FileNameJoin @ {$smdrExecutableDirectory, "calc_fit"},
             If[StringQ @ input, {"-i", input}, Nothing],
             If[StringQ @ output, {"-o", output}, Nothing],
             If[StringQ @ error, {"-e", error}, Nothing],
             If[NumericQ @ q, {"-Q", ToString[q, FortranForm]}, Nothing]
            } // Flatten
          ],
          "StandardOutput" -> "\n" <> calcfit @ "StandardOutput"
        ]
      ] // Throw
    ];
    
    calcfit = If[StringQ @ output, Import[output, "Text"], calcfit @ "StandardOutput"];
    
    Association[
      "RenormalizationScheme" -> "MSbar",
      "RenormalizationScale" -> StringCases[calcfit, "Q = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "g" -> StringCases[calcfit, "g = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "gp" -> StringCases[calcfit, "gp = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "g3" -> StringCases[calcfit, "g3 = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yt" -> StringCases[calcfit, "yt = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yc" -> StringCases[calcfit, "yc = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yu" -> StringCases[calcfit, "yu = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yb" -> StringCases[calcfit, "yb = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "ys" -> StringCases[calcfit, "ys = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yd" -> StringCases[calcfit, "yd = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "ytau" -> StringCases[calcfit, "ytau = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "ymu" -> StringCases[calcfit, "ymu = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "ye" -> StringCases[calcfit, "ye = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "lambda" -> StringCases[calcfit, "lambda = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "m2" -> StringCases[calcfit, "m2 = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "v" -> StringCases[calcfit, "VEV = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "Lambda" -> StringCases[calcfit, "Lambda = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]]
    ]
    
  ] // Catch
  
]


(* ::Subsubsubsection::Closed:: *)
(*SMDRCalcFit*)


Options @ SMDRCalcFit = Options @ iSMDRCalcFit;


SyntaxInformation @ SMDRCalcFit = {"ArgumentsPattern" -> {_., OptionsPattern[]}, "OptionNames" -> Keys @ Options @ SMDRCalcFit}


SMDRCalcFit[opts : OptionsPattern[]] := If[
  smdrAvailableQ, 
  iSMDRCalcFit[opts], 
  With[{make = BuildSMDR[]}, If[smdrAvailableQ, iSMDRCalcFit[opts], make]]
]


SMDRCalcFit[
  args_Association /; ContainsAll[Keys @ args, Keys @ smdrInputParameters[]], 
  opts : OptionsPattern[]
] := With[
  {
   inputfile = FileTemplateApply[
     FileTemplate[
       FileNameJoin @ {RunningCouplings`Developer`$PackageDirectory, "include", "SM", "calc_fit-input.template"},
       InsertionFunction -> (ToString[DecimalForm @ #] &)
     ], 
     args
   ]
  },
 
  With[
    {output = SMDRCalcFit[Sequence @@ (Normal @ <|opts, "InputFile" -> inputfile|>)]},
    DeleteFile @ inputfile;
    output
  ]
  
]


(* ::Subsection:: *)
(*calc_RGrun wrapper*)


(* ::Subsubsection::Closed:: *)
(*input parameter names*)


$calcRGrunInputParameterNames = {"RenormalizationScale", "g", "gp", "g3", "g", "yt", "yc", "yu", "yb", "ys", "yd", "ytau", "ymu", "ye", "lambda", "v", "m2", "Lambda"};


(* ::Subsubsection:: *)
(*SMDRCalcRGrun*)


(* ::Subsubsubsection::Closed:: *)
(*iSMDRCalcRGrun*)


Options @ iSMDRCalcRGrun = {"InputFile" -> Automatic, "LoopOrder" -> Automatic};


iSMDRCalcRGrun[scale_, opts : OptionsPattern[]] := With[
  {
   input = OptionValue @ "InputFile",
   loops = OptionValue @ "LoopOrder",
   dir = Directory[]
  },
  
  Module[{calcRGrun},
  
    SetDirectory[DirectoryName @ input];

    calcRGrun = RunProcess[
      {
       FileNameJoin @ {$smdrExecutableDirectory, "calc_RGrun"},
       ToString[scale, FortranForm],
       If[StringQ @ input, {"-i", FileNameJoin @ {".", FileNameTake @ input}}, Nothing],
       If[IntegerQ @ loops && 1 <= loops <= 5, {"-l", ToString @ loops}, Nothing]
      } // Flatten
    ];
    
    SetDirectory @ dir;
    
    If[
      calcRGrun @ "ExitCode" != 0, 
      Failure[
        "ExternalProgramError", 
        Association[
          "ExternalCommand" -> StringRiffle[
            {
             FileNameJoin @ {$smdrExecutableDirectory, "calc_RGrun"},
             scale,
             If[StringQ @ input, {"-i", input}, Nothing],
             If[IntegerQ @ loops && 1 <= loops <= 5, {"-l", loops}, Nothing]
            } // Flatten
          ],
          "StandardError" -> "\n" <> calcRGrun @ "StandardError"
        ]
      ] // Throw
    ];
    
    calcRGrun = StringDelete[calcRGrun @ "StandardOutput", StartOfString ~~ __ ~~ "FINAL MODEL PARAMETERS:"];
    
    Association[
      "RenormalizationScheme" -> "MSbar",
      "RenormalizationScale" -> StringCases[calcRGrun, "Q = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "g" -> StringCases[calcRGrun, "g = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "gp" -> StringCases[calcRGrun, "gp = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "g3" -> StringCases[calcRGrun, "g3 = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yt" -> StringCases[calcRGrun, "yt = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yc" -> StringCases[calcRGrun, "yc = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yu" -> StringCases[calcRGrun, "yu = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yb" -> StringCases[calcRGrun, "yb = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "ys" -> StringCases[calcRGrun, "ys = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "yd" -> StringCases[calcRGrun, "yd = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "ytau" -> StringCases[calcRGrun, "ytau = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "ymu" -> StringCases[calcRGrun, "ymu = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "ye" -> StringCases[calcRGrun, "ye = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "lambda" -> StringCases[calcRGrun, "lambda = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "m2" -> StringCases[calcRGrun, "m2 = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "v" -> StringCases[calcRGrun, "VEV = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]],
      "Lambda" -> StringCases[calcRGrun, "Lambda = " ~~ val : fortranNumberString :> toPreciseNumber @ val][[1]]
    ]
    
  ] // Catch
  
]


(* ::Subsubsubsection::Closed:: *)
(*SMDRCalcRGrun*)


Options @ SMDRCalcRGrun = Options @ iSMDRCalcRGrun;


SyntaxInformation @ SMDRCalcRGrun = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> Keys @ Options @ SMDRCalcRGrun}


SMDRCalcRGrun[scale_? NumericQ, opts : OptionsPattern[]] := If[
  smdrAvailableQ, 
  iSMDRCalcRGrun[scale, opts], 
  With[{make = BuildSMDR[]}, If[smdrAvailableQ, iSMDRCalcRGrun[scale, opts], make]]
]


SMDRCalcRGrun[
  scale_? NumericQ, 
  args_Association /; ContainsAll[Keys @ args, $calcRGrunInputParameterNames], 
  opts : OptionsPattern[]
] := With[
  {
   inputfile = FileTemplateApply[
     FileTemplate[
       FileNameJoin @ {RunningCouplings`Developer`$PackageDirectory, "include", "SM", "calc_RGrun-input.template"},
       InsertionFunction -> (ToString[DecimalForm @ #] &)
     ], 
     args
   ]
  },
 
  With[
    {output = SMDRCalcRGrun[scale, Sequence @@ (Normal @ <|opts, "InputFile" -> inputfile|>)]},
    DeleteFile @ inputfile;
    output
  ]
  
]


(* ::Section:: *)
(*mixing*)


(* ::Subsection:: *)
(*iMixingMatrix*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ iMixingMatrix = {
  "MajoranaPhaseConvention" -> Automatic, 
  "NeutrinoMassHierarchy" -> Automatic,
  "NeutrinoMassOperator" -> Automatic
};


(* ::Subsubsection::Closed:: *)
(*general*)


iMixingMatrix[{spec_String, "StandardParametrization"}, {\[Theta]12_, \[Theta]23_, \[Theta]13_, \[Delta]_}, opts : OptionsPattern[]] := With[
  {
   s12 = Sin @ \[Theta]12,
   c12 = Sqrt[1 - (Sin[\[Theta]12]^2)],
   s23 = Sin @ \[Theta]23,
   c23 = Sqrt[1 - (Sin[\[Theta]23]^2)],
   s13 = Sin @ \[Theta]13,
   c13 = Sqrt[1 - (Sin[\[Theta]13]^2)]
  },
 
  {
   {c12 c13, s12 c13, s13 Exp[-I \[Delta]]},
   {-s12 c23 - c12 s23 s13 Exp[I \[Delta]], c12 c23 - s12 s23 s13 Exp[I \[Delta]], s23 c13},
   {s12 s23 - c12 c23 s13 Exp[I \[Delta]], -c12 s23 - s12 c23 s13 Exp[I \[Delta]], c23 c13}
  }
 
]


(* ::Subsubsection:: *)
(*CKM*)


(* ::Subsubsubsection::Closed:: *)
(*Wolfenstein parametrization*)


iMixingMatrix[{"CKM", "WolfensteinParametrization"}, {\[Lambda]_, A_, \[Rho]bar_, \[Eta]bar_}, opts : OptionsPattern[]] := {
 {
  1 - (\[Lambda]^2) / 2 - (\[Lambda]^4) / 8, 
  \[Lambda], 
  A (\[Lambda]^3) (\[Rho]bar - I \[Eta]bar)
 },
 {
  -\[Lambda] + (A^2) (\[Lambda]^5) (1 - 2(\[Rho]bar + I \[Eta]bar)) / 2, 
  1 - (\[Lambda]^2) / 2 - (1 + 4 (A^2)) (\[Lambda]^4) / 8, 
  A (\[Lambda]^2)
 },
 {
  A (\[Lambda]^3) (1 - (\[Rho]bar + I \[Eta]bar)), 
  - A (\[Lambda]^2) + A (\[Lambda]^4) (1 - 2 (\[Rho]bar + I \[Eta]bar)) / 2, 
  1 - (A^2) (\[Lambda]^4) / 2
 }
}


(* ::Subsubsubsection::Closed:: *)
(*from parameter association*)


iMixingMatrix["CKM", params_Association] := With[
  {parametrization = Lookup[params, "Parametrization", "Wolfenstein"]},
  
  Switch[
    parametrization,
    
    "Standard",
    iMixingMatrix[{"CKM", "StandardParametrization"}, Lookup[params, {"Theta12", "Theta23", "Theta13", "Delta"}, 0]],
    
    "Wolfenstein",
    iMixingMatrix[{"CKM", "WolfensteinParametrization"}, Lookup[params, {"Lambda", "A", "RhoBar", "EtaBar"}, 0]]
  ]
  
]


(* ::Subsubsubsection::Closed:: *)
(*from mass matrices*)


mem : iMixingMatrix[
  "CKM", 
  {mu_? numericMatrixQ, md_? numericMatrixQ}, 
  opts : OptionsPattern[]
] := mem = With[
  {
   rotation = Switch[
     $config @ "DoubletSingularVectors",
     
     Left | "Left",
     Identity,
     
     Right | "Right",
     ConjugateTranspose
   ]
  },
  
  With[
    {
     vu = orderColumns[SingularValueDecomposition[rotation @ mu, Tolerance -> 0][[1]], Reverse],
     vd = orderColumns[SingularValueDecomposition[rotation @ md, Tolerance -> 0][[1]], Reverse]
    },
    
    ConjugateTranspose[vu].vd
    
  ]
]


(* ::Subsubsection:: *)
(*PMNS*)


(* ::Subsubsubsection::Closed:: *)
(*standard parametrization with Majorana phase*)


iMixingMatrix[{"PMNS", "StandardParametrization"}, {\[Theta]12_, \[Theta]23_, \[Theta]13_, \[Delta]_, \[Phi]1_, \[Phi]2_}, opts : OptionsPattern[]] := With[
  {
   p = Switch[
     OptionValue @ "MajoranaPhaseConvention" /. Automatic -> "PDG",
     
     "PDG",
     DiagonalMatrix[{1, Exp[I \[Phi]1 / 2], Exp[I \[Phi]2 / 2]}],
     
     (* arxiv:hep-ph/0608111 *)
     "MohapatraRodejohann", 
     DiagonalMatrix[{1, Exp[I \[Phi]1], Exp[I (\[Phi]2 + \[Delta])]}],
     
     _,
     IdentityMatrix @ 3
   ]
  },

  iMixingMatrix[{"PMNS", "StandardParametrization"}, {\[Theta]12, \[Theta]23, \[Theta]13, \[Delta]}, opts].p

]


(* ::Subsubsubsection::Closed:: *)
(*from parameter association*)


iMixingMatrix["PMNS", params_Association] := With[
  {
   parametrization = Lookup[params, "Parametrization", "Standard"],
   convention = Lookup[params, "MajoranaPhaseConvention", "PDG"],
   phases = Switch[
     Lookup[params, "MajoranaPhaseConvention", "PDG"],
     
     "PDG",
     {"Alpha21", "Alpha31"},
     
     "MohapatraRodejohann",
     {"Alpha", "Beta"}
   ]
  },
  
  Switch[
    parametrization,
    
    "Standard",
    iMixingMatrix[{"PMNS", "StandardParametrization"}, Lookup[params, {"Theta12", "Theta23", "Theta13", "Delta", Sequence @@ phases}, 0], "MajoranaPhaseConvention" -> convention]
  ]
  
]


(* ::Subsubsubsection::Closed:: *)
(*from mass matrices*)


mem : iMixingMatrix[
  "PMNS", 
  {me_? MatrixQ, mv_? SymmetricMatrixQ}, 
  opts : OptionsPattern[]
] := mem = With[
  {
   order = Switch[
     OptionValue @ "NeutrinoMassHierarchy" /. Automatic -> "Normal",
     
     "Normal",
     {3, 2, 1},
     
     "Inverted",
     {2, 1, 3},
     
     Alternatives @@ (Permutations @ {1, 2, 3}),
     Reverse[OptionValue @ "NeutrinoMassHierarchy"]
   ],
   rotation = Switch[
     $config @ "DoubletSingularVectors",
     
     Left | "Left",
     Identity,
     
     Right | "Right",
     ConjugateTranspose
   ]
  },

  With[ 
    {
     ve = orderColumns[SingularValueDecomposition[rotation @ me, Tolerance -> 0][[1]], Reverse],
     vv = orderColumns[TakagiFactorization[mv, Tolerance -> 0][[1]], order] // Conjugate
    },
  
    ConjugateTranspose[ve].vv
 
  ] 
]


mem : iMixingMatrix[
  "PMNS", 
  {me_? numericMatrixQ, mv_? numericMatrixQ}, 
  opts : OptionsPattern[]
] := mem = With[
  {
   order = Switch[
     OptionValue @ "NeutrinoMassHierarchy" /. Automatic -> "Normal",
     
     "Normal",
     {3, 2, 1},
     
     "Inverted",
     {2, 1, 3},
     
     Alternatives @@ (Permutations @ {1, 2, 3}),
     Reverse[OptionValue @ "NeutrinoMassHierarchy"]
   ],
   rotation = Switch[
     $config @ "DoubletSingularVectors",
     
     Left | "Left",
     Identity,
     
     Right | "Right",
     ConjugateTranspose
   ]
  },

  With[ 
    {
     ve = orderColumns[SingularValueDecomposition[rotation @ me, Tolerance -> 0][[1]], Reverse],
     vv = Switch[
       OptionValue @ "NeutrinoMassOperator" /. Automatic -> "Majorana",
       
       "Majorana",
       orderColumns[TakagiFactorization[Symmetrize @ mv, Tolerance -> 0][[1]], order] // Conjugate,
       
       "Dirac",
       orderColumns[SingularValueDecomposition[rotation @ mv, Tolerance -> 0][[1]], order]
     ]
    },
  
    ConjugateTranspose[ve].vv
 
  ] 
]


(* ::Subsection:: *)
(*MixingMatrix*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ MixingMatrix = Normal[<|"Parametrization" -> Automatic, Options @ iMixingMatrix|> // KeySort];


SyntaxInformation @ MixingMatrix = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> Keys @ Options @ MixingMatrix};


(* ::Subsubsection::Closed:: *)
(*general*)


MixingMatrix[
  spec_String, 
  params_List, 
  opts : OptionsPattern[]
] := NAroundApply[
  iMixingMatrix[
    {spec, (OptionValue @ "Parametrization" /. Automatic -> "Standard") <> "Parametrization"}, 
    {##}, 
    FilterRules[{opts}, Options @ iMixingMatrix]
  ] &,
  params
]


MixingMatrix[
  spec_String, 
  params_Association
] := NAroundApply[iMixingMatrix[spec, AssociationThread[Keys @ params, {##}]] &, Values @ params]


MixingMatrix[
  spec_String, 
  m : {m1_? MatrixQ, m2_? MatrixQ} /; AllTrue[Flatten @ {Dimensions @ m1, Dimensions @ m1}, (# == 3) &], 
  opts : OptionsPattern[]
] := NAroundApply[
  iMixingMatrix[spec, ArrayReshape[{##}, {2, 3, 3}], FilterRules[{opts}, Options @ iMixingMatrix]] &,
  Flatten @ m
]


(* ::Subsection:: *)
(*iMixingParameters*)


(* ::Subsubsection::Closed:: *)
(*general*)


iMixingParameters[{spec_String, "StandardParametrization"}, u_? MatrixQ] := With[
  {
   \[Theta]12 = ArcSin[Abs[u[[1, 2]]] / Sqrt[1 - (Abs[u[[1, 3]]]^2)]],
   \[Theta]23 = ArcSin[Abs[u[[2, 3]]] / Sqrt[1 - (Abs[u[[1, 3]]]^2)]],
   \[Theta]13 = ArcSin @ Abs[u[[1, 3]]]
  },

  Association[
    "Parametrization" -> "Standard",
    "Theta12" -> \[Theta]12,
    "Theta23" -> \[Theta]23,
    "Theta13" -> \[Theta]13,
    (*
    "SineSquaredTheta12" \[Rule] (Abs[u\[LeftDoubleBracket]1, 2\[RightDoubleBracket]]^2) / (1 - Abs[u\[LeftDoubleBracket]1, 3\[RightDoubleBracket]]^2),
    "SineSquaredTheta23" \[Rule] (Abs[u\[LeftDoubleBracket]2, 3\[RightDoubleBracket]]^2) / (1 - Abs[u\[LeftDoubleBracket]1, 3\[RightDoubleBracket]]^2),
    "SineSquaredTheta13" \[Rule] (Abs[u\[LeftDoubleBracket]1, 3\[RightDoubleBracket]]^2),
    *)
    
    With[
      {
       s\[Delta] = 8 Im[u[[2, 3]] Conjugate[u[[1, 3]]] u[[1, 2]] Conjugate[u[[2, 2]]]] / (Cos[\[Theta]13] Sin[2 \[Theta]12] Sin[2 \[Theta]23] Sin[2 \[Theta]13]),
       c\[Delta] = ((Cos[\[Theta]12]^2) (Cos[\[Theta]23]^2) + (Sin[\[Theta]12]^2) (Sin[\[Theta]23]^2) (Sin[\[Theta]13]^2) - (Abs[u[[2, 2]]]^2)) / (2 Cos[\[Theta]12] Cos[\[Theta]23] Sin[\[Theta]12] Sin[\[Theta]23] Sin[\[Theta]13])
      },
      
    "Delta" -> NAroundApply[Mod[ArcTan[##], 2\[Pi]] &, {c\[Delta], s\[Delta]}]
    
    ]
  ]
  
]


(* ::Subsubsection::Closed:: *)
(*CKM*)


iMixingParameters[{"CKM", "WolfensteinParametrization"}, v_? MatrixQ] := Association[
  "Parametrization" -> "Wolfenstein",
  "Lambda" -> Abs[v[[1, 2]]] /  Sqrt[Abs[v[[1, 1]]]^2 + Abs[v[[1, 2]]]^2],
  "A" -> (Abs[v[[2, 3]]] / (Abs[v[[1, 2]]]^2)) Sqrt[Abs[v[[1, 1]]]^2 + Abs[v[[1, 2]]]^2],
  "RhoBar" -> Re[- (v[[1, 1]] Conjugate[v[[1, 3]]]) / (v[[2, 1]] Conjugate[v[[2, 3]]])],
  "EtaBar" -> Im[- (v[[1, 1]] Conjugate[v[[1, 3]]]) / (v[[2, 1]] Conjugate[v[[2, 3]]])]
] 


iMixingParameters["CKM", v_? MatrixQ]["JarlskogInvariant"] := Im[v[[2, 2]] Conjugate[v[[1, 2]]] v[[1, 1]] Conjugate[v[[2, 1]]]]


(* ::Subsubsection::Closed:: *)
(*PMNS*)


(* ::Input:: *)
(*iMixingParameters[{"PMNS", "StandardParametrization"}, u_? MatrixQ] := Association[*)
(*    *)
(*  iMixingParameters[{"PMNS", "StandardParametrization"}, u],*)
(*  *)
(*  Switch[*)
(*    OptionValue @ "MajoranaPhaseConvention" /. Automatic -> "PDG",*)
(*  *)
(*    "PDG",*)
(*    {*)
(*     "MajoranaPhaseConvention" -> "PDG",*)
(*     "Alpha12" -> aroundMod[(Arg[u[[1, 2]]] /. Derivative[1][Arg] -> 0) + 2\[Pi], 2\[Pi]],*)
(*     "Alpha13" -> aroundMod[(Arg[u[[2, 3]]] /. Derivative[1][Arg] -> 0) + 2\[Pi], 2\[Pi]]*)
(*    },*)
(*    *)
(*    "MohapatraRodejohann",*)
(*    {*)
(*     "MajoranaPhaseConvention" -> "MohapatraRodejohann",*)
(*     "Alpha" -> aroundMod[(Arg[u[[1, 2]]] /. Derivative[1][Arg] -> 0) + 2\[Pi], 2\[Pi]] / 2,*)
(*     "Beta" -> aroundMod[(Arg[u[[1, 3]]] /. Derivative[1][Arg] -> 0) + 2\[Pi], 2\[Pi]] / 2*)
(*    },*)
(*      *)
(*    _,*)
(*    Nothing*)
(*  ]*)
(*    *)
(*]*)


iMixingParameters["PMNS", u_? MatrixQ]["JarlskogInvariant"] := Im[u[[2, 3]] Conjugate[u[[1, 3]]] u[[1, 2]] Conjugate[u[[2, 2]]]]


(* ::Subsection:: *)
(*MixingParameters*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ MixingParameters = Normal[<|"Parametrization" -> Automatic, Options @ MixingMatrix|> // KeySort];


SyntaxInformation @ MixingParameters = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}, "OptionNames" -> Keys @ Options @ MixingParameters};


(* ::Subsubsection::Closed:: *)
(*general*)


MixingParameters[
  spec_String, 
  u_? MatrixQ, 
  opts : OptionsPattern[]
] := iMixingParameters[{spec, (OptionValue @ "Parametrization" /. Automatic -> "Standard") <> "Parametrization"}, u]


MixingParameters[
  spec_String, 
  m : {m1_? MatrixQ, m2_? MatrixQ} /; AllTrue[Flatten @ {Dimensions @ m1, Dimensions @ m1}, (# == 3) &], 
  opts : OptionsPattern[]
] := MixingParameters[spec, MixingMatrix[spec, m, FilterRules[{opts}, Options @ MixingMatrix]], opts]


MixingParameters[spec_String, params_Association, opts : OptionsPattern[]] := If[
  SameQ[params @ "Parametrization", OptionValue @ "Parametrization" /. Automatic -> "Standard"],
  params,
  MixingParameters[spec, MixingMatrix[spec, params], opts]
]


(* ::Section:: *)
(*renormalization*)


(* ::Subsection:: *)
(*flavor mixing data*)


(* ::Subsubsection::Closed:: *)
(*import*)


Import[FileNameJoin @ {RunningCouplings`Developer`$PackageDirectory, "include", "SM", "flavormixing.wl"}, "Package"];


(* ::Subsubsection::Closed:: *)
(*CKM*)


$ckm @ "WolfensteinLambda" = Around[$ckm["WolfensteinLambda", "CentralValue"], Abs @ $ckm["WolfensteinLambda", "Uncertainty"]];
$ckm @ "WolfensteinA" = Around[$ckm["WolfensteinA", "CentralValue"], Abs @ $ckm["WolfensteinA", "Uncertainty"]];
$ckm @ "WolfensteinRhoBar" = Around[$ckm["WolfensteinRhoBar", "CentralValue"], Abs @ $ckm["WolfensteinRhoBar", "Uncertainty"]];
$ckm @ "WolfensteinEtaBar" = Around[$ckm["WolfensteinEtaBar", "CentralValue"], Abs @ $ckm["WolfensteinEtaBar", "Uncertainty"]];


$ckm["WolfensteinParametrization"] = MixingMatrix["CKM", "WolfensteinParametrization"][
  {
   $ckm @ "WolfensteinLambda", 
   $ckm @ "WolfensteinA", 
   $ckm @ "WolfensteinRhoBar", 
   $ckm @ "WolfensteinEtaBar"
  }
];


(* ::Subsubsection::Closed:: *)
(*PMNS*)


Options @ $pmns = {"NeutrinoMassHierarchy" -> Automatic};


Scan[
  (
   $pmns[#, "SineSquaredTheta12"] = Around[$pmns[#, "SineSquaredTheta12", "CentralValue"], Abs @ $pmns[#, "SineSquaredTheta12", "Uncertainty"]];
   $pmns[#, "SineSquaredTheta23"] = Around[$pmns[#, "SineSquaredTheta23", "CentralValue"], Abs @ $pmns[#, "SineSquaredTheta23", "Uncertainty"]];
   $pmns[#, "SineSquaredTheta13"] = Around[$pmns[#, "SineSquaredTheta13", "CentralValue"], Abs @ $pmns[#, "SineSquaredTheta13", "Uncertainty"]];
   $pmns[#, "Theta12"] = Around[$pmns[#, "Theta12", "CentralValue"], Abs @ $pmns[#, "Theta12", "Uncertainty"]];
   $pmns[#, "Theta23"] = Around[$pmns[#, "Theta23", "CentralValue"], Abs @ $pmns[#, "Theta23", "Uncertainty"]];
   $pmns[#, "Theta13"] = Around[$pmns[#, "Theta13", "CentralValue"], Abs @ $pmns[#, "Theta13", "Uncertainty"]];
   $pmns[#, "Delta"] = Around[$pmns[#, "Delta", "CentralValue"], Abs @ $pmns[#, "Delta", "Uncertainty"]];
  ) &,
  {"NormalHierarchy", "InvertedHierarchy"}
]


Scan[
  (
   $pmns[#, "StandardParametrization"] = MixingMatrix["PMNS", "StandardParametrization"][
     {
      $pmns[#, "Theta12"],
      $pmns[#, "Theta23"],
      $pmns[#, "Theta13"],
      $pmns[#, "Delta"]
     }
   ]
  ) &,
  {"NormalHierarchy", "InvertedHierarchy"}
]


$pmns[para_String, opts : OptionsPattern[]] := $pmns[(OptionValue["NeutrinoMassHierarchy"] /. Automatic -> "Normal") <> "Hierarchy", para]


(* ::Subsection::Closed:: *)
(*neutrino mass data*)


(* ::Subsubsection::Closed:: *)
(*import*)


Import[FileNameJoin @ {RunningCouplings`Developer`$PackageDirectory, "include", "SM", "neutrino.wl"}, "Package"];


$neutrino["NormalHierarchy", "DeltaMassSquared21"] = Around[
  $neutrino["NormalHierarchy", "DeltaMassSquared21", "CentralValue"], 
  $neutrino["NormalHierarchy", "DeltaMassSquared21", "Uncertainty"] // Abs
];

$neutrino["NormalHierarchy", "DeltaMassSquared31"] = Around[
  $neutrino["NormalHierarchy", "DeltaMassSquared31", "CentralValue"], 
  $neutrino["NormalHierarchy", "DeltaMassSquared31", "Uncertainty"] // Abs
];

$neutrino["InvertedHierarchy", "DeltaMassSquared21"] = Around[
  $neutrino["InvertedHierarchy", "DeltaMassSquared21", "CentralValue"], 
  $neutrino["InvertedHierarchy", "DeltaMassSquared21", "Uncertainty"] // Abs
];

$neutrino["InvertedHierarchy", "DeltaMassSquared32"] = Around[
  $neutrino["InvertedHierarchy", "DeltaMassSquared32", "CentralValue"], 
  $neutrino["InvertedHierarchy", "DeltaMassSquared32", "Uncertainty"] // Abs
];


(* ::Subsubsection::Closed:: *)
(*neutrinoMasses*)


Options @ neutrinoMasses = {"LightestNeutrinoMass" -> Automatic, "NeutrinoMassHierarchy" -> Automatic};


neutrinoMasses[opts : OptionsPattern[]] := With[
  {
   m0 = OptionValue @ "LightestNeutrinoMass" /. Automatic -> 0.0,
   hierarchy = OptionValue @ "NeutrinoMassHierarchy" /. Automatic -> "Normal"
  },

  Switch[
    hierarchy,
  
   "Inverted",
    Sqrt[
      {
       m0^2 - ($neutrino["NormalHierarchy", "DeltaMassSquared32"]) + ($neutrino["NormalHierarchy", "DeltaMassSquared21"] 10^(-18)),
       m0^2 - ($neutrino["NormalHierarchy", "DeltaMassSquared32"] 10^(-18)),
       m0^2
      }
    ],
  
    "Normal",
    Sqrt[
      {
       m0^2, 
       m0^2 + ($neutrino["NormalHierarchy", "DeltaMassSquared21"] 10^(-18)), 
       m0^2 + ($neutrino["NormalHierarchy", "DeltaMassSquared31"] 10^(-18))
      }
    ]
  ] 

]


(* ::Subsection:: *)
(*beta function data*)


(* ::Subsubsection::Closed:: *)
(*import*)


mem : ibetaFunction[modelname_] := mem = Module[
  {array},
  
  Block[
    {$Context = $context},
    
    array = FileNames[
             "beta-*", 
             FileNameJoin @ {RunningCouplings`Developer`$PackageDirectory, "include", modelname, "rge"}
           ];
  
    array = Join @@ (Import[#, "Package"] & /@ array);
    array = ReplaceRepeated[
      array,
      
      With[
        {
         rotation = Switch[
           $config @ "DoubletSingularVectors",
           
           Left | "Left",
           Identity,
           
           Right | "Right",
           ConjugateTranspose
         ]
        },
         
         {
          MatMul[x__] :> rotation[Dot @ x],
          trace[x__] :> Tr @ rotation[Dot @ x],
          Adj -> ConjugateTranspose,
          Tp -> Transpose,
          conj -> Conjugate
         }
         
      ]
    ];
    array = array /. ((Symbol[#] -> Symbol[#][t]) & /@ array[[All, 1]]);
  
    Association[(#1 -> {##2} Log[10] / ((16 \[Pi]^2)^Range[Length @ {##2}])) & @@@ array]
  
  ]
]


(* ::Subsubsection::Closed:: *)
(*betaFunction*)


betaFunction[modelname_, param_] := Total @ ibetaFunction[modelname][param]


betaFunction[modelname_, param_, n_Integer] := Total[ibetaFunction[modelname][param][[1 ;; n]]]


(* ::Subsection:: *)
(*model matching*)


(* ::Subsubsection:: *)
(*matchModels*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ matchModels = {
  "CKMParametrization" -> Automatic,
  "FlavorMixing" -> Automatic,
  "LightestNeutrinoMass" -> Automatic, 
  "NeutrinoMassOperator" -> Automatic,
  (* "NeutrinoMajoranaPhases" \[Rule] Automatic, *) 
  "NeutrinoMassHierarchy" -> Automatic,
  "PMNSParametrization" -> Automatic,
  "TanBeta" -> Automatic
};


(* ::Subsubsubsection::Closed:: *)
(*SM (SMDR) -> SM (NDSolve) matching*)


matchModels[
  {"SM", "SMDR"} -> {"SM", "NDSolve"}, 
  params_, 
  opts : OptionsPattern[]
] := With[
  {
   operator = OptionValue @ "NeutrinoMassOperator" /. Automatic -> "Majorana",
   mixing = OptionValue @ "FlavorMixing" /. Automatic -> All
   (* phases = Replace[OptionValue @ "NeutrinoMajoranaPhases", Automatic \[Rule] {"Convention" \[Rule] "PDG", "Alpha21" \[Rule] 0, "Alpha31" \[Rule] 0}] *)
  },
  
  
  Association[
  
    KeyTake[params, {"RenormalizationScale"}],
    "g1" -> Sqrt[5 / 3] params["gp"],
    "g2" -> params @ "g",
    KeyTake[params, {"g3", "yu", "yc", "yt", "yd", "ys", "yb", "ye", "ymu", "ytau"}],

    Switch[
      operator,
  
      "Dirac",
      Thread[{"y1", "y2", "y3"} -> (Sqrt[2] / params["v"]) neutrinoMasses @ FilterRules[{opts}, Options @ neutrinoMasses]],
      
      "Majorana",
      Thread[{"kappa1", "kappa2", "kappa3"} -> (2 / (params["v"]^2)) neutrinoMasses @ FilterRules[{opts}, Options @ neutrinoMasses]],
    
      None,
      Nothing
    ],
    
    With[
      {
       ckm = MixingParameters[
         "CKM",
         Association[
           "Parametrization" -> "Wolfenstein",
           "Lambda" -> $ckm @ "WolfensteinLambda",
           "A" -> $ckm @ "WolfensteinA",
           "RhoBar" -> $ckm @ "WolfensteinRhoBar",
           "EtaBar" -> $ckm @ "WolfensteinEtaBar"
         ],
         "Parametrization" -> (OptionValue @ "CKMParametrization" /. Automatic -> "Wolfenstein")
       ],
       pmns = MixingParameters[
         "PMNS",
         Association[
           "Parametrization" -> "Standard", 
           "Theta12" -> $pmns["Theta12", FilterRules[{opts}, Options @ $pmns]],
           "Theta23" -> $pmns["Theta23", FilterRules[{opts}, Options @ $pmns]],
           "Theta13" -> $pmns["Theta13", FilterRules[{opts}, Options @ $pmns]],
           (*
           "SineSquaredTheta12" \[Rule] $pmns["SineSquaredTheta12", FilterRules[{opts}, Options @ $pmns]],
           "SineSquaredTheta23" \[Rule] $pmns["SineSquaredTheta23", FilterRules[{opts}, Options @ $pmns]],
           "SineSquaredTheta13" \[Rule] $pmns["SineSquaredTheta13", FilterRules[{opts}, Options @ $pmns]],
           *)
           "Delta" -> $pmns["Delta", FilterRules[{opts}, Options @ $pmns]]
           (*
           If[
             operator === "Majorana",
             {
              "MajoranaPhaseConvention" \[Rule] ("Convention" /. phases) /. Automatic \[Rule] "PDG",
              KeyDrop[phases, "Convention"] /. Automatic \[Rule] 0
             },
             Nothing
           ]
           *)
         ],
         "Parametrization" -> (OptionValue @ "PMNSParametrization" /. Automatic -> "Standard")
       ]
      },
    
      Switch[
        mixing,
      
        "Quark",
        "CKM" -> ckm,
      
        "Lepton",
        "PMNS" -> pmns,
      
        All | True,
        {
         "CKM" -> ckm,
         "PMNS" -> pmns
        },
      
        None | False,
        Nothing
      ]
      
    ],
    
    "lambda" -> 2 params["lambda"],
    "m2" -> params @ "m2"
    
  ]
  
]


(* ::Subsubsubsection::Closed:: *)
(*SM (NDSolve) -> MSSM (NDSolve) matching*)


matchModels[
  {"SM", "NDSolve"} -> {"MSSM", "NDSolve"}, 
  params_, 
  opts : OptionsPattern[]
] := With[
  {
   operator = OptionValue @ "NeutrinoMassOperator" /. Automatic -> "Majorana",
   mixing = OptionValue @ "FlavorMixing" /. Automatic -> All,
   tan\[Beta] = OptionValue @ "TanBeta" /. Automatic -> 10,
   g1 = params["g1"],
   g2 = params["g2"],
   g3 = params["g3"]
  },
  
  With[
    {
      sin\[Beta] = tan\[Beta] / Sqrt[1 + tan\[Beta]^2],
      cos\[Beta] = 1 / Sqrt[1 + tan\[Beta]^2]
    },
  
    Association[
  
      KeyTake[params, {"RenormalizationScale"}], 
      "SupersymmetryScale" -> params @ "RenormalizationScale",
      "TanBeta" -> tan\[Beta],
    
      (* gauge coupling Overscript[MS, _] \[Rule] Overscript[DR, _] regularization conversion from arXiv:hep-ph/9308222 Eq. 2.4 *)
      "g1" -> g1,
      "g2" -> g2 / (1 - (g1^2) / (48 \[Pi]^2)),
      "g3" -> g3 / (1 - (g3^2) / (32 \[Pi]^2)),
    
      (* Yukawa coupling Overscript[MS, _] \[Rule] Overscript[DR, _] regularization conversion from arXiv:hep-ph/9308222 Eq. 2.8 *)
      "yu" -> (params["yu"] / sin\[Beta]) / (1 - (g1^2) / (1920 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2) + (g3^2) / (12 \[Pi]^2)),
      "yc" -> (params["yc"] / sin\[Beta]) / (1 - (g1^2) / (1920 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2) + (g3^2) / (12 \[Pi]^2)),
      "yt" -> (params["yt"] / sin\[Beta]) / (1 - (g1^2) / (1920 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2) + (g3^2) / (12 \[Pi]^2)),
      "yd" -> (params["yd"] / cos\[Beta]) / (1 - 13 (g1^2) / (1920 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2) + (g3^2) / (12 \[Pi]^2)),
      "ys" -> (params["ys"] / cos\[Beta]) / (1 - 13 (g1^2) / (1920 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2) + (g3^2) / (12 \[Pi]^2)),
      "yb" -> (params["yb"] / cos\[Beta]) / (1 - 13 (g1^2) / (1920 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2) + (g3^2) / (12 \[Pi]^2)),
      "ye" -> (params["ye"] / cos\[Beta]) / (1 + 9 (g1^2) / (640 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2)),
      "ymu" -> (params["ymu"] / cos\[Beta]) / (1 + 9 (g1^2) / (640 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2)),
      "ytau" -> (params["ytau"] / cos\[Beta]) / (1 + 9 (g1^2) / (640 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2)),

      Switch[
        operator,
  
        "Dirac",
        {
         (* Yukawa coupling Overscript[MS, _] \[Rule] Overscript[DR, _] regularization conversion from arXiv:hep-ph/9308222 Eq. 2.8 *)
         "y1" -> (params["y1"] / sin\[Beta]) / (1 - 3 (g1^2) / (640 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2)),
         "y2" -> (params["y2"] / sin\[Beta]) / (1 - 3 (g1^2) / (640 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2)),
         "y3" -> (params["y3"] / sin\[Beta]) / (1 - 3 (g1^2) / (640 \[Pi]^2) - 3 (g2^2) / (128 \[Pi]^2))
        },
      
        "Majorana",
        {
         "kappa1" -> params["kappa1"] / (sin\[Beta]^2),
         "kappa2" -> params["kappa2"] / (sin\[Beta]^2),
         "kappa3" -> params["kappa3"] / (sin\[Beta]^2)
        },
    
        None,
        Nothing
      ],
      
      With[
        {
         ckm = MixingParameters[
           "CKM", 
           params @ "CKM", 
           "Parametrization" -> (OptionValue @ "CKMParametrization" /. Automatic -> "Wolfenstein")
         ],
         pmns = MixingParameters[
           "PMNS", 
           params @ "PMNS", 
           "Parametrization" -> (OptionValue @ "PMNSParametrization" /. Automatic -> "Standard")
         ]
        },
    
        Switch[
          mixing,
      
          "Quark",
          "CKM" -> ckm,
      
          "Lepton",
          "PMNS" -> pmns,
      
          All | True,
          {
           "CKM" -> ckm,
           "PMNS" -> pmns
          },
      
          None | False,
          Nothing
        ]
        
      ]
    
    ]
  
  ]
]


(* ::Subsection:: *)
(*diagonalization and mixing utilities*)


(* ::Subsubsection::Closed:: *)
(*yukawaSingularValueInterpolation*)


yukawaSingularValueInterpolation[y_InterpolatingFunction, ordering_: {1, 2, 3}] := With[
  {
   coords = First[InterpolatingFunctionCoordinates @ y],
   sv = Transpose[Reverse[PadLeft[SingularValueList[#, Tolerance -> 0], 3]] & /@ InterpolatingFunctionValuesOnGrid[y]],
   order = First[InterpolatingFunctionInterpolationOrder @ y],
   position = PositionIndex @ ordering
  },
          
  {
   Interpolation[{coords, sv[[position[1] // First]]} // Transpose, InterpolationOrder -> order],
   Interpolation[{coords, sv[[position[2] // First]]} // Transpose, InterpolationOrder -> order],
   Interpolation[{coords, sv[[position[3] // First]]} // Transpose, InterpolationOrder -> order]
  }
          
]


(* ::Subsubsection::Closed:: *)
(*mixingParametersInterpolation*)


Options @ mixingParametersInterpolation = Normal[<|"TanBeta" -> Automatic, Options @ MixingParameters|> // KeySort];


mixingParametersInterpolation[
 spec_String, 
 {m1_InterpolatingFunction, m2_InterpolatingFunction},
 opts : OptionsPattern[]
] := With[
  {
   coords = First[InterpolatingFunctionCoordinates @ m1],
   matrix = Merge[
     MapThread[
       MixingParameters[spec, {##}, FilterRules[{opts}, Options @ MixingParameters]] &,
       {
        InterpolatingFunctionValuesOnGrid @ m1,
        m2 /@ First[InterpolatingFunctionCoordinates @ m1]
       }
     ],
     Identity
   ] // Normal,
   order = First[InterpolatingFunctionInterpolationOrder @ m1]
  },
          
  Association[
    matrix /. {
      points : {__Real} :> Interpolation[{coords, points} // Transpose, InterpolationOrder -> order],
      strings : {__String} :> First @ strings
    }
  ]
          
]


(* ::Subsection:: *)
(*nRenormalizeModel*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ nRenormalizeModel = Normal[<|"LoopOrder" -> Automatic, Options @ NDSolve|>];


(* ::Subsubsection:: *)
(*SM*)


(* ::Subsubsubsection::Closed:: *)
(*NDSolve*)


mem : nRenormalizeModel[{"SM", "NDSolve"}, args_Association, opts : OptionsPattern[]] := mem = With[
  {
   tmin = 0,
   tin = Log10[args @ "RenormalizationScale"],
   tmax = 18,
   
   loops = OptionValue @ "LoopOrder" /. Automatic -> 2,
   
   dirac = ContainsAll[Keys @ args, {"y1", "y2", "y3"}],
   majorana = ContainsAll[Keys @ args, {"kappa1", "kappa2", "kappa3"}],
   qflavor = KeyExistsQ[args, "CKM"],
   lflavor = KeyExistsQ[args, "PMNS"],
   
   rargs = args /. n_? NumericQ :> toRoundedNumber @ n
  },
  
  Block[
    {$Context = $context},
  
    MapAt[
      ToString,
      NDSolve[
        {
         
         {
          g1'[t] == betaFunction["SM", "g1", loops],
          g2'[t] == betaFunction["SM", "g2", loops],
          g3'[t] == betaFunction["SM", "g3", loops],
          Yu'[t] == betaFunction["SM", "Yu", loops],
          Yd'[t] == betaFunction["SM", "Yd", loops],
          Ye'[t] == betaFunction["SM", "Ye", loops],
          
          If[
            dirac, 
            Yv'[t] == betaFunction["SM", "Yv", loops], 
            Nothing
          ],
          
          If[
            majorana, 
            kappa'[t] == betaFunction["SM", "kappa", loops], 
            Nothing
          ],
          
          lambda'[t] == betaFunction["SM", "lambda", loops],
          m2'[t] == betaFunction["SM", "m2", loops]
         } /. {
           If[dirac, kappa[t] -> ConstantArray[0, {3, 3}], Nothing],
           If[majorana, Yv[t] -> ConstantArray[0, {3, 3}], Nothing]
         },
         
         {
          g1[tin] == rargs @ "g1",
          g2[tin] == rargs @ "g2",
          g3[tin] == rargs @ "g3",
          
          With[
            {
             ckm = If[qflavor, MixingMatrix["CKM", rargs @ "CKM"], IdentityMatrix @ 3],
             pmns = If[lflavor, MixingMatrix["PMNS", rargs @ "PMNS"], IdentityMatrix @ 3]
            },
          
            {
             Switch[
               $config @ "DoubletSingularVectors",
               
               Left | "Left",
               {
                Yu[tin] == ConjugateTranspose[ckm].DiagonalMatrix[Lookup[rargs, {"yu", "yc", "yt"}]],
                If[
                  dirac,
                  Yv[tin] == pmns.DiagonalMatrix[Lookup[rargs, {"y1", "y2", "y3"}]],
                  Nothing
                ]
               },
               
               Right | "Right",
               {
                Yu[tin] == DiagonalMatrix[Lookup[rargs, {"yu", "yc", "yt"}]].ckm,
                If[
                  dirac,
                  Yv[tin] == DiagonalMatrix[Lookup[rargs, {"y1", "y2", "y3"}]].ConjugateTranspose[pmns],
                  Nothing
                ]
               }
             ],
             
             If[
               majorana,
               kappa[tin] == Conjugate[pmns].DiagonalMatrix[Lookup[rargs, {"kappa1", "kappa2", "kappa3"}]].ConjugateTranspose[pmns],
               Nothing
             ]
            }
          ],
          
          Yd[tin] == DiagonalMatrix[Lookup[rargs, {"yd", "ys", "yb"}]],
          Ye[tin] == DiagonalMatrix[Lookup[rargs, {"ye", "ymu", "ytau"}]],
          
          lambda[tin] == rargs @ "lambda",
          m2[tin] == rargs @ "m2"
         }
         
         (*
         If[
           majorana,
           WhenEvent[
             Mod[t, 1], 
             kappa[t] \[Rule] Symmetrize[kappa[t], Symmetric @ {1, 2}],
             "DetectionMethod" \[Rule] "Sign",
             "LocationMethod" \[Rule] "StepEnd",
             "IntegrateEvent" \[Rule] False
           ],
           Nothing
         ]
         *)
          
        } // Flatten,
        
        {
         g1,
         g2,
         g3,
         Yu,
         Yd,
         Ye,
         If[dirac, Yv, Nothing],
         If[majorana, kappa, Nothing],
         lambda,
         m2
        },
        
        {t, tmin, tmax},
        
        FilterRules[{opts}, Options @ NDSolve]
      ] // First,
      {All, 1}
    ] // Association
        
  ]
]


(* ::Subsubsection:: *)
(*MSSM*)


(* ::Subsubsubsection::Closed:: *)
(*NDSolve*)


mem : nRenormalizeModel[{"MSSM", "NDSolve"}, args_Association, opts : OptionsPattern[]] := mem = With[
  {
   tmin = 0,
   tin = Log10[args @ "RenormalizationScale"],
   tmax = 18,
   
   loops = OptionValue @ "LoopOrder" /. Automatic -> 2,
   
   dirac = ContainsAll[Keys @ args, {"y1", "y2", "y3"}],
   majorana = ContainsAll[Keys @ args, {"kappa1", "kappa2", "kappa3"}],
   qflavor = KeyExistsQ[args, "CKM"],
   lflavor = KeyExistsQ[args, "PMNS"],
   
   rargs = args /. n_? NumericQ :> toRoundedNumber @ n
  },
  
  Block[
    {$Context = $context},
  
    MapAt[
      ToString,
      NDSolve[
        {
         
         {
          g1'[t] == betaFunction["MSSM", "g1", loops],
          g2'[t] == betaFunction["MSSM", "g2", loops],
          g3'[t] == betaFunction["MSSM", "g3", loops],
          Yu'[t] == betaFunction["MSSM", "Yu", loops],
          Yd'[t] == betaFunction["MSSM", "Yd", loops],
          Ye'[t] == betaFunction["MSSM", "Ye", loops],
           
          If[
            dirac, 
            Yv'[t] == betaFunction["MSSM", "Yv", loops], 
            Nothing
          ],
          
          If[
            majorana, 
            kappa'[t] == betaFunction["MSSM", "kappa", loops], 
            Nothing
          ]
         } /. {
           If[dirac, kappa[t] -> ConstantArray[0, {3, 3}], Nothing],
           If[majorana, Yv[t] -> ConstantArray[0, {3, 3}], Nothing]
         },
         
         {
          g1[tin] == rargs @ "g1",
          g2[tin] == rargs @ "g2",
          g3[tin] == rargs @ "g3",
          
          With[
            {
             ckm = If[qflavor, MixingMatrix["CKM", rargs @ "CKM"], IdentityMatrix @ 3],
             pmns = If[lflavor, MixingMatrix["PMNS", rargs @ "PMNS"], IdentityMatrix @ 3]
            },
          
            {
             Switch[
               $config @ "DoubletSingularVectors",
               
               Left | "Left",
               {
                Yu[tin] == ConjugateTranspose[ckm].DiagonalMatrix[Lookup[rargs, {"yu", "yc", "yt"}]],
                If[
                  dirac,
                  Yv[tin] == pmns.DiagonalMatrix[Lookup[rargs, {"y1", "y2", "y3"}]],
                  Nothing
                ]
               },
               
               Right | "Right",
               {
                Yu[tin] == DiagonalMatrix[Lookup[rargs, {"yu", "yc", "yt"}]].ckm,
                If[
                  dirac,
                  Yv[tin] == DiagonalMatrix[Lookup[rargs, {"y1", "y2", "y3"}]].ConjugateTranspose[pmns],
                  Nothing
                ]
               }
             ],
             
             If[
               majorana,
               kappa[tin] == Conjugate[pmns].DiagonalMatrix[Lookup[rargs, {"kappa1", "kappa2", "kappa3"}]].ConjugateTranspose[pmns],
               Nothing
             ]
            }
          ],
          
          Yd[tin] == DiagonalMatrix[Lookup[rargs, {"yd", "ys", "yb"}]],
          Ye[tin] == DiagonalMatrix[Lookup[rargs, {"ye", "ymu", "ytau"}]]
         }
         
         (*
         If[
           majorana,
           WhenEvent[
             Mod[t, 1], 
             kappa[t] \[Rule] Symmetrize[kappa[t], Symmetric @ {1, 2}],
             "DetectionMethod" \[Rule] "Sign",
             "LocationMethod" \[Rule] "StepEnd",
             "IntegrateEvent" \[Rule] False
           ],
           Nothing
         ]
         *)
        
        } // Flatten,
        
        {
         g1,
         g2,
         g3,
         Yu,
         Yd,
         Ye,
         If[dirac, Yv, Nothing],
         If[majorana, kappa, Nothing]
        },
        
        {t, tmin, tmax},
        
        FilterRules[{opts}, Options @ NDSolve],
        
        If[
          majorana, 
          StepMonitor :> (kappa[t] = Symmetrize[kappa[t], Symmetric @ {1, 2}]), 
          Sequence
        ]
      ] // First,
      {All, 1}
    ] // Association
        
  ]
]


(* ::Subsection:: *)
(*iRenormalizeModel*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ iRenormalizeModel = Options @ nRenormalizeModel;


(* ::Subsubsection:: *)
(*SM*)


(* ::Subsubsubsection::Closed:: *)
(*NDSolve*)


mem : iRenormalizeModel[{"SM", "NDSolve"}, args_Association, opts : OptionsPattern[]] := mem = With[
  {
   dirac = ContainsAll[Keys @ args, {"y1", "y2", "y3"}],
   majorana = ContainsAll[Keys @ args, {"kappa1", "kappa2", "kappa3"}],
   qflavor = KeyExistsQ[args, "CKM"],
   lflavor = KeyExistsQ[args, "PMNS"],
   
   ordering = Which[
     ContainsAll[Keys @ args, {"y1", "y2", "y3"}],
     Lookup[args, {"y1", "y2", "y3"}],
     
     ContainsAll[Keys @ args, {"kappa1", "kappa2", "kappa3"}],
     Lookup[args, {"kappa1", "kappa2", "kappa3"}]
   ] // Ordering
  },
  
  Block[
    {$Context = $context},
  
    With[
      {sol = nRenormalizeModel[{"SM", "NDSolve"}, args, opts]},
      
      Association[
        KeyTake[sol, {"g1", "g2", "g3"}],
        
        Thread[{"yu", "yc", "yt"} -> yukawaSingularValueInterpolation[sol @ "Yu"]],
        Thread[{"yd", "ys", "yb"} -> yukawaSingularValueInterpolation[sol @ "Yd"]],
        Thread[{"ye", "ymu", "ytau"} -> yukawaSingularValueInterpolation[sol @ "Ye"]],
        
        If[
          dirac,
          Thread[{"y1", "y2", "y3"} -> yukawaSingularValueInterpolation[sol @ "Yv", ordering]],
          Nothing
        ],
        
        If[
          majorana,
          Thread[{"kappa1", "kappa2", "kappa3"} -> yukawaSingularValueInterpolation[sol @ "kappa", ordering]],
          Nothing
        ],
        
        If[
          qflavor,
          "CKM" -> mixingParametersInterpolation[
            "CKM", 
            {sol @ "Yu", sol @ "Yd"}, 
            "Parametrization" -> args["CKM", "Parametrization"]
          ],
          Nothing
        ],
        
        If[
          lflavor,
          "PMNS" -> mixingParametersInterpolation[
            "PMNS", 
            {sol @ "Ye", sol @ Which[dirac, "Yv", majorana, "kappa"]}, 
            "NeutrinoMassHierarchy" -> ordering,
            "NeutrinoMassOperator" -> Which[dirac, "Dirac", majorana, "Majorana"],
            "Parametrization" -> args["PMNS", "Parametrization"]
          ],
          Nothing
        ],
        
        KeyTake[sol, {"lambda", "m2"}]
      ]
      
    ] 
  ]
]


(* ::Subsubsection:: *)
(*MSSM*)


(* ::Subsubsubsection::Closed:: *)
(*NDSolve*)


mem : iRenormalizeModel[{"MSSM", "NDSolve"}, args_Association, opts : OptionsPattern[]] := mem = With[
  {
   dirac = ContainsAll[Keys @ args, {"y1", "y2", "y3"}],
   majorana = ContainsAll[Keys @ args, {"kappa1", "kappa2", "kappa3"}],
   qflavor = KeyExistsQ[args, "CKM"],
   lflavor = KeyExistsQ[args, "PMNS"],
   
   ordering = Which[
     ContainsAll[Keys @ args, {"y1", "y2", "y3"}],
     Lookup[args, {"y1", "y2", "y3"}],
     
     ContainsAll[Keys @ args, {"kappa1", "kappa2", "kappa3"}],
     Lookup[args, {"kappa1", "kappa2", "kappa3"}]
   ] // Ordering
  },
  
  Block[
    {$Context = $context},
  
    With[
      {sol = nRenormalizeModel[{"MSSM", "NDSolve"}, args, opts]},
      
      Association[
        KeyTake[sol, {"g1", "g2", "g3"}],
        
        Thread[{"yu", "yc", "yt"} -> yukawaSingularValueInterpolation[sol @ "Yu"]],
        Thread[{"yd", "ys", "yb"} -> yukawaSingularValueInterpolation[sol @ "Yd"]],
        Thread[{"ye", "ymu", "ytau"} -> yukawaSingularValueInterpolation[sol @ "Ye"]],
        
        If[
          dirac,
          Thread[{"y1", "y2", "y3"} -> yukawaSingularValueInterpolation[sol @ "Yv", ordering]],
          Nothing
        ],
        
        If[
          majorana,
          Thread[{"kappa1", "kappa2", "kappa3"} -> yukawaSingularValueInterpolation[sol @ "kappa", ordering]],
          Nothing
        ],
        
        If[
          qflavor,
          "CKM" -> mixingParametersInterpolation[
            "CKM", 
            {sol @ "Yu", sol @ "Yd"}, 
            "Parametrization" -> args["CKM", "Parametrization"]
          ],
          Nothing
        ],
        
        If[
          lflavor,
          "PMNS" -> mixingParametersInterpolation[
            "PMNS",
            {sol @ "Ye", sol @ Which[dirac, "Yv", majorana, "kappa"]}, 
            "NeutrinoMassHierarchy" -> ordering,
            "NeutrinoMassOperator" -> Which[dirac, "Dirac", majorana, "Majorana"],
            "Parametrization" -> args["PMNS", "Parametrization"]
          ],
          Nothing
        ]
      ]
      
    ] 
  ]
]


(* ::Subsection:: *)
(*RenormalizeModel*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ RenormalizeModel = Normal[<|"ErrorTolerance" -> Automatic, "LoopOrder" -> Automatic, Options @ iRenormalizeModel|> // KeySort];


SyntaxInformation @ RenormalizeModel = {"ArgumentsPattern" -> {{_, _}, _., _, OptionsPattern[]}, "OptionNames" -> Keys @ Options @ RenormalizeModel};


(* ::Subsubsection:: *)
(*SM*)


(* ::Subsubsubsection::Closed:: *)
(*SMDR*)


RenormalizeModel[
  {"SM", "SMDR"}, 
  "TopQuarkPoleMassScale", 
  args_Association /; ContainsAll[Keys @ args, Keys @ smdrInputParameters[]], 
  opts : OptionsPattern[]
] := SMDRCalcFit[args, FilterRules[{opts}, Options @ SMDRCalcFit]]


RenormalizeModel[
  {"SM", "SMDR"}, 
  scale_? NumericQ, 
  args_Association /; ContainsAll[Keys @ args, $calcRGrunInputParameterNames], 
  opts : OptionsPattern[]
] := SMDRCalcRGrun[scale, args, FilterRules[{opts}, Options @ SMDRCalcRGrun]]


(* ::Subsubsubsection::Closed:: *)
(*NDSolve*)


RenormalizeModel[
  {"SM", "NDSolve"}, 
  args_Association /; ContainsAll[
    Keys @ args, 
    {"RenormalizationScale", "g1", "g2", "g3", "yu", "yc", "yt", "yd", "ys", "yb", "ye", "ymu", "ytau", "lambda", "m2"}
  ], 
  opts : OptionsPattern[]
] := iRenormalizeModel[{"SM", "NDSolve"}, args, FilterRules[{opts}, Options @ iRenormalizeModel]]


RenormalizeModel[
  {"SM", "NDSolve"}, 
  scale_? NumericQ, 
  args_Association /; ContainsAll[
    Keys @ args, 
    {"RenormalizationScale", "g1", "g2", "g3", "yu", "yc", "yt", "yd", "ys", "yb", "ye", "ymu", "ytau", "lambda", "m2"}
  ], 
  opts : OptionsPattern[]
] := With[
  {sol = RenormalizeModel[{"SM", "NDSolve"}, args, FilterRules[{opts}, Options @ RenormalizeModel]]},
  
  Association[
    "RenormalizationScheme" -> "MSbar", 
    "RenormalizationScale" -> scale,
    (*Normal[sol] /. if_InterpolatingFunction \[RuleDelayed] Chop[Re @ if[Log10 @ scale], $config @ "NumericTolerance"]*)
    With[
     {pos = Position[sol, _InterpolatingFunction]},
  
      ReplacePart[
        sol,
        Thread[pos -> Chop[Re @ Through[Extract[sol, pos][Log10 @ scale]], $config @ "NumericTolerance"]]
      ]
    
    ]
  ]
  
]


(* ::Subsubsection:: *)
(*MSSM*)


(* ::Subsubsubsection::Closed:: *)
(*NDSolve*)


RenormalizeModel[
  {"MSSM", "NDSolve"}, 
  args_Association /; ContainsAll[
    Keys @ args, 
    {"RenormalizationScale", "g1", "g2", "g3", "yu", "yc", "yt", "yd", "ys", "yb", "ye", "ymu", "ytau"}
  ], 
  opts : OptionsPattern[]
] := iRenormalizeModel[{"MSSM", "NDSolve"}, args, FilterRules[{opts}, Options @ iRenormalizeModel]]


RenormalizeModel[
  {"MSSM", "NDSolve"}, 
  scale_? NumericQ, 
  args_Association /; ContainsAll[
    Keys @ args, 
    {"RenormalizationScale", "g1", "g2", "g3", "yu", "yc", "yt", "yd", "ys", "yb", "ye", "ymu", "ytau"}
  ], 
  opts : OptionsPattern[]
] := With[
  {sol = RenormalizeModel[{"MSSM", "NDSolve"}, args, FilterRules[{opts}, Options @ RenormalizeModel]]},
  
  Association[
    "RenormalizationScheme" -> "DRbar", 
    "RenormalizationScale" -> scale,
    KeyTake[args, {"SupersymmetryScale", "TanBeta"}],
    (*Normal[sol] /. if_InterpolatingFunction \[RuleDelayed] Chop[Re @ if[Log10 @ scale], $config @ "NumericTolerance"]*)
    With[
     {pos = Position[sol, _InterpolatingFunction]},
  
      ReplacePart[
        sol,
        Thread[pos -> Chop[Re @ Through[Extract[sol, pos][Log10 @ scale]], $config @ "NumericTolerance"]]
      ]
    
    ]
  ]
  
]


(* ::Section:: *)
(*RunningCouplings*)


(* ::Subsection::Closed:: *)
(*uRunningCouplings*)


Options @ uRunningCouplings = Options @ RenormalizeModel;


uRunningCouplings[
  model_,
  scale_? NumericQ,
  inputparams_Association,
  opts : OptionsPattern[]
] := With[
  {
   iparams = flattenAssociation[inputparams /. a_Around :> a @ "Value"],
   keys = Keys @ Select[inputparams // flattenAssociation, MatchQ[#, _Around] &],
   quadsum = {Total[Select[#, Negative]^2], Total[Select[#, Positive]^2]} &
  },
  
  With[
    {
     centralvalues = RenormalizeModel[
       model, 
       scale, 
       iparams // nestAssociation, 
       FilterRules[{opts}, Options @ RenormalizeModel]
     ] // flattenAssociation,
     
     uncertainties = KeyValueMap[
       {
        RenormalizeModel[
          model, 
          scale, 
          nestAssociation[<|iparams, #1 -> Min[#2 @ "Interval"]|>], 
          FilterRules[{opts}, Options @ RenormalizeModel]
        ] // flattenAssociation,
        RenormalizeModel[
          model, 
          scale, 
          nestAssociation[<|iparams, #1 -> Max[#2 @ "Interval"]|>], 
          FilterRules[{opts}, Options @ RenormalizeModel]
        ] // flattenAssociation 
       } &,
       Select[inputparams // flattenAssociation, MatchQ[#, _Around] &]
     ] // Flatten     
    },
    
    Association[
      centralvalues,
      
      Merge[
        {
         KeyValueMap[(#1 -> toPreciseNumber[#2, iparams[#1]]) &, KeyTake[centralvalues, keys]],
         Merge[KeyTake[uncertainties, keys], Identity]
        },
        
        Around[
          #[[1]], 
          Sqrt[quadsum[toRoundedNumber[#[[2]]] - toRoundedNumber[#[[1]]]]] /. {x_, y_} /; x == y :> x
        ] &
      ]
    ] // nestAssociation
    
  ]
]


(* ::Input:: *)
(*uRunningCouplings[*)
(*  model_,*)
(*  scale_? NumericQ,*)
(*  inputparams_Association,*)
(*  opts : OptionsPattern[]*)
(*] := With[*)
(*  {*)
(*   iparams = flattenAssociation[inputparams /. a_Around :> a @ "Value"],*)
(*   keys = Keys @ Select[inputparams // flattenAssociation, MatchQ[#, _Around] &],*)
(*   qmixing = KeyExistsQ[inputparams, "CKM"],*)
(*   lmixing = KeyExistsQ[inputparams, "PMNS"],*)
(*   quadsum = {Total[Select[#, Negative]^2], Total[Select[#, Positive]^2]} &,*)
(*   quadmax = {Max[0, Select[#, Negative]^2], Max[0, Select[#, Positive]^2]} &*)
(*  },*)
(*  *)
(*  With[*)
(*    {*)
(*     centralvalues = RenormalizeModel[model, scale, iparams, FilterRules[{opts}, Options @ RenormalizeModel]] // flattenAssociation,*)
(*     *)
(*     couplings = KeyValueMap[*)
(*       {*)
(*        RenormalizeModel[*)
(*          model, *)
(*          scale, *)
(*          nestAssociation[<|iparams, #1 -> Min[#2 @ "Interval"]|>], *)
(*          FilterRules[{opts}, Options @ RenormalizeModel]*)
(*        ] // flattenAssociation,*)
(*        RenormalizeModel[*)
(*          model, *)
(*          scale, *)
(*          nestAssociation[<|iparams, #1 -> Max[#2 @ "Interval"]|>], *)
(*          FilterRules[{opts}, Options @ RenormalizeModel]*)
(*        ] // flattenAssociation *)
(*       } &,*)
(*       Select[inputparams, MatchQ[#, _Around] &]*)
(*     ] // Flatten,*)
(*     *)
(*     ckm = If[*)
(*       qmixing,*)
(*       *)
(*       With[*)
(*         {param = Select[inputparams @ "CKM", MatchQ[#, _Around] &]},*)
(*         *)
(*         KeyValueMap[*)
(*           KeyDrop[*)
(*             {*)
(*              RenormalizeModel[*)
(*                model, *)
(*                scale, *)
(*                nestAssociation[<|iparams, {"CKM", #1} -> Min[#2 @ "Interval"]|>], *)
(*                FilterRules[{opts}, Options @ RenormalizeModel]*)
(*              ] // flattenAssociation,*)
(*              RenormalizeModel[*)
(*                model, *)
(*                scale, *)
(*                nestAssociation[<|iparams, {"CKM", #1} -> Max[#2 @ "Interval"]|>], *)
(*                FilterRules[{opts}, Options @ RenormalizeModel]*)
(*              ] // flattenAssociation*)
(*             },*)
(*             Key[{"CKM", #}] & /@ Complement[Keys @ param, {#1}]*)
(*           ] &,*)
(*           param*)
(*         ] // Flatten*)
(*       *)
(*       ],*)
(*       {}*)
(*     ],*)
(*     *)
(*     pmns = If[*)
(*       lmixing,*)
(*       *)
(*       With[*)
(*         {param = Select[inputparams @ "PMNS", MatchQ[#, _Around] &]},*)
(*         *)
(*         KeyValueMap[*)
(*           KeyDrop[*)
(*             {*)
(*              RenormalizeModel[*)
(*                model, *)
(*                scale, *)
(*                nestAssociation[<|iparams, {"PMNS", #1} -> Min[#2 @ "Interval"]|>], *)
(*                FilterRules[{opts}, Options @ RenormalizeModel]*)
(*              ] // flattenAssociation,*)
(*              RenormalizeModel[*)
(*                model, *)
(*                scale, *)
(*                nestAssociation[<|iparams, {"PMNS", #1} -> Max[#2 @ "Interval"]|>], *)
(*                FilterRules[{opts}, Options @ RenormalizeModel]*)
(*              ] // flattenAssociation*)
(*             },*)
(*             Key[{"PMNS", #}] & /@ Complement[Keys @ param, {#1}]*)
(*           ] &,*)
(*           param*)
(*         ] // Flatten*)
(*       *)
(*       ],*)
(*       *)
(*       {}*)
(*     ]*)
(*    },*)
(*    *)
(*    Association[*)
(*      centralvalues,*)
(**)
(*      Merge[*)
(*        {*)
(*         KeyValueMap[(#1 -> toPreciseNumber[#2, iparams[#1]]) &, KeyTake[centralvalues, keys]],*)
(*         Merge[KeyTake[couplings, keys], Identity],*)
(*         Merge[KeyTake[ckm, keys], Identity],*)
(*         Merge[KeyTake[pmns, keys], Identity]*)
(*        },*)
(*        *)
(*        Around[*)
(*          #[[1]], *)
(*          Sqrt[*)
(*            {*)
(*             quadsum[toRoundedNumber[#[[2]]] - toRoundedNumber[#[[1]]]],*)
(*             If[qmixing, quadmax[toRoundedNumber[#[[3]]] - toRoundedNumber[#[[1]]]], Nothing],*)
(*             If[lmixing, quadmax[toRoundedNumber[#[[4]]] - toRoundedNumber[#[[1]]]], Nothing]*)
(*            } // Total*)
(*          ] /. {x_, y_} /; x == y :> x*)
(*        ] &*)
(*      ]*)
(*    ] // nestAssociation*)
(*    *)
(*  ]*)
(*]*)


uRunningCouplings[
  model_,
  scale_? NumericQ,
  {inputlow_Association, inputparams_Association, inputhigh_Association},
  opts : OptionsPattern[]
] := With[
  {
   iparams = inputparams /. a_Around :> a @ "Value",
   keys = Keys @ Select[inputparams // flattenAssociation, MatchQ[#, _Around] &],
   quadmax = {Max[0, Select[#, Negative]^2], Max[0, Select[#, Positive]^2]} &
  },
  
  With[
    {
     params = uRunningCouplings[model, scale, inputparams, opts] // flattenAssociation,
     
     thresholds = Flatten[
       {
        RenormalizeModel[
          model, 
          scale, 
          inputlow /. a_Around :> a @ "Value", 
          FilterRules[{opts}, Options @ RenormalizeModel]
        ] // flattenAssociation,
        RenormalizeModel[
          model, 
          scale, 
          inputhigh /. a_Around :> a @ "Value", 
          FilterRules[{opts}, Options @ RenormalizeModel]
        ] // flattenAssociation 
       }
     ]
    },
    
    Association[
      params,
      Merge[
        {
         KeyTake[params, keys],
         Merge[KeyTake[thresholds, keys], Identity]
        },
        Around[
          #[[1]]["Value"], 
          #[[1]]["Uncertainty"] + Sqrt[quadmax[toRoundedNumber[#[[2]]] - toRoundedNumber[#[[1]]]]] /. {x_, y_} /; x == y :> x
        ] &
        
      ]
    ] // nestAssociation
    
  ]
]


(* ::Input:: *)
(*uRunningCouplings[*)
(*  model_,*)
(*  scale_,*)
(*  inputparams_,*)
(*  opts : OptionsPattern[]*)
(*] := With[*)
(*  {*)
(*   iparams = flattenAssociation @ inputparams,*)
(*   matching = OptionValue @ "MatchingScaleUncertainty" /. Automatic -> False*)
(*  },*)
(*  *)
(*  With[*)
(*    {*)
(*     params = NAroundApply[*)
(*       flattenAssociation[*)
(*         RenormalizeModel[*)
(*           model, *)
(*           scale, *)
(*           nestAssociation @ AssociationThread[Keys @ iparams, {##}],*)
(*           FilterRules[{opts}, Options @ RenormalizeModel]*)
(*         ]*)
(*       ] &,*)
(*       Values @ iparams*)
(*     ] // Association,*)
(*      *)
(*     thresholds = If[*)
(*       matching,*)
(*       NAroundApply[*)
(*         flattenAssociation[*)
(*           RenormalizeModel[*)
(*             model, *)
(*             scale, *)
(*             nestAssociation @ AssociationThread[Keys[iparams // KeySort], {##}],*)
(*             FilterRules[{opts}, Options @ RenormalizeModel]*)
(*           ]*)
(*         ] &,*)
(*         Values[*)
(*           Association[*)
(*             iparams /. a_Around :> a @ "Value",*)
(*              *)
(*             "RenormalizationScale" -> Around[*)
(*               inputparams @ "RenormalizationScale",*)
(*               {inputparams["RenormalizationScale"] / 2, inputparams["RenormalizationScale"]}*)
(*             ]*)
(*             *)
(*           ] // KeySort*)
(*         ]*)
(*       ] // Association,*)
(*       <||>*)
(*     ]*)
(*    },*)
(*    *)
(*    Association[*)
(*      params,*)
(*      Merge[*)
(*        KeyTake[{params, thresholds}, Keys @ Select[iparams, MatchQ[#, _Around] &]],*)
(*        Around[#[[1]]["Value"], Plus @@ Through[Cases[#, _Around] @ "Uncertainty"]] &*)
(*      ]*)
(*    ] // nestAssociation*)
(*    *)
(*  ]*)
(*]*)


(* ::Subsection:: *)
(*iRunningCouplings*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ iRunningCouplings = Normal[
  Association[
    "ThresholdUncertaintyFactor" -> Automatic, 
    "SupersymmetryScale" -> Automatic, 
    Options @ matchModels, 
    Options @ uRunningCouplings
  ] // KeySort
];


(* ::Subsubsection:: *)
(*SM*)


(* ::Subsubsubsection::Closed:: *)
(*SMDR*)


mem : iRunningCouplings[
  {"SM", "SMDR"}, 
  "TopQuarkPoleMassScale",
  opts : OptionsPattern[]
] := mem = Monitor[

  With[
    {inputfile = FileNameJoin @ {RunningCouplings`Developer`$PackageDirectory, "include", "SM", "smdr-runningcouplings-mt.wl"}},
  
    If[
    
      FileExistsQ @ inputfile, 
      
      Import[inputfile, "Package"],
      
      With[
        {
         output = uRunningCouplings[
           {"SM", "SMDR"}, 
           "TopQuarkPoleMassScale", 
           smdrInputParameters[],
           FilterRules[{opts}, Options @ uRunningCouplings]  
         ]
        },
      
        Export[
          FileNameJoin @ {RunningCouplings`Developer`$PackageDirectory, "include", "SM", "smdr-runningcouplings-mt.wl"},
          output
        ];
      
        output
      
      ]
    
    ]
  
  ],
    
  Row @ {"Calculating MSbar parameters at top-quark mass scale ", ProgressIndicator[Appearance -> "Ellipsis"]}
    
]


mem : iRunningCouplings[
  {"SM", "SMDR"}, 
  scale_? NumericQ, 
  opts : OptionsPattern[]
] := mem = uRunningCouplings[
  {"SM", "SMDR"}, 
  scale, 
  iRunningCouplings[{"SM", "SMDR"}, "TopQuarkPoleMassScale", opts], 
  FilterRules[{opts}, Options @ uRunningCouplings]
]


(* ::Subsubsubsection::Closed:: *)
(*NDSolve*)


mem : iRunningCouplings[
  {"SM", "NDSolve"},
  scale_? NumericQ, 
  opts : OptionsPattern[]
] := mem = uRunningCouplings[
  {"SM", "NDSolve"}, 
  scale, 
  matchModels[
    {"SM", "SMDR"} -> {"SM", "NDSolve"}, 
    iRunningCouplings[{"SM", "SMDR"}, 91.19, opts], 
    FilterRules[{opts}, Options @ matchModels]
  ], 
  FilterRules[{opts}, Options @ uRunningCouplings]
]


(* ::Subsubsection:: *)
(*MSSM*)


(* ::Subsubsubsection::Closed:: *)
(*NDSolve*)


mem : iRunningCouplings[
  {"MSSM", "NDSolve"},
  scale_? NumericQ, 
  opts : OptionsPattern[]
] := mem = With[
  {
   msusy = OptionValue @ "SupersymmetryScale" /. Automatic -> Min[scale, 10^4],
   factor = OptionValue @ "ThresholdUncertaintyFactor" /. Automatic -> 2
  },
  
  uRunningCouplings[
    {"MSSM", "NDSolve"}, 
    scale, 
    {
     matchModels[
       {"SM", "NDSolve"} -> {"MSSM", "NDSolve"}, 
       iRunningCouplings[{"SM", "NDSolve"}, msusy / factor, opts], 
       FilterRules[{opts}, Options @ matchModels]
     ],
     matchModels[
       {"SM", "NDSolve"} -> {"MSSM", "NDSolve"}, 
       iRunningCouplings[{"SM", "NDSolve"}, msusy, opts], 
       FilterRules[{opts}, Options @ matchModels]
     ],
     matchModels[
       {"SM", "NDSolve"} -> {"MSSM", "NDSolve"}, 
       iRunningCouplings[{"SM", "NDSolve"}, factor msusy, opts], 
       FilterRules[{opts}, Options @ matchModels]
     ]
    },
    FilterRules[{opts}, Options @ uRunningCouplings]
  ]

]


(* ::Subsection:: *)
(*RunningCouplings*)


(* ::Subsubsubsection::Closed:: *)
(*info*)


Options @ RunningCouplings = Normal[<|"RenormalizationMethod" -> Automatic, Options @ iRunningCouplings|> // KeySort];


SyntaxInformation @ RunningCouplings = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> Keys @ Options @ RunningCouplings};


(* ::Subsubsubsection::Closed:: *)
(*Around behaviors*)


RunningCouplings[{model_, "Values"}, scale_, opts : OptionsPattern[]] := With[
  {couplings = RunningCouplings[model, scale, opts]},

  With[
    {pos = Position[couplings, _Around]},
  
    ReplacePart[
      couplings,
      Thread[pos -> (Extract[couplings, pos] /. a_Around :> a @ "Value")]
    ]
    
  ]
]


RunningCouplings[{model_, "Uncertainties"}, scale_, opts : OptionsPattern[]] := With[
  {couplings = RunningCouplings[model, scale, opts]},

  With[
    {pos = Position[couplings, _Around]},
  
    ReplacePart[
      couplings,
      Thread[pos -> (Extract[couplings, pos] /. a_Around :> a @ "Uncertainty")]
    ]
    
  ]
]


RunningCouplings[{model_, "Invervals"}, scale_, opts : OptionsPattern[]] := With[
  {couplings = RunningCouplings[model, scale, opts]},

  With[
    {pos = Position[couplings, _Around]},
  
    ReplacePart[
      couplings,
      Thread[pos -> (Extract[couplings, pos] /. a_Around :> a @ "Interval")]
    ]
    
  ]
]


RunningCouplings[{model_, "Ranges"}, scale_, opts : OptionsPattern[]] := With[
  {couplings = RunningCouplings[model, scale, opts]},

  With[
    {pos = Position[couplings, _Around]},
  
    ReplacePart[
      couplings,
      Thread[pos -> (Extract[couplings, pos] /. a_Around :> MinMax[a @ "Inverval"])]
    ]
    
  ]
]


(* ::Subsubsection::Closed:: *)
(*SMDR*)


RunningCouplings["SM", "TopQuarkPoleMassScale", opts : OptionsPattern[]] := iRunningCouplings[{"SM", "SMDR"}, 
  "TopQuarkPoleMassScale", 
  FilterRules[{opts}, Options @ iRunningCouplings]
]


(* ::Subsubsection::Closed:: *)
(*NDSolve*)


RunningCouplings[
  model : ("SM" | "MSSM"), 
  scale_? NumericQ, 
  opts : OptionsPattern[]
] := iRunningCouplings[
  {model, OptionValue @ "RenormalizationMethod" /. Automatic -> "NDSolve"}, 
  scale, 
  FilterRules[{opts}, Options @ iRunningCouplings]
]


(* ::Section:: *)
(*running quantities*)


(* ::Subsection::Closed:: *)
(*RunningGaugeCoupling*)


Options @ RunningGaugeCoupling = Options @ RunningCouplings;


SyntaxInformation @ RunningGaugeCoupling = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> Keys[Options @ RunningGaugeCoupling]};


RunningGaugeCoupling[{model_, 1}, scale_? NumericQ, opts : OptionsPattern[]] := RunningCouplings[model, scale, opts]["g1"]
RunningGaugeCoupling[{model_, 2}, scale_? NumericQ, opts : OptionsPattern[]] := RunningCouplings[model, scale, opts]["g2"]
RunningGaugeCoupling[{model_, 3}, scale_? NumericQ, opts : OptionsPattern[]] := RunningCouplings[model, scale, opts]["g3"]


(* ::Subsection::Closed:: *)
(*RunningYukawaCoupling*)


Options @ RunningYukawaCoupling = Options @ RunningCouplings;


SyntaxInformation @ RunningYukawaCoupling = {"ArgumentsPattern"->{_, _, OptionsPattern[]}, "OptionNames" -> Keys[Options @ RunningYukawaCoupling]};


RunningYukawaCoupling[{model_, "E"}, scale_? NumericQ, opts : OptionsPattern[]] := Lookup[RunningCouplings[model, scale, opts], {"ye", "ymu", "ytau"}]
RunningYukawaCoupling[{model_, "U"}, scale_? NumericQ, opts : OptionsPattern[]] := Lookup[RunningCouplings[model, scale, opts], {"yu", "yc", "yt"}]
RunningYukawaCoupling[{model_, "D"}, scale_? NumericQ, opts : OptionsPattern[]] := Lookup[RunningCouplings[model, scale, opts], {"yd", "ys", "yb"}]
RunningYukawaCoupling[{model_, "V"}, scale_? NumericQ, opts : OptionsPattern[]] := Lookup[RunningCouplings[model, scale, opts, "NeutrinoMassOperator" -> "Dirac"], {"y1", "y2", "y3"}]


RunningYukawaCoupling[
  {model_, spec : ("E" | "U" | "D" | "V"), i : (1 | 2 | 3)}, 
  scale_? NumericQ, 
  opts : OptionsPattern[]
] := RunningYukawaCoupling[{model, spec}, scale, opts][[i]]


RunningYukawaCoupling[{model_, "e"}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "E", 1}, scale, opts]
RunningYukawaCoupling[{model_, ("mu" | "\[Mu]")}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "E", 2}, scale, opts]
RunningYukawaCoupling[{model_, ("tau" | "\[Tau]")}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "E", 3}, scale, opts]


RunningYukawaCoupling[{model_, "u"}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "U", 1}, scale, opts]
RunningYukawaCoupling[{model_, "c"}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "U", 2}, scale, opts]
RunningYukawaCoupling[{model_, "t"}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "U", 3}, scale, opts]


RunningYukawaCoupling[{model_, "d"}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "D", 1}, scale, opts]
RunningYukawaCoupling[{model_, "s"}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "D", 2}, scale, opts]
RunningYukawaCoupling[{model_, "b"}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "D", 3}, scale, opts]


RunningYukawaCoupling[{model_, ("v1" | "\[ScriptV]1")}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "V", 1}, scale, opts]
RunningYukawaCoupling[{model_, ("v2" | "\[ScriptV]2")}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "V", 2}, scale, opts]
RunningYukawaCoupling[{model_, ("v3" | "\[ScriptV]3")}, scale_? NumericQ, opts : OptionsPattern[]] := RunningYukawaCoupling[{model, "V", 3}, scale, opts]


(* ::Subsection::Closed:: *)
(*RunningWienbergOperator*)


Options @ RunningWeinbergOperator = Options @ RunningCouplings;


SyntaxInformation @ RunningWeinbergOperator = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> Keys[Options @ RunningWeinbergOperator]};


RunningWeinbergOperator[model_, scale_? NumericQ, opts : OptionsPattern[]] := Lookup[RunningCouplings[model, scale, opts], {"kappa1", "kappa2", "kappa3"}]


RunningWeinbergOperator[{model_, i : (1 | 2 | 3)}, scale_? NumericQ, opts : OptionsPattern[]] := RunningWeinbergOperator[model, scale, opts, "NeutrinoMassOperator" -> "Majorana"][[i]]


(* ::Subsection::Closed:: *)
(*RunningMixingParameter*)


Options @ RunningMixingParameter = Options @ RunningCouplings;


SyntaxInformation @ RunningMixingParameter = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}, "OptionNames" -> Keys[Options @ RunningMixingParameter]};


RunningMixingParameter[{model_, spec : ("CKM" | "PMNS")}, scale_? NumericQ, opts : OptionsPattern[]] := RunningCouplings[model, scale, opts][spec]


RunningMixingParameter[{model_, spec : ("CKM" | "PMNS"), param_String}, scale_? NumericQ, opts : OptionsPattern[]] := RunningCouplings[model, scale, opts][spec, param]


(* ::Chapter:: *)
(*end package*)


(* ::Subsubsection::GrayLevel[0]::Closed:: *)
(*end private context*)


End[];


(* ::Subsubsection::GrayLevel[0]::Closed:: *)
(*end package context*)


EndPackage[];
