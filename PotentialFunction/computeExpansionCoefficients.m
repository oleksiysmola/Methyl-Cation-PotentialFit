(* D3h(M) Symmetry Operations *)
(* Note that technically there are 6 symmetry operations but by choosing sin rho as our dihedral 
function, this degree of freedom now transforms as A1' and so trivially satisfies the required 
properties *)
numberOfSymmetryOperations = 3;
Subscript[S, E] := IdentityMatrix[6];
Subscript[S, 23] := {
{1, 0, 0, 0, 0, 0},
{0, 0, 1, 0, 0, 0},
{0, 1, 0, 0, 0, 0},
{0, 0, 0, 1, 0, 0},
{0, 0, 0, 0, -1, 0},
{0, 0, 0, 0, 0, 1}
};
Subscript[S, 123] := {
{0, 1, 0, 0, 0, 0},
{0, 0, 1, 0, 0, 0},
{1, 0, 0, 0, 0, 0},
{0, 0, 0, -1/2, Sqrt[3]/2, 0},
{0, 0, 0, -Sqrt[3]/2, -1/2, 0},
{0, 0, 0, 0, 0, 1}
};
MatrixForm[Subscript[S, 123]]; 
MatrixForm[Subscript[S, 23]]; 
symmetryOperations = {Subscript[S, E], Subscript[S, 123], Subscript[S, 23]};

(* Setup Coefficients *)
maxOrder = 6;
maxMultiMode = 2;
numberOfRigidModes = 6;
(* 1st Order *)
powers = ConstantArray[0, numberOfRigidModes];
firstOrderCoefficients = ConstantArray[0, numberOfRigidModes];
firstOrderVariables = {};
For[j = 1, j <= numberOfRigidModes, j++,
  powers[[j]] = powers[[j]] + 1;
  powersString = "";
  For[k = 1, k <= numberOfRigidModes, k++,
      powersString = powersString <> ToString[powers[[k]]];
   ];
   firstOrderCoefficients[[j]] = 0 + Subscript[f, powersString];
   firstOrderVariables = Append[firstOrderVariables, 2];
   firstOrderVariables = Append[firstOrderVariables, Subscript[f, powersString]];
   powers = ConstantArray[0, {numberOfRigidModes}];
];
MatrixForm[firstOrderCoefficients];
Dimensions[firstOrderCoefficients];
MatrixForm[firstOrderVariables];
MatrixForm[firstOrderCoefficients[[1]]];
Dimensions[firstOrderCoefficients[[1]]];
firstOrderVariables = DeleteDuplicates[firstOrderVariables];
firstOrderVariables = Delete[firstOrderVariables, 1];
(* Print["First order coefficients"]
Print[firstOrderCoefficients]
Print[firstOrderVariables]
Print[firstOrderVariables[[2]][[2]]] *)
(* 2nd Order *)
powers = ConstantArray[0, {numberOfRigidModes}];
secondOrderCoefficients = ConstantArray[0, {numberOfRigidModes, numberOfRigidModes}];
secondOrderVariables = {};

For[j = 1, j <= numberOfRigidModes, j++,
  For[l = 1, l <= numberOfRigidModes, l++,
   powers[[j]] = powers[[j]] + 1;
   powers[[l]] = powers[[l]] + 1;
   powersString = "";
   For[k = 1, k <= numberOfRigidModes, k++,
       powersString = powersString <> ToString[powers[[k]]];
    ];
    secondOrderCoefficients[[j, l]] = 0 + Subscript[f, powersString];
    secondOrderVariables = Append[secondOrderVariables, 2];
    secondOrderVariables = Append[secondOrderVariables, 0 + Subscript[f, powersString]];
    powers = ConstantArray[0, {numberOfRigidModes}];
  ];
];
secondOrderVariables = DeleteDuplicates[secondOrderVariables];
MatrixForm[secondOrderCoefficients];
MatrixForm[secondOrderVariables];
secondOrderVariables = Delete[secondOrderVariables, 1];
Print["Second order coefficients"]
Print[secondOrderCoefficients]
Print[secondOrderVariables]
Print[secondOrderVariables[[2]][[2]]]
(* 3rd Order *)
powers = ConstantArray[0, {numberOfRigidModes}];
thirdOrderCoefficients = ConstantArray[0, {numberOfRigidModes, numberOfRigidModes, numberOfRigidModes}];
thirdOrderVariables = {};
multiModeVector = ConstantArray[0, {numberOfRigidModes}];
multiMode = 0;

For[j = 1, j <= numberOfRigidModes, j++,
  For[l = 1, l <= numberOfRigidModes, l++,
     For[m = 1, m <= numberOfRigidModes, m++,
      powers[[j]] = powers[[j]] + 1;
      powers[[l]] = powers[[l]] + 1;
      powers[[m]] = powers[[m]] + 1;
      powersString = "";
      For[k = 1, k <= numberOfRigidModes, k++,
         If[k =/= 4, If[k =/= 5, If[powers[[k]] > 0, multiMode = multiMode + 1,],],];
          powersString = powersString <> ToString[powers[[k]]];
       ];
       If[multiMode < maxMultiMode,
          fourthOrderCoefficients[[j, l, m, n]] = 0 + Subscript[f, powersString];
          fourthOrderVariables = Append[fourthOrderVariables, 2];
          fourthOrderVariables = Append[fourthOrderVariables, 0 + Subscript[f, powersString]];
          powers = ConstantArray[0, {numberOfRigidModes}];,
        ];
     ]:
  ];
];
thirdOrderVariables = DeleteDuplicates[thirdOrderVariables];
MatrixForm[thirdOrderCoefficients];
MatrixForm[thirdOrderVariables];
Print["Third order coefficients"]
thirdOrderVariables = Delete[thirdOrderVariables, 1];
Print[thirdOrderCoefficients]
Print[thirdOrderVariables]
Print[thirdOrderVariables[[2]][[2]]]
(* 4th Order *)
powers = ConstantArray[0, {numberOfRigidModes}];
fourthOrderCoefficients = ConstantArray[0, {maxFourierMode + 1, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes}];
fourthOrderVariables = {};
multiModeVector = ConstantArray[0, {numberOfRigidModes}];
multiMode = 0;

For[j = 1, j <= numberOfRigidModes, j++,
  For[l = 1, l <= numberOfRigidModes, l++,
     For[m = 1, m <= numberOfRigidModes, m++,
       powers[[j]] = powers[[j]] + 1;
       powers[[l]] = powers[[l]] + 1;
       powers[[m]] = powers[[m]] + 1;
       powersString = "";
       For[k = 1, k <= numberOfRigidModes, k++,
           powersString = powersString <> ToString[powers[[k]]];
        ];
        thirdOrderCoefficients[[j, l, m]] = 0 + Subscript[f, powersString];
        thirdOrderVariables = Append[thirdOrderVariables, 2];
        thirdOrderVariables = Append[thirdOrderVariables, 0 + Subscript[f, powersString]];
        powers = ConstantArray[0, {numberOfRigidModes}];
      ];
  ];
];
fourthOrderVariables = DeleteDuplicates[fourthOrderVariables];
MatrixForm[fourthOrderCoefficients];
MatrixForm[fourthOrderVariables];
(* 5th Order *)
powers = ConstantArray[0, {numberOfRigidModes}];
fifthOrderCoefficients = ConstantArray[0, {maxFourierMode + 1, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes}];
fifthOrderVariables = {};
multiModeVector = ConstantArray[0, {numberOfRigidModes}];
multiMode = 0;

For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  For[j = 1, j <= numberOfRigidModes, j++,
     For[k = 1, k <= numberOfRigidModes, k++,
        For[l = 1, l <= numberOfRigidModes, l++,
            For[n = 1, n <= numberOfRigidModes, n++,
               For[o = 1, o <= numberOfRigidModes, o++,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  multiMode = 0;
                  powers[[j]] = powers[[j]] + 1;
                  powers[[k]] = powers[[k]] + 1;
                  powers[[l]] = powers[[l]] + 1;
                  powers[[n]] = powers[[n]] + 1;
                  powers[[o]] = powers[[o]] + 1;
                  powersString = "";
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                        multiModeVector[[m]] = 0];
                     powersString = powersString <> ToString[powers[[m]]];
	                 ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];
                  If[multiMode <= maxMultiMode, 
                     labelCos = Subscript[f, powersString, \[Alpha]];
	                  labelSin = Subscript[h, powersString, \[Alpha]];
	                  If[Count[powers, 5] == 1, 
                        fifthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o]] = {Subscript[f, powersString, \[Alpha]], Subscript[h, powersString, \[Alpha]]},
                        If[Count[powers, 4] == 1,
                           fifthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o]] = {Subscript[f, powersString, \[Alpha]]/5, Subscript[h, powersString, \[Alpha]]/5},
                           If[Count[powers, 3] == 1,
                              fifthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o]] = {Subscript[f, powersString, \[Alpha]]/10, Subscript[h, powersString, \[Alpha]]/10},
                              If[Count[powers, 1] == 5, 
                                 fifthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o]] = {Subscript[f, powersString, \[Alpha]]/(5!), Subscript[h, powersString, \[Alpha]]/(5!)};,
                                 fifthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o]] = {0, 0};
                                 ]; 
                              ];
                           ];
                        ];
                        fifthOrderVariables = Append[fifthOrderVariables, Subscript[f, powersString, \[Alpha]]];
	                     fifthOrderVariables = Append[fifthOrderVariables, Subscript[h, powersString, \[Alpha]]];, 
                        fifthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o]] = {0, 0};
                  ];
                  powers = ConstantArray[0, {numberOfRigidModes}];
               ];
            ];
        ];
     ];
   ];
];
fifthOrderVariables = DeleteDuplicates[fifthOrderVariables];
MatrixForm[fifthOrderCoefficients];
MatrixForm[fifthOrderVariables];
(* 6th Order *)
powers = ConstantArray[0, {numberOfRigidModes}];
sixthOrderCoefficients = ConstantArray[0, {maxFourierMode + 1, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes}];
sixthOrderVariables = {};
multiModeVector = ConstantArray[0, {numberOfRigidModes}];
multiMode = 0;
Print["Defining 6th order coefficients...."]
For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  For[j = 1, j <= numberOfRigidModes, j++,
     For[k = 1, k <= numberOfRigidModes, k++,
        For[l = 1, l <= numberOfRigidModes, l++,
            For[n = 1, n <= numberOfRigidModes, n++,
               For[o = 1, o <= numberOfRigidModes, o++,
                  For[p = 1, p <= numberOfRigidModes, p++,
                     multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                     multiMode = 0;
                     powers[[j]] = powers[[j]] + 1;
                     powers[[k]] = powers[[k]] + 1;
                     powers[[l]] = powers[[l]] + 1;
                     powers[[n]] = powers[[n]] + 1;
                     powers[[o]] = powers[[o]] + 1;
                     powers[[p]] = powers[[p]] + 1;
                     powersString = "";
                     For[m = 1, m <= numberOfRigidModes, m++,
                        If[powers[[m]] > 0, 
                           If[m < numberOfRigidModes - 1, 
                              multiModeVector[[m]] = 1, 
                              multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
                        powersString = powersString <> ToString[powers[[m]]];
	                    ];
                     multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];
                     If[multiMode <= maxMultiMode, 
                        labelCos = Subscript[f, powersString, \[Alpha]];
	                     labelSin = Subscript[h, powersString, \[Alpha]];
	                     If[Count[powers, 6] == 1, 
                           sixthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o, p]] = {Subscript[f, powersString, \[Alpha]], Subscript[h, powersString, \[Alpha]]},
                           If[Count[powers, 5] == 1,
                              sixthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o, p]] = {Subscript[f, powersString, \[Alpha]]/6, Subscript[h, powersString, \[Alpha]]/6},
                              If[Count[powers, 4] == 1,
                                 sixthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o, p]] = {Subscript[f, powersString, \[Alpha]]/15, Subscript[h, powersString, \[Alpha]]/15},
                                 If[Count[powers, 3] == 2, 
                                    sixthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o, p]] = {Subscript[f, powersString, \[Alpha]]/20, Subscript[h, powersString, \[Alpha]]/20};,
                                    sixthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o, p]] = {Subscript[f, powersString, \[Alpha]]/(6!), Subscript[h, powersString, \[Alpha]]/(6!)};
                                    ]; 
                                 ];
                              ];
                           ];
                           sixthOrderVariables = Append[sixthOrderVariables, Subscript[f, powersString, \[Alpha]]];
	                        sixthOrderVariables = Append[sixthOrderVariables, Subscript[h, powersString, \[Alpha]]];, 
                           sixthOrderCoefficients[[\[Alpha] + 1, j, k, l, n, o, p]] = {0, 0};
                     ];
                     powers = ConstantArray[0, {numberOfRigidModes}];
                  ];
               ];
            ];
        ];
     ];
   ];
];
sixthOrderVariables = DeleteDuplicates[sixthOrderVariables];
MatrixForm[sixthOrderCoefficients];
MatrixForm[sixthOrderVariables];
(* Solve equations *)
(* 0th Order *)
zeroOrderEquations = ConstantArray[0, {numberOfSymmetryOperations, maxFourierMode + 1}];
For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  torsionSymmetryOperations = {T[\[Alpha]], Subscript[T, 123][\[Alpha]], Subscript[T, 132][\[Alpha]], Subscript[T, 12][\[Alpha]], 
  Subscript[T, 23][\[Alpha]], Subscript[T, 13][\[Alpha]]};
  For[operationNumber = 1, operationNumber <= numberOfSymmetryOperations, operationNumber++,
      zeroOrderEquations[[operationNumber, \[Alpha] + 1]] = zeroOrderCoefficients[[\[Alpha] + 1]] - Dot[zeroOrderCoefficients[[\[Alpha] + 1]], torsionSymmetryOperations[[operationNumber]]];
      ]
]
zeroOrderEquations = Flatten[zeroOrderEquations];
zeroOrderEquations = Thread[zeroOrderEquations == 0];
zeroOrderSolutions = Solve[zeroOrderEquations, zeroOrderVariables];
zeroOrderSolutions
Print["Printing 0th order solutions..."]
For[fourierCycle = 0, fourierCycle <= 7, fourierCycle++, 
   For[i = 0, i <= Dimensions[zeroOrderVariables][[1]], i++,
      variable = zeroOrderVariables[[i]];
      If[variable == (variable /. zeroOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], " 0 0 0  ", variable[[3]] + fourierCycle*3
         ];,
         variable
      ];
   ];
];
(* 1st Order *)
firstOrderEquations = ConstantArray[0, {numberOfSymmetryOperations, maxFourierMode + 1, numberOfRigidModes}];
For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  torsionSymmetryOperations = {T[\[Alpha]], Subscript[T, 123][\[Alpha]], Subscript[T, 132][\[Alpha]], Subscript[T, 12][\[Alpha]], 
  Subscript[T, 23][\[Alpha]], Subscript[T, 13][\[Alpha]]};
  For[operationNumber = 1, operationNumber <= numberOfSymmetryOperations, operationNumber++,
      firstOrderEquations[[operationNumber, \[Alpha] + 1]] = firstOrderCoefficients[[\[Alpha] + 1]] - Transpose[TensorContract[TensorProduct[Dot[firstOrderCoefficients[[\[Alpha] + 1]], torsionSymmetryOperations[[operationNumber]]], symmetryOperations[[operationNumber]]], {1, 3}]];
      ]
]
Dimensions[firstOrderEquations];
firstOrderEquations = Flatten[firstOrderEquations];
firstOrderEquations = Thread[firstOrderEquations == 0];
firstOrderSolutions = Solve[firstOrderEquations, firstOrderVariables];
firstOrderSolutions
(* firstOrderSolutions = ParallelTable[Solve[firstOrderEquations[[i]], firstOrderVariables], {i, 1, Length[firstOrderEquations]}]; *)
(* firstOrderSolutions = ParallelMap[Solve[#, firstOrderVariables] &, firstOrderEquations]; *)
Print["Printing 1st order solutions..."]
For[fourierCycle = 0, fourierCycle <= 7, fourierCycle++,
   For[i = 0, i <= Dimensions[zeroOrderVariables][[1]], i++,
      variable = zeroOrderVariables[[i]];
      If[variable == (variable /. zeroOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   1 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
         variable
      ];
   ];
   For[i = 0, i <= Dimensions[firstOrderVariables][[1]], i++,
      variable = firstOrderVariables[[i]];
      If[variable == (variable /. firstOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3
         ];,
         variable
      ];
   ];
];
(* 2nd Order *)
secondOrderEquations = ConstantArray[0, {numberOfSymmetryOperations, maxFourierMode + 1, numberOfRigidModes, numberOfRigidModes}];
For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  torsionSymmetryOperations = {T[\[Alpha]], Subscript[T, 123][\[Alpha]], Subscript[T, 132][\[Alpha]], Subscript[T, 12][\[Alpha]], 
  Subscript[T, 23][\[Alpha]], Subscript[T, 13][\[Alpha]]};
  For[operationNumber = 1, operationNumber <= numberOfSymmetryOperations, operationNumber++,
      secondOrderEquations[[operationNumber, \[Alpha] + 1]] = secondOrderCoefficients[[\[Alpha] + 1]] - TensorContract[TensorProduct[symmetryOperations[[operationNumber]], symmetryOperations[[operationNumber]], Dot[secondOrderCoefficients[[\[Alpha] + 1]], torsionSymmetryOperations[[operationNumber]]]], {{1, 5}, {3, 6}}];
      ]
]
secondOrderEquations = Flatten[secondOrderEquations];
secondOrderEquations = Thread[secondOrderEquations == 0];
(* secondOrderSolutions = ParallelTable[Solve[secondOrderEquations[[i]], secondOrderVariables], {i, 1, Length[secondOrderEquations]}]; *)
secondOrderSolutions = Solve[secondOrderEquations, secondOrderVariables];
(* secondOrderSolutions = ParallelMap[Solve[#, secondOrderVariables] &, secondOrderEquations]; *)
secondOrderSolutions
Print["Printing 2rd order solutions..."]
For[fourierCycle = 0, fourierCycle <= 7, fourierCycle++,
   For[i = 0, i <= Dimensions[zeroOrderVariables][[1]], i++,
      variable = zeroOrderVariables[[i]];
      If[variable == (variable /. zeroOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   2 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 2 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 2 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   1 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   1 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
         variable
      ];
   ];
   For[i = 0, i <= Dimensions[firstOrderVariables][[1]], i++,
      variable = firstOrderVariables[[i]];
      If[variable == (variable /. firstOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   1 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
         variable
      ];
   ];
   For[i = 0, i <= Dimensions[secondOrderVariables][[1]], i++,
      variable = secondOrderVariables[[i]];
      If[variable == (variable /. secondOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3
         ];,
         variable
      ];
   ];
];
(* 3rd Order *)
Print["Applying symmetry operations to 3rd order coefficients..."]
thirdOrderEquations = ConstantArray[0, {numberOfSymmetryOperations, maxFourierMode + 1, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes}];
For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  torsionSymmetryOperations = {T[\[Alpha]], Subscript[T, 123][\[Alpha]], Subscript[T, 132][\[Alpha]], Subscript[T, 12][\[Alpha]], 
  Subscript[T, 23][\[Alpha]], Subscript[T, 13][\[Alpha]]};
  For[operationNumber = 1, operationNumber <= numberOfSymmetryOperations, operationNumber++,
      transformedCoefficients = Dot[thirdOrderCoefficients[[\[Alpha] + 1]], torsionSymmetryOperations[[operationNumber]]];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 3}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 4}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 5}];
      thirdOrderEquations[[operationNumber, \[Alpha] + 1]] = thirdOrderCoefficients[[\[Alpha] + 1]] - transformedCoefficients;
      ]
]
Print["Done!"]
Print["Solving 3rd order equations..."]
thirdOrderEquations = Flatten[thirdOrderEquations];
thirdOrderEquations = Thread[thirdOrderEquations == 0];
thirdOrderSolutions = Solve[thirdOrderEquations, thirdOrderVariables];
Print["Done!"]
thirdOrderSolutions
Print["Printing 3rd order solutions..."]
For[fourierCycle = 0, fourierCycle <= 7, fourierCycle++,
   For[i = 0, i <= Dimensions[zeroOrderVariables][[1]], i++,
      variable = zeroOrderVariables[[i]];
      If[variable == (variable /. zeroOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   3 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 3 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 3 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   2 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   2 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   1 2 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   1 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 2 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 2 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 2 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         If[maxMultiMode >= 3,
            Print["* ", variable[[1]], "   1 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
         ],
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[firstOrderVariables][[1]], i++,
      variable = firstOrderVariables[[i]];
      If[variable == (variable /. firstOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   2 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 2 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 2 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
         If[maxMultiMode >= 3,
            Print["* ", variable[[1]], "   1 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
            Print["* ", variable[[1]], "   1 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
            Print["* ", variable[[1]], "   0 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
            variable;
            ];,
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[secondOrderVariables][[1]], i++,
      variable = secondOrderVariables[[i]];
      If[variable == (variable /. secondOrderSolutions)[[1]],
         If[Count[Characters[variable[[2]]], "2"] == 1, 
            Print["* ", variable[[1]], "   1 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
            Print["* ", variable[[1]], "   0 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
            Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
            If[maxMultiMode >= 3, 
               Print["* ", variable[[1]], "   1 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
               Print["* ", variable[[1]], "   0 1 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];
               Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 1 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
               variable
               ];
            ];,
         variable
      ];
   ];
   For[i = 0, i <= Dimensions[thirdOrderVariables][[1]], i++,
      variable = thirdOrderVariables[[i]];
      If[variable == (variable /. thirdOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3
         ];,
         variable
      ];
   ];
];
(* 4th Order *)
Print["Applying symmetry operations to 4th order coefficients..."]
fourthOrderEquations = ConstantArray[0, {numberOfSymmetryOperations, maxFourierMode + 1, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes}];
For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  torsionSymmetryOperations = {T[\[Alpha]], Subscript[T, 123][\[Alpha]], Subscript[T, 132][\[Alpha]], Subscript[T, 12][\[Alpha]], 
  Subscript[T, 23][\[Alpha]], Subscript[T, 13][\[Alpha]]};
  For[operationNumber = 1, operationNumber <= numberOfSymmetryOperations, operationNumber++,
      transformedCoefficients = Dot[fourthOrderCoefficients[[\[Alpha] + 1]], torsionSymmetryOperations[[operationNumber]]];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 3}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 4}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 5}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 6}];
      fourthOrderEquations[[operationNumber, \[Alpha] + 1]] = fourthOrderCoefficients[[\[Alpha] + 1]] - transformedCoefficients;
      ]
]
Print["Done!"]
Print["Solving 4th order equations..."]
fourthOrderEquations = Flatten[fourthOrderEquations];
fourthOrderEquations = Thread[fourthOrderEquations == 0];
fourthOrderSolutions = Solve[fourthOrderEquations, fourthOrderVariables];
Print["Done!"]
fourthOrderSolutions
Print["Printing 4th order solutions..."]
maxPower = 4;
For[fourierCycle = 0, fourierCycle <= 7, fourierCycle++,
   For[i = 0, i <= Dimensions[zeroOrderVariables][[1]], i++,
      variable = zeroOrderVariables[[i]];
      If[variable == (variable /. zeroOrderSolutions)[[1]],
         For[j = 0, j <= maxPower, j++,
            For[k = 0, k <= maxPower, k++,
               For[l = 0, l <= maxPower, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[firstOrderVariables][[1]], i++,
      variable = firstOrderVariables[[i]];
      If[variable == (variable /. firstOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 1, j++,
            For[k = 0, k <= maxPower - 1, k++,
               For[l = 0, l <= maxPower - 1, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 1,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[secondOrderVariables][[1]], i++,
      variable = secondOrderVariables[[i]];
      If[variable == (variable /. secondOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 2, j++,
            For[k = 0, k <= maxPower - 2, k++,
               For[l = 0, l <= maxPower - 2, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 2,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[thirdOrderVariables][[1]], i++,
      variable = thirdOrderVariables[[i]];
      If[variable == (variable /. thirdOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 3, j++,
            For[k = 0, k <= maxPower - 3, k++,
               For[l = 0, l <= maxPower - 3, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 3,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[fourthOrderVariables][[1]], i++,
      variable = fourthOrderVariables[[i]];
      If[variable == (variable /. fourthOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
         variable;
      ];
   ];
];
(* 5th Order *)
Print["Applying symmetry operations to 5th order coefficients..."]
fifthOrderEquations = ConstantArray[0, {numberOfSymmetryOperations, maxFourierMode + 1, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes}];
For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  torsionSymmetryOperations = {T[\[Alpha]], Subscript[T, 123][\[Alpha]], Subscript[T, 132][\[Alpha]], Subscript[T, 12][\[Alpha]], 
  Subscript[T, 23][\[Alpha]], Subscript[T, 13][\[Alpha]]};
  For[operationNumber = 1, operationNumber <= numberOfSymmetryOperations, operationNumber++,
      transformedCoefficients = Dot[fifthOrderCoefficients[[\[Alpha] + 1]], torsionSymmetryOperations[[operationNumber]]];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 3}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 4}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 5}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 6}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 7}];
      fifthOrderEquations[[operationNumber, \[Alpha] + 1]] = fifthOrderCoefficients[[\[Alpha] + 1]] - transformedCoefficients;
      ]
]
Print["Done!"]
Print["Solving 5th order equations..."]
fifthOrderEquations = Flatten[fifthOrderEquations];
fifthOrderEquations = Thread[fifthOrderEquations == 0];
fifthOrderSolutions = Solve[fifthOrderEquations, fifthOrderVariables];
Print["Done!"]
fifthOrderSolutions
Print["Printing 5th order solutions..."]
maxPower = 5;
For[fourierCycle = 0, fourierCycle <= 7, fourierCycle++,
   For[i = 0, i <= Dimensions[zeroOrderVariables][[1]], i++,
      variable = zeroOrderVariables[[i]];
      If[variable == (variable /. zeroOrderSolutions)[[1]],
         For[j = 0, j <= maxPower, j++,
            For[k = 0, k <= maxPower, k++,
               For[l = 0, l <= maxPower, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[firstOrderVariables][[1]], i++,
      variable = firstOrderVariables[[i]];
      If[variable == (variable /. firstOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 1, j++,
            For[k = 0, k <= maxPower - 1, k++,
               For[l = 0, l <= maxPower - 1, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 1,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[secondOrderVariables][[1]], i++,
      variable = secondOrderVariables[[i]];
      If[variable == (variable /. secondOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 2, j++,
            For[k = 0, k <= maxPower - 2, k++,
               For[l = 0, l <= maxPower - 2, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 2,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[thirdOrderVariables][[1]], i++,
      variable = thirdOrderVariables[[i]];
      If[variable == (variable /. thirdOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 3, j++,
            For[k = 0, k <= maxPower - 3, k++,
               For[l = 0, l <= maxPower - 3, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 3,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[fourthOrderVariables][[1]], i++,
      variable = fourthOrderVariables[[i]];
      If[variable == (variable /. fourthOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 4, j++,
            For[k = 0, k <= maxPower - 4, k++,
               For[l = 0, l <= maxPower - 4, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 4,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[fifthOrderVariables][[1]], i++,
      variable = fifthOrderVariables[[i]];
      If[variable == (variable /. fifthOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
         variable;
      ];
   ];
];
(* 6th Order *)
Print["Applying symmetry operations to 6th order coefficients..."]
sixthOrderEquations = ConstantArray[0, {numberOfSymmetryOperations, maxFourierMode + 1, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes, numberOfRigidModes}];
For[\[Alpha] = 0, \[Alpha] <= maxFourierMode, \[Alpha]++,
  torsionSymmetryOperations = {T[\[Alpha]], Subscript[T, 123][\[Alpha]], Subscript[T, 132][\[Alpha]], Subscript[T, 12][\[Alpha]], 
  Subscript[T, 23][\[Alpha]], Subscript[T, 13][\[Alpha]]};
  For[operationNumber = 1, operationNumber <= numberOfSymmetryOperations, operationNumber++,
      transformedCoefficients = Dot[sixthOrderCoefficients[[\[Alpha] + 1]], torsionSymmetryOperations[[operationNumber]]];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 3}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 4}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 5}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 6}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 7}];
      transformedCoefficients = TensorContract[TensorProduct[symmetryOperations[[operationNumber]], transformedCoefficients], {1, 8}];
      sixthOrderEquations[[operationNumber, \[Alpha] + 1]] = sixthOrderCoefficients[[\[Alpha] + 1]] - transformedCoefficients;
      ]
]
Print["Done!"]
Print["Solving 6th order equations..."]
sixthOrderEquations = Flatten[sixthOrderEquations];
sixthOrderEquations = Thread[sixthOrderEquations == 0];
sixthOrderSolutions = Solve[sixthOrderEquations, sixthOrderVariables];
Print["Done!"]
sixthOrderSolutions
Print["Printing 6th order solutions..."]
maxPower = 6;
For[fourierCycle = 0, fourierCycle <= 7, fourierCycle++,
   For[i = 0, i <= Dimensions[zeroOrderVariables][[1]], i++,
      variable = zeroOrderVariables[[i]];
      If[variable == (variable /. zeroOrderSolutions)[[1]],
         For[j = 0, j <= maxPower, j++,
            For[k = 0, k <= maxPower, k++,
               For[l = 0, l <= maxPower, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[firstOrderVariables][[1]], i++,
      variable = firstOrderVariables[[i]];
      If[variable == (variable /. firstOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 1, j++,
            For[k = 0, k <= maxPower - 1, k++,
               For[l = 0, l <= maxPower - 1, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 1,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[secondOrderVariables][[1]], i++,
      variable = secondOrderVariables[[i]];
      If[variable == (variable /. secondOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 2, j++,
            For[k = 0, k <= maxPower - 2, k++,
               For[l = 0, l <= maxPower - 2, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 2,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[thirdOrderVariables][[1]], i++,
      variable = thirdOrderVariables[[i]];
      If[variable == (variable /. thirdOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 3, j++,
            For[k = 0, k <= maxPower - 3, k++,
               For[l = 0, l <= maxPower - 3, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 3,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[fourthOrderVariables][[1]], i++,
      variable = fourthOrderVariables[[i]];
      If[variable == (variable /. fourthOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 4, j++,
            For[k = 0, k <= maxPower - 4, k++,
               For[l = 0, l <= maxPower - 4, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 4,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[fifthOrderVariables][[1]], i++,
      variable = fifthOrderVariables[[i]];
      If[variable == (variable /. fifthOrderSolutions)[[1]],
         For[j = 0, j <= maxPower - 5, j++,
            For[k = 0, k <= maxPower - 5, k++,
               For[l = 0, l <= maxPower - 5, l++,
               nonSymmetricPowers = {j, k, l};
               If[Total[nonSymmetricPowers] == maxPower - 5,
                  multiModeVector = ConstantArray[0, {numberOfRigidModes}];
                  powers = ToExpression[Characters[variable[[2]]]];
                  For[m = 1, m <= numberOfRigidModes, m++,
                     If[powers[[m]] > 0, 
                        If[m < numberOfRigidModes - 1, 
                           multiModeVector[[m]] = 1, 
                           multiModeVector[[m]] = 2];, 
                           multiModeVector[[m]] = 0];
   	              ];
                  multiMode = Total[multiModeVector[[1 ;; numberOfRigidModes - 2]]] + If[Total[multiModeVector[[numberOfRigidModes - 1 ;; numberOfRigidModes]]] > 0, 1, 0];  
                  multiMode = multiMode + If[j > 0, 1, 0] + If[k > 0, 1, 0] + If[l > 0, 1, 0];
                  If[multiMode <= maxMultiMode,
                      Print["* ", variable[[1]], "   ", j, " ", k, " ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " ", l, " ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
                      variable;
                     ];      
                  ];
               ];
            ];
         ];
         variable;
      ];
   ];
   For[i = 0, i <= Dimensions[sixthOrderVariables][[1]], i++,
      variable = sixthOrderVariables[[i]];
      If[variable == (variable /. sixthOrderSolutions)[[1]], 
         Print["* ", variable[[1]], "   0 0 ", StringTake[variable[[2]], {1, 1}], " ", StringTake[variable[[2]], {2, 2}], " ", StringTake[variable[[2]], {3, 3}], " 0 ", StringTake[variable[[2]], {4, 4}], " ", StringTake[variable[[2]], {5, 5}], " ", StringTake[variable[[2]], {6, 6}], " ", StringTake[variable[[2]], {7, 7}], " ", StringTake[variable[[2]], {8, 8}], "  ", variable[[3]] + fourierCycle*3];,
         variable;
      ];
   ];
];