(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["PeterBurbery`LinearAlgebraPaclet`"];


(* ::Text:: *)
(*Declare your public symbols here:*)


ConsistentMatrixQ;CofactorMatrix;DiagonalizeMatrix;


Begin["`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Text:: *)
(*Define your public and private symbols here:*)


ConsistentMatrixQ // ClearAll
ConsistentMatrixQ[mat_?MatrixQ] := 
 MatrixRank[Map[Most][mat]] == MatrixRank[mat]
 
 CofactorMatrix//ClearAll
 CofactorMatrix[mat_?SquareMatrixQ]:=Transpose[Adjugate[mat]]
 
 DiagonalizeMatrix//ClearAll
 DiagonalizeMatrix[matrix_?DiagonalizableMatrixQ]:=Last[JordanDecomposition[matrix]]


(* ::Section:: *)
(*Package Footer*)


End[];
EndPackage[];
