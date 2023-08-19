(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["PeterBurbery`NewLinearAlgebraPaclet`"];


(* ::Text:: *)
(*Declare your public symbols here:*)


PeterBurbery`NewLinearAlgebraPaclet`UlamMatrix;

PeterBurbery`NewLinearAlgebraPaclet`Antidiagonal;


Begin["`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Text:: *)
(*Define your public and private symbols here:*)


UlamMatrix//ClearAll

UlamMatrix::usage="UlamMatrix[n] returns the n by n Ulam matrix.";

Options[UlamMatrix]={WorkingPrecision->Infinity};

UlamMatrix[n_Integer?Positive,OptionsPattern[]]:=Module[{d,i,j,k,q,res},
res=ConstantArray[0,{n,n}];
j=Ceiling[n/2];
i=j+Mod[n+1,2];
res[[i,j]]=1;
d=k=1;
Do[
q=Range[Min[p,n-1]];
j+=d q;k+=q;
res[[i,j]]=k;
If[p==n,Break[]];
j=j[[p]];k=k[[p]];
i-=d q;k+=q;
res[[i,j]]=k;
i=i[[p]];k=k[[p]];
d=-d,
{p,n}];
N[res,OptionValue[WorkingPrecision]]]

ClearAll[Antidiagonal]

Antidiagonal::usage="Antidiagonal[m] gives the list of elements on the leading antidiagonal of the matrix m.\nAntidiagonal[m, k] gives the elements on the kth antidiagonal of m.";

Antidiagonal[m_?MatrixQ]:=Diagonal[Reverse[m,2]]
Antidiagonal[m_?MatrixQ,k_Integer]:=Diagonal[Reverse[m,2],k]


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
