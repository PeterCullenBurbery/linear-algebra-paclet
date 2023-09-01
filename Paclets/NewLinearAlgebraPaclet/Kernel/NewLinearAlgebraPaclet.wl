(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["PeterBurbery`NewLinearAlgebraPaclet`"];


(* ::Text:: *)
(*Declare your public symbols here:*)


PeterBurbery`NewLinearAlgebraPaclet`Antidiagonal;

PeterBurbery`NewLinearAlgebraPaclet`DeTriangularizableMatrixQ;

PeterBurbery`NewLinearAlgebraPaclet`DeTriangularizeMatrix;

PeterBurbery`NewLinearAlgebraPaclet`PyramidMatrix;

PeterBurbery`NewLinearAlgebraPaclet`DesymmetrizedMatrix;

PeterBurbery`NewLinearAlgebraPaclet`UlamMatrix;


Begin["`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Text:: *)
(*Define your public and private symbols here:*)


ClearAll[Antidiagonal]

Antidiagonal::usage="Antidiagonal[m] gives the list of elements on the leading antidiagonal of the matrix m.\nAntidiagonal[m, k] gives the elements on the kth antidiagonal of m.";

Antidiagonal[m_?MatrixQ]:=Diagonal[Reverse[m,2]]
Antidiagonal[m_?MatrixQ,k_Integer]:=Diagonal[Reverse[m,2],k]

DeTriangularizableMatrixQ//ClearAll

DeTriangularizableMatrixQ[nonMatrix_]:=False

DeTriangularizableMatrixQ[matrix_?(LowerTriangularMatrixQ[#, (*the main diagonal*) 0] &)]:=True

DeTriangularizableMatrixQ[matrix_?(UpperTriangularMatrixQ[#,(*the main diagonal*) 0]&)]:=True

DeTriangularizableMatrixQ::usage="DeTriangularizableMatrixQ[matrix] gives True if matrix is a lower triangular matrix or an upper triangular matrix and False otherwise. The function will return True if the matrix is detriangularizable, and False otherwise."

DeTriangularizeMatrix // ClearAll

DeTriangularizeMatrix::usage = 
  "DeTriangularizeMatrix[matrix] detriangularizes the upper triangular \
or lower triangular matrix matrix into a symmetric matrix.";

DeTriangularizeMatrix[matrix_?(LowerTriangularMatrixQ[#] &)] := 
 LowerTriangularize[matrix, -1] + UpperTriangularize[Transpose[matrix]]

DeTriangularizeMatrix[matrix_?(UpperTriangularMatrixQ[#] &)] := 
 UpperTriangularize[matrix, 1] + LowerTriangularize[Transpose[matrix]]

PyramidMatrix//ClearAll

PyramidMatrix::usage="PyramidMatrix[n] makes a pyramid matrix of size n by n.";

OddPyramidMatrix//ClearAll

OddPyramidMatrix[n_Integer?IntegerQ]/;n>=1:=Join[Join[(NestList[Join[Drop[DeleteDuplicates[#],-1],ConstantArray[DeleteDuplicates[#][[-2]],n-Length[DeleteDuplicates[#]]+1]]&,Range[n],n-1]//Reverse),Dataset[Transpose[NestList[Join[Drop[DeleteDuplicates[#],-1],ConstantArray[DeleteDuplicates[#][[-2]],n-Length[DeleteDuplicates[#]]+1]]&,Range[n],n-1]]][All,2;;]//Normal,2],Normal[Dataset[Join[(NestList[Join[Drop[DeleteDuplicates[#],-1],ConstantArray[DeleteDuplicates[#][[-2]],n-Length[DeleteDuplicates[#]]+1]]&,Range[n],n-1]//Reverse),Dataset[Transpose[NestList[Join[Drop[DeleteDuplicates[#],-1],ConstantArray[DeleteDuplicates[#][[-2]],n-Length[DeleteDuplicates[#]]+1]]&,Range[n],n-1]]][All,2;;]//Normal,2]][Reverse][2;;]]]

EvenPyramidMatrix//ClearAll

EvenPyramidMatrix[n_Integer?IntegerQ]/;n>=1:=Join[Join[(NestList[Join[Drop[DeleteDuplicates[#],-1],ConstantArray[DeleteDuplicates[#][[-2]],n-Length[DeleteDuplicates[#]]+1]]&,Range[n],n-1]//Reverse),Transpose[NestList[Join[Drop[DeleteDuplicates[#],-1],ConstantArray[DeleteDuplicates[#][[-2]],n-Length[DeleteDuplicates[#]]+1]]&,Range[n],n-1]],2],Join[NestList[Join[Drop[DeleteDuplicates[#],-1],ConstantArray[DeleteDuplicates[#][[-2]],n-Length[DeleteDuplicates[#]]+1]]&,Range[n],n-1],Transpose[NestList[Join[Drop[DeleteDuplicates[#],-1],ConstantArray[DeleteDuplicates[#][[-2]],n-Length[DeleteDuplicates[#]]+1]]&,Range[n],n-1]]//Reverse,2]]

PyramidMatrix[n_?OddQ]:=OddPyramidMatrix[Quotient[n,2]+1]

PyramidMatrix[n_?EvenQ]:=EvenPyramidMatrix[n/2]

DesymmetrizedMatrix // ClearAll

DesymmetrizedMatrix[matrix_?MatrixQ] := 
 SparseArray[matrix - Symmetrize[matrix]]

DesymmetrizedMatrix::usage = 
  "DesymmetrizedMatrix[matrix] gives a matrix where the parts of the matrix that stop it from being symmetric are left and the parts that are symmetric become 0, in effect returning a desymmetrized matrix by desymmetrizing the input matrix matrix.";

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


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
