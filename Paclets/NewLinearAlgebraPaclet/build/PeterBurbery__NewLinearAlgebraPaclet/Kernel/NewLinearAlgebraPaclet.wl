(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["PeterBurbery`NewLinearAlgebraPaclet`"];


(* ::Text:: *)
(*Declare your public symbols here:*)


PeterBurbery`NewLinearAlgebraPaclet`Antidiagonal;

PeterBurbery`NewLinearAlgebraPaclet`AntidiagonalMatrix;

PeterBurbery`NewLinearAlgebraPaclet`DesymmetrizedMatrix;

PeterBurbery`NewLinearAlgebraPaclet`DeTriangularizableMatrixQ;

PeterBurbery`NewLinearAlgebraPaclet`DeTriangularizeMatrix;

PeterBurbery`NewLinearAlgebraPaclet`MatrixSymmetrizability;

PeterBurbery`NewLinearAlgebraPaclet`PyramidMatrix;

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

ClearAll[AntidiagonalMatrix]
AntidiagonalMatrix::usage="AntidiagonalMatrix[list] gives a matrix with the elements of list on the leading antidiagonal, and 0 elsewhere.\nAntidiagonalMatrix[list,k] gives a matrix with the elements of list on the kth antidiagonal.\nAntidiagonalMatrix[list,k,n] pads with 0s to create an n by n matrix.";

AntidiagonalMatrix[list:(_List|_SparseArray)]:=AntidiagonalMatrix[list,0]
AntidiagonalMatrix[list:(_List|_SparseArray),k_Integer]:=Reverse[DiagonalMatrix[list,k],2]
AntidiagonalMatrix[list:(_List|_SparseArray),k_Integer,n_Integer]:=AntidiagonalMatrix[list,k,{n,n}]
AntidiagonalMatrix[list:(_List|_SparseArray),k_Integer,mn:{m_Integer,n_Integer}]:=Reverse[DiagonalMatrix[list,k,mn],2]

DeTriangularizableMatrixQ//ClearAll

DeTriangularizableMatrixQ[nonMatrix_]:=False

DeTriangularizableMatrixQ[matrix_?(LowerTriangularMatrixQ[#, (*the main diagonal*) 0] &)]:=True

DeTriangularizableMatrixQ[matrix_?(UpperTriangularMatrixQ[#,(*the main diagonal*) 0]&)]:=True

DeTriangularizableMatrixQ::usage="DeTriangularizableMatrixQ[matrix] gives True if matrix is a lower triangular matrix or an upper triangular matrix and False otherwise. The function will return True if the matrix is detriangularizable, and False otherwise.";

DeTriangularizeMatrix // ClearAll

DeTriangularizeMatrix::usage = 
  "DeTriangularizeMatrix[matrix] detriangularizes the upper triangular \
or lower triangular matrix matrix into a symmetric matrix.";

DeTriangularizeMatrix[matrix_?(LowerTriangularMatrixQ[#] &)] := 
 LowerTriangularize[matrix, -1] + UpperTriangularize[Transpose[matrix]]

DeTriangularizeMatrix[matrix_?(UpperTriangularMatrixQ[#] &)] := 
 UpperTriangularize[matrix, 1] + LowerTriangularize[Transpose[matrix]]

MatrixSymmetrizability//ClearAll

MatrixSymmetrizability[matrix_?SquareMatrixQ]:=Mean[Boole/@MapThread[Equal,{Flatten[MapThread[Take[#1,#2]&,{Rest@Transpose[UpperTriangularize[matrix,1]],Range[First[Dimensions[UpperTriangularize[matrix,1]]]-1]}]],Flatten[MapThread[Take[#1,#2]&,{Rest@LowerTriangularize[matrix,1],Range[First[Dimensions[LowerTriangularize[matrix,1]]]-1]}]]}]]

MatrixSymmetrizability::usage="MatrixSymmetrizability[matrix] returns 1 if matrix is completely symmetric and 0 if matrix has no symmetry other than on the main diagonal. The closer the value is to 1, the more symmetralizable the matrix is.";

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
