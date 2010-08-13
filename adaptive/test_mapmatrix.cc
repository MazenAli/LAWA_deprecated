/*
 * test_mapmatrix.cc
 *
 *  Created on: 12.08.2010
 *      Author: sebastian
 */

#include <iostream>
#include <adaptive/mapmatrix.h>
#include <lawa/lawa.h>

typedef double T;

using namespace lawa;
using namespace std;

typedef Basis<T,Primal,Interval,Dijkema> MyBasis;
typedef HelmholtzOperator1d<T,MyBasis> MyBilinearForm;
typedef Preconditioner<T,Index1d,MyBasis, MyBilinearForm > MyPreconditioner;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > FullColMatrixT;

int main() {
	int d=2, d_=4;
	MyBilinearForm bil(0.5);
	MapMatrix<T,Index1d,MyBilinearForm,MyBilinearForm,MyPreconditioner> A(bil);
	Index1d index1(2,2,XBSpline), index2(2,2,XWavelet);
	cout << A(index1,index2) << endl;

	Coefficients<Lexicographical,T,Index1d> res(d,d_), v(d,d_);
	v[index1] = 2.; v[index2] = 3.;
	IndexSet<Index1d> LambdaRow(d,d_);
	LambdaRow.insert(index1); LambdaRow.insert(index2); LambdaRow.insert(Index1d(3,4,XWavelet));
	res = mv(LambdaRow, A, v);
	cout << "v = " << v << endl;
	cout << "Ab = " << res << endl;

	flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A_flens;
	A_flens = A.toFlensSparseMatrix(LambdaRow,LambdaRow);
	FullColMatrixT A_dense;
	densify(NoTrans,A_flens,A_dense);
	cout << A_dense << endl;

	Preconditioner<T,Index1d,MyBasis,MyBilinearForm> prec;

	return 0;
}
