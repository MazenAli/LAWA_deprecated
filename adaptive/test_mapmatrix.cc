/*
 * test_mapmatrix.cc
 *
 *  Created on: 12.08.2010
 *      Author: sebastian
 */

#include <iostream>
#include <adaptive/mapmatrix.h>
#include <adaptive/referencesolutions.h>
#include <adaptive/problem.h>
#include <adaptive/rhs.h>
#include <lawa/lawa.h>

typedef double T;

using namespace lawa;
using namespace std;

typedef Basis<T,Primal,R,CDF> MyBasis;
typedef HelmholtzOperator1d<T,MyBasis> MyBilinearForm;
typedef Preconditioner<T,Index1d,MyBasis, MyBilinearForm > MyPreconditioner;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > FullColMatrixT;


int main() {
	int d=2, d_=2;
	T c = 0.;
	MyBasis basis(d,d_,0);
	MyBilinearForm bil(basis,c);
	MapMatrix<T,Index1d,MyBilinearForm,MyBilinearForm,MyPreconditioner> A(bil);
	Index1d index1(0,0,XBSpline), index2(0,0,XWavelet);
	Wavelet<T,Primal,R,CDF> d_psi1(d,d_,1), d_psi2(basis,1);
	Integral<T, Gauss, Wavelet<T,Primal,R,CDF>, Wavelet<T,Primal,R,CDF> > dd_integral_ww1(d_psi1,d_psi1), dd_integral_ww2(d_psi2,d_psi2);
	cout << "d=" << basis.d  << ", d_=" << basis.d_ << endl;
	cout << "wavlet.deriv=" << d_psi2.deriv << endl;
	cout << dd_integral_ww1(0,0,0,0) << " " << dd_integral_ww2(0,0,0,0) << " " << A(index2,index2) << endl;

	Coefficients<Lexicographical,T,Index1d> res(d,d_), v(d,d_);
	v[index1] = 2.; v[index2] = 3.;
	IndexSet<Index1d> LambdaRow(d,d_);
	LambdaRow.insert(index1); LambdaRow.insert(index2); LambdaRow.insert(Index1d(3,4,XWavelet));
	res = mv(LambdaRow, A, v);
	cout << "v = " << v << endl;
	cout << "Av = " << res << endl;


	flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A_flens;
	A_flens = A.toFlensSparseMatrix(LambdaRow,LambdaRow);
	FullColMatrixT A_dense;
	densify(NoTrans,A_flens,A_dense);
	cout << A_dense << endl;

	Preconditioner<T,Index1d,MyBasis,MyBilinearForm> prec;
	ReferenceSolution1d<T,MyBasis,MyBilinearForm> refsol;
	refsol.setExample(2, bil, R);

	cout << refsol.sing_pts << endl;
	cout << refsol.deltas << endl;

	Problem1d<T,MyBasis,MyBilinearForm>(basis, 1, bil, R);

	RHS<T,Index1d,MyBasis,MyPreconditioner> F(basis,refsol.rhs,refsol.sing_pts, refsol.deltas);
	Coefficients<Lexicographical,T,Index1d> f_Lambda = F(LambdaRow);
	cout << f_Lambda << endl;


	return 0;
}
