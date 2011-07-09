/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <iostream>
#include <lawa/lawa.h>
#include <applications/unbounded_domains/referencesolutions/referencesolutions.h>
#include <applications/unbounded_domains/parameters/parameters.h>

#define ROW_SIZE 4*8192
#define COL_SIZE 4*2048

typedef double T;
using namespace lawa;
using namespace std;

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >        SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >      DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                            DenseVectorT;

typedef IndexSet<Index1D>::const_iterator                               const_set1d_it;

//Basis definitions
typedef Basis<T,Orthogonal,R,Multi>                                     MWBasis1D;
typedef TensorBasis2D<Adaptive, MWBasis1D, MWBasis1D>                   MWBasis2D;

//Operator definitions
typedef HelmholtzOperator2D<T, MWBasis2D>                               HelmholtzOp2D;
typedef DiagonalMatrixPreconditioner2D<T,MWBasis2D, HelmholtzOp2D >     Preconditioner2D;

typedef AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D, Preconditioner2D>   MA_MW;
typedef AdaptiveHelmholtzOperator2D<T, MWBasis2D, Preconditioner2D>     MA;

typedef Parameters2D<T, MWBasis2D, HelmholtzOp2D>                       ParametersMW2D;

typedef APPLY_Helmholtz_MW_2D<T, MWBasis2D, ParametersMW2D, MA_MW>      APPLY2D;

int main(int argc, char *argv[]) {
    int d=2;
    int j0=-2;
    int J=atoi(argv[1]);
    MWBasis1D           mwbasis1d(d,j0);
    MWBasis2D           mwbasis2d(mwbasis1d,mwbasis1d);
    HelmholtzOp2D       helmholtz_op(mwbasis2d,1.);
    Preconditioner2D    prec(helmholtz_op);

    ParametersMW2D params(mwbasis2d, helmholtz_op);

    MA_MW               A_MW(mwbasis2d,1.,prec);
    MA                  A(mwbasis2d,1.,prec);

    APPLY2D Apply2d(mwbasis2d, params, A_MW);

    IndexSet<Index1D> Lambda1D;
    for (int k=-5; k<=5; ++k) {
        Lambda1D.insert(Index1D(j0,k,XBSpline));
    }

    for (int j=j0; j<=J; ++j) {
        for (int k=-5; k<=5; ++k) {
            Lambda1D.insert(Index1D(j,k,XWavelet));
        }
    }

    int N = Lambda1D.size() * Lambda1D.size();
    DenseVectorT x(N), y1(N), y2(N), diff(N);
    int count = 1;
    IndexSet<Index2D> Lambda;
    Coefficients<Lexicographical,T,Index2D> x_coeff, y1_coeff, y2_coeff, diff_coeff;
    for (const_set1d_it row_x=Lambda1D.begin(); row_x!=Lambda1D.end(); ++row_x) {
        for (const_set1d_it row_y=Lambda1D.begin(); row_y!=Lambda1D.end(); ++row_y) {
            T val = 1.;
            Lambda.insert(Index2D(*row_x,*row_y));
            x_coeff[Index2D(*row_x,*row_y)] = val;
            x(count) = val;
            ++count;
        }
    }


    Timer time;
    time.start();
    y1_coeff = A_MW.mv(Lambda, x_coeff);
    time.stop();
    cout << "Elapsed time for mv for Lamda.size()=" << Lambda.size() << ": " << time.elapsed() << endl;
    time.start();
    y1_coeff = A_MW.mv(Lambda, x_coeff);
    time.stop();
    cout << "Elapsed time for mv for Lamda.size()=" << Lambda.size() << ": " << time.elapsed() << endl;
    //cout << y1_coeff << endl;
    y2_coeff = mv(Lambda, A_MW, x_coeff);
    //cout << y2_coeff << endl;
    diff_coeff = y1_coeff - y2_coeff;

    cout << "ell_2 norm of difference: " << diff_coeff.norm(2.) << endl;

    time.start();
    SparseMatrixT       A_sparse1(Lambda.size(),Lambda.size());
    A_MW.toFlensSparseMatrix(Lambda, Lambda, A_sparse1,5);
    y1 = A_sparse1*x;
    time.stop();
    cout << "Elapsed time for assembly for Lamda.size()=" << Lambda.size() << ": " << time.elapsed() << endl;
    //DenseMatrixT A_dense1;
    //densify(cxxblas::NoTrans,A_sparse1,A_dense1);
    //cout << A_dense1 << endl;

    time.start();
    SparseMatrixT       A_sparse2(Lambda.size(),Lambda.size());
    toFlensSparseMatrix(A, Lambda, Lambda, A_sparse2);
    y2 = A_sparse2*x;
    time.stop();
    cout << "Standard approach: " << time.elapsed() << endl;
    //DenseMatrixT A_dense2;
    //densify(cxxblas::NoTrans,A_sparse2,A_dense2);
    //cout << A_dense2 << endl;


    diff = y1-y2;
    cout << "ell2-norm of error: " << sqrt(diff*diff) << endl;



    return 0.;
}
