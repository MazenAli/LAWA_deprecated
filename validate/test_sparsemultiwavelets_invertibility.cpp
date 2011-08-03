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
#include <vector>
#include <lawa/lawa.h>

typedef double T;
using namespace lawa;
using namespace std;

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef Basis<T,Primal,R,SparseMulti> Basis1D;

typedef IntegralF<Gauss, Basis1D>                                   IntegralF1D;

template <typename T>
void
computeEV(DenseMatrixT &A, T &cB, T &CB);

template <typename T>
void
computeSV(DenseMatrixT &A, T &cB, T &CB);

int shift = 1;

T
wtphi1(T x) {
    x -= (T)shift;
    if (0.<=x && x<=2.) {
        return 1.;
    }
    else return 0.;
}

T
wtphi2(T x) {
    x -= (T)shift;
    if (0.<=x && x<=2.) {
        return (x-1.);
    }
    else return 0.;
}

T
wtphi3(T x) {
    x -= (T)shift;
    if (0.<=x && x<=2.) {
        return (x-1.)*(x-1.);
    }
    else return 0.;
}

T
wtphi4(T x) {
    x -= (T)shift;
    if (0.<=x && x<=2.) {
        return (x-1.)*(x-1.)*(x-1.);
    }
    else return 0.;
}

int main (int argc, char *argv[]) {

    cout.precision(6);

    Basis1D             basis(4,0);

    DenseVectorT        singPts;
    Function<T> wtphi1Fct(wtphi1,singPts);
    Function<T> wtphi2Fct(wtphi2,singPts);
    Function<T> wtphi3Fct(wtphi3,singPts);
    Function<T> wtphi4Fct(wtphi4,singPts);

    IntegralF1D integral_wtphi1_phi(wtphi1Fct,basis);
    IntegralF1D integral_wtphi2_phi(wtphi2Fct,basis);
    IntegralF1D integral_wtphi3_phi(wtphi3Fct,basis);
    IntegralF1D integral_wtphi4_phi(wtphi4Fct,basis);

    int n=500;
    int N = n*4;
    DenseMatrixT A(N,N);
    cout << "N = " << N << endl;

    for (int k_row=0; k_row<=N-1; ++k_row) {
        for (int k_col=0; k_col<=N-1; ++k_col) {
            shift = 2*(k_col / 4);
            //cout << "phi_support = " <<  basis.mra.phi.support(0,k_row) << endl;
            //cout << "wtphi_support = " <<  Support<T>(0.+shift,2.+shift) << endl;
            T val = 0.;
            if (k_col % 4==0) {
                val = integral_wtphi1_phi(0,k_row,XBSpline,0);
            }
            else if (k_col % 4==1) {
                val = integral_wtphi2_phi(0,k_row,XBSpline,0);
            }
            else if (k_col % 4==2) {
                val = integral_wtphi3_phi(0,k_row,XBSpline,0);
            }
            else {
                val = integral_wtphi4_phi(0,k_row,XBSpline,0);
            }
            //cout << "val = " << val << endl << endl;
            A(k_row+1,k_col+1) = val;
        }
    }

    if (N<=16) {
        cout << A << endl;
    }
    T cB, CB;
    computeSV(A, cB, CB);
    cout << " " << cB << " " << " " << CB << endl;

    return 0;
}

template <typename T>
void
computeEV(DenseMatrixT &A, T &cB, T &CB) {

    DenseMatrixT U(A.numRows(),A.numRows()), V(A.numCols(),A.numCols());
    int N = A.numRows();
    DenseVectorT wr(N), wi(N);
    DenseMatrixT vl,vr;
    ev(false, false, A, wr, wi, vl, vr);
    cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
    for (int i=1; i<=wr.lastIndex(); ++i) {
        cB = std::min(cB,wr(i));
        CB = std::max(CB,wr(i));
    }
}

template <typename T>
void
computeSV(DenseMatrixT &A, T &cB, T &CB) {
    DenseMatrixT U(A.numRows(),A.numRows()), V(A.numCols(),A.numCols());
    DenseVectorT s(A.numCols());
    int iterations = svd(A,s,U,V);
    CB = s(s.firstIndex());
    cB = s(s.lastIndex());
}
