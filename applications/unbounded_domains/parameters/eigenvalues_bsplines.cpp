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

typedef Basis<T,Primal,R,CDF> Basis1D;

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;


//Operator definitions
typedef PDEOperator1D<T, Basis1D>          PDEOp;
typedef CompressionPDE1D<T, Basis1D>       Compression;
typedef NoPreconditioner<T,Index1D>        Preconditioner;

//MapMatrix definition
typedef MapMatrix<T,Index1D,PDEOp, Compression,Preconditioner>   MA;

template <typename T>
void
computeEV(DenseMatrixT &A, T &cB, T &CB);

template <typename T>
void
computeSV(DenseMatrixT &A, T &cB, T &CB);

IndexSet<Index1D>
LambdaForEigenvalues(const Basis1D &basis, int k_left, int k_right);

typedef IndexSet<Index1D>::const_iterator const_set_it;

int main (int argc, char *argv[]) {
    if (argc != 5) {
        cout << "usage " << argv[0] << " d jmin maxK H1" << endl;
        exit(1);
    }

    int d=     atoi(argv[1]);
    int jmin = atoi(argv[2]);
    int maxK = atoi(argv[3]);
    int H1   = atoi(argv[4]);
    cout.precision(16);
    T reaction = 1., convection=0., diffusion=0.;
    if (H1==1)                      diffusion=1.;

    Basis1D         basis(d,d,jmin);
    PDEOp           pdeop(basis,reaction,convection,diffusion);
    Preconditioner  prec;
    Compression     compression(basis);
    MA              A(pdeop,prec,compression);

    cout << "Chosen parameters: d=" << d << ", jmin=" << jmin << ", maxK=" << maxK
                                         << ", H1=" << H1 << endl;

    std::stringstream filename;
    filename << "eigenvalues_BSplines_in_H" << H1 << "_" << jmin << "_" << d <<  ".dat";
    std::ofstream file_eigenvalues(filename.str().c_str());
    file_eigenvalues.precision(16);

    for (int k=0; k<=maxK; k+=50) {
        IndexSet<Index1D> Lambda = LambdaForEigenvalues(basis,-k,k);
        int N = Lambda.size();
        cout << "Size of Lambda: " << N << endl;
        SparseMatrixT A_flens(N,N);
        toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
        DenseMatrixT A_dense;
        densify(cxxblas::NoTrans,A_flens,A_dense);

        T cB, CB;
        computeEV(A_dense, cB, CB);
        file_eigenvalues << " " << Lambda.size()
                         << " " << cB << " " << " " << CB << " " << CB/cB
                         << " " << sqrt(cB) << " " << sqrt(CB) << " " << sqrt(CB/cB) << endl;
        cout             << " " << Lambda.size()
                         << " " << cB << " " << " " << CB << " " << CB/cB
                         << " " << sqrt(cB) << " " << sqrt(CB) << " " << sqrt(CB/cB) << endl;
    }

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

IndexSet<Index1D>
LambdaForEigenvalues(const Basis1D &basis, int k_left, int k_right)
{
    IndexSet<Index1D> Lambda;
    for (int k=k_left; k<=k_right; ++k) {
        Lambda.insert(Index1D(basis.j0,k,XBSpline));
    }
    return Lambda;

}
