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

typedef Basis<T,Primal,R,SparseMulti> Basis1D;

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

//Operator definitions
//typedef HelmholtzOperator1D<T, Basis1D>                             PDEOp1D;
typedef IdentityOperator1D<T, Basis1D>                              PDEOp1D;
typedef NoCompression<T, Index1D, Basis1D>                          Compression1D;
//typedef CompressionPDE1D<T, Basis1D>                                Compression1D;
typedef DiagonalMatrixPreconditioner1D<T,Basis1D,PDEOp1D>     Preconditioner1D;

//MapMatrix definition
typedef MapMatrix<T,Index1D,PDEOp1D,
                  Compression1D,Preconditioner1D>               MA;

template <typename T>
void
computeEV(DenseMatrixT &A, T &cB, T &CB);

template <typename T>
void
computeSV(DenseMatrixT &A, T &cB, T &CB);

template<typename T>
IndexSet<Index1D>
LambdaForEigenvalues(int jmin, int jmax, const Basis1D &basis, T radius);

typedef IndexSet<Index1D>::const_iterator const_set_it;

int main (int argc, char *argv[]) {
    if (argc != 6) {
        cout << "usage " << argv[0] << " d c radius jmin max_level" << endl;
        exit(1);
    }

    int d=atoi(argv[1]);
    T   c=atof(argv[2]);
    T   radius    =atof(argv[3]);
    int jmin =atoi(argv[4]);
    int max_level =atoi(argv[5]);
    cout.precision(3);

    Basis1D             basis(d,jmin);
    //PDEOp1D           pde_op(basis,c);
    PDEOp1D             pde_op(basis);
    Preconditioner1D    prec(pde_op);
    Compression1D       compression(basis);
    MA                  A(pde_op,prec,compression);

    cout << "Chosen parameters: d=" << d << ", c=" << c << ", radius="<< radius
         << ", jmin=" << jmin << ", max_level=" << max_level << endl;

    std::stringstream filename;
    filename << "eigenvalues_mw_helmholtz_" << d << "_" << jmin << "_" << c <<  ".dat";
    std::ofstream file_eigenvalues(filename.str().c_str());
    for (int jmax=0; jmax<=max_level; jmax+=1) {
        for (T r=1.; r<=radius; r+=1.) {
            IndexSet<Index1D> Lambda = LambdaForEigenvalues(jmin, jmax, basis, r);
            int N = Lambda.size();
            cout << "Size of Lambda: " << N << endl;
            SparseMatrixT A_flens(N,N);
            toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
            DenseMatrixT A_dense;
            densify(cxxblas::NoTrans,A_flens,A_dense);

            cout << A_dense << endl;

            T cB, CB;
            computeEV(A_dense, cB, CB);
            file_eigenvalues << " " << jmax << " " << r
                             << " " << cB << " " << " " << CB << endl;
            cout             << " " << jmax << " " << r
                             << " " << cB << " " << " " << CB << endl;
        }
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

template<typename T>
IndexSet<Index1D>
LambdaForEigenvalues(int jmin, int jmax, const Basis1D &basis, T radius)
{
    IndexSet<Index1D> Lambda;
    int k_left, k_right;
    int numWavelets = (int)basis.psi._numSplines;
    int numScaling = (int)basis.mra.phi._numSplines;

    int wavelet_count = 0;

    for (int j=jmin; j<=jmax; ++j) {
        k_left = 1;//std::floor(-pow2i<T>(j)*radius-basis.psi.max_support().l2);
        k_right = 4;//std::ceil(pow2i<T>(j)*radius-basis.psi.max_support().l1);
        for (int k_help=k_left; k_help<=k_right; ++k_help) {
            for (int k=(k_help-1)*numWavelets; k<=k_help*numWavelets-1; ++k) {
                Lambda.insert(Index1D(j,k,XWavelet));
                ++wavelet_count;
            }
        }
    }
/*
    int scaling_count = 0;
    k_left  = int(std::floor(-pow2i<T>(jmin)*radius-basis.mra.phi.max_support().l2));
    k_right = int(std::ceil(  pow2i<T>(jmin)*radius-basis.mra.phi.max_support().l1));
    for (int k_help=k_left; k_help<=k_right; ++k_help) {
        for (int k=(k_help-1)*numScaling+1; k<=k_help*numScaling; ++k) {
            Lambda.insert(Index1D(jmin,k,XBSpline));
            ++scaling_count;
        }
    }
    cout << "   -> Current index set: " << scaling_count << " scaling indices, "
         << wavelet_count << " wavelet indices." << endl;
*/
    return Lambda;
}
