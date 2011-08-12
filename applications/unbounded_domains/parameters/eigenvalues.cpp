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
#include <unistd.h>
#include <lawa/lawa.h>
#include <applications/unbounded_domains/referencesolutions/referencesolutions.h>

#define ROW_SIZE 4*8192
#define COL_SIZE 4*2048

typedef double T;
using namespace lawa;
using namespace std;

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >        SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >      DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                            DenseVectorT;

//Iterator definitions
typedef IndexSet<Index1D>::const_iterator                               const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator         const_coeff1d_it;
typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator           const_coeff1d_abs_it;

//Basis definitions
typedef Basis<T,Primal,R,CDF>                                           CDF_Basis1D;
typedef Basis<T,Orthogonal,R,Multi>                                     MW_Basis1D;
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;

//Operator definitions
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>            CDF_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,R,Multi>      MW_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,SparseMulti>    SparseMW_MA;


IndexSet<Index1D>
LambdaForEigenvalues(const CDF_Basis1D &basis, int jmin, int jmax, bool w_XBSpline, T radius);

void
computelargestEV(const SparseMatrixT &A, T &CB);

void
computesmallestEV(const SparseMatrixT &A, T &cB);

void
computeEV(DenseMatrixT &A, T &cB, T &CB);

void
computeSV(DenseMatrixT &A, T &cB, T &CB);

int main (int argc, char *argv[]) {
    if (argc!=6) {
        cout << "usage " << argv[0] << " basistype d d_ jmin c" << endl; exit(1);
    }
    cout.precision(10);

    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    T   c=atof(argv[5]);
    int j0; bool w_XBSpline;

    if   (strcmp(argv[4],"-inf")==0) { j0=0;             w_XBSpline=false; }
    else                             { j0=atoi(argv[4]); w_XBSpline=true; }

    int max_level  = 10;
    T   max_radius = 10.;

    std::stringstream filename;
    filename << "eigenvalues_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                                         << argv[4] << "_" << argv[5] << ".dat";
    std::ofstream file_eigenvalues(filename.str().c_str());
    file_eigenvalues.precision(16);

    /*
    DenseVectorT x(5);
    x = 4., 1., -2., 3., 5.;
    int i;
    cxxblas::iamax<int>((int)x.length(),x.engine().data(),(int)1,i);
    cout << "i = " << i << endl;
    */

    if (strcmp(argv[1],"CDF")==0) {
        CDF_Basis1D             CDF_basis(d,d_,j0);
        CDF_MA                  CDF_A(CDF_basis,w_XBSpline,c);
        if (w_XBSpline) {
            for (int jmax=std::max(0,CDF_basis.j0); jmax<=max_level; jmax+=1) {
                for (T r=0.25; r<=max_radius; r+=1.) {
                    IndexSet<Index1D> Lambda = LambdaForEigenvalues(CDF_basis, CDF_basis.j0, jmax,
                                                                    w_XBSpline, r);
                    int N = Lambda.size();
                    cout << "Size of Lambda: " << N << endl;
                    SparseMatrixT A(N,N);
                    CDF_A.toFlensSparseMatrix(Lambda, Lambda, A);
                    DenseMatrixT A_dense;
                    densify(cxxblas::NoTrans,A,A_dense);

                    T cB, CB;
                    computeEV(A_dense, cB, CB);
                    T cB2, CB2;
                    computesmallestEV(A,cB2);
                    computelargestEV(A,CB2);

                    file_eigenvalues << " " << jmax << " " << r
                                     << " " << cB << " " << " " << CB << endl;
                    cout             << " " << jmax << " " << r
                                     << " " << cB << " " << cB2 << " " << CB  << " " << CB2 << endl;
                }
                file_eigenvalues << endl;
            }
        }
        else {
            for (int jmax=0; jmax<=max_level; jmax+=1) {
                for (T r=1.; r<=max_radius; r+=1.) {
                    IndexSet<Index1D> Lambda = LambdaForEigenvalues(CDF_basis, -100, jmax,
                                                                    w_XBSpline, r);
                    int N = Lambda.size();
                    cout << "Size of Lambda: " << N << endl;
                    SparseMatrixT A(N,N);
                    CDF_A.toFlensSparseMatrix(Lambda, Lambda, A);
                    DenseMatrixT A_dense;
                    densify(cxxblas::NoTrans,A,A_dense);

                    T cB, CB;
                    computeEV(A_dense, cB, CB);

                    file_eigenvalues << " " << jmax << " " << r
                                     << " " << cB << " " << " " << CB << endl;
                    cout             << " " << jmax << " " << r
                                     << " " << cB << " " << " " << CB << endl;

                }
            }
        }
    }
    else if (strcmp(argv[1],"MW")==0) {
        MW_Basis1D             MW_basis(d,j0);
        MW_MA                  MW_A(MW_basis,w_XBSpline,c);
    }
    else if (strcmp(argv[1],"SparseMW")==0) {
        SparseMW_Basis1D       SparseMW_basis(d,j0);
        SparseMW_MA            SparseMW_A(SparseMW_basis,w_XBSpline,c);

    }

    return 0;

}

IndexSet<Index1D>
LambdaForEigenvalues(const CDF_Basis1D &basis, int jmin, int jmax, bool w_XBSpline, T r)
{
    const BSpline<T,Primal,R,CDF> phi = basis.mra.phi;
    const Wavelet<T,Primal,R,CDF> psi = basis.psi;
    IndexSet<Index1D> Lambda;
    int k_left, k_right;
    if (w_XBSpline) {
        for (int j=jmin; j<=jmax; ++j) {
            k_left = std::floor(-pow2i<T>(j)*r-psi.support(0,0).l2);
            k_right =std::ceil(pow2i<T>(j)*r-psi.support(0,0).l1);
            for (int k=k_left; k<=k_right; ++k) {
                Lambda.insert(Index1D(j,k,XWavelet));
            }
        }

        k_left  = std::min(int(std::floor(-pow2i<T>(jmin)*r-phi.support(0,0).l2)),-20);
        k_right = std::min(int(std::ceil(  pow2i<T>(jmin)*r-phi.support(0,0).l1)),20);
        for (int k=k_left; k<=k_right; ++k) {
            Lambda.insert(Index1D(jmin,k,XBSpline));
        }
    }
    else {
        for (int j=jmin; j<=jmax; ++j) {
            int k_left, k_right;
            if (j>=-6) {
                k_left = std::floor(-pow2i<T>(j)*r-psi.support(0,0).l2);
                k_right = std::ceil(pow2i<T>(j)*r-psi.support(0,0).l1);
                if (k_left>=0) cout << "j=" << j << ", k_left=" << k_left << endl;
                if (k_right<=0) cout << "j=" << j << ", k_right=" << k_right << endl;
            }
            else {
                k_left = -psi.d-psi.d_;
                k_right = psi.d+psi.d_;
            }
            for (int k=k_left; k<=k_right; ++k) {
                Lambda.insert(Index1D(j,k,XWavelet));
            }
        }
    }
    return Lambda;
}

void
computelargestEV(const SparseMatrixT &A, T &CB)
{
    int N = A.numRows();
    DenseVectorT x(N);
    for (int i=1; i<=N; ++i) {
        x(i) = 1.;
    }
    cout << "powerMethod started." << endl;
    lawa::powerMethod(A,(T)1e-12,CB,x);
    cout << "powerMethod finished." << endl;
}

void
computesmallestEV(const SparseMatrixT &A, T &cB)
{
    int N = A.numRows();
    DenseVectorT x(N);
    for (int i=1; i<=N; ++i) {
        x(i) = 1.;
    }
    cout << "inversePowerMethod started." << endl;
    lawa::inversePowerMethod(A,(T)1e-12,cB,x);
    cout << "inversePowerMethod finished." << endl;
}

void
computeEV(DenseMatrixT &A, T &cB, T &CB) {
    DenseMatrixT U(A.numRows(),A.numRows()), V(A.numCols(),A.numCols());
    int N = A.numRows();
    DenseVectorT wr(N), wi(N), x(N);
    DenseMatrixT vl,vr;
    ev(false, false, A, wr, wi, vl, vr);
    T lambda;
    cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
    for (int i=1; i<=wr.lastIndex(); ++i) {
        cB = std::min(cB,wr(i));
        CB = std::max(CB,wr(i));
    }
}

void
computeSV(DenseMatrixT &A, T &cB, T &CB) {
    DenseMatrixT U(A.numRows(),A.numRows()), V(A.numCols(),A.numCols());
    DenseVectorT s(A.numCols());
    int iterations = svd(A,s,U,V);
    CB = s(s.firstIndex());
    cB = s(s.lastIndex());
}
