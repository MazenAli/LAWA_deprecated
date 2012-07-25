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

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >            SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >          DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                                DenseVectorT;

//Iterator definitions
typedef IndexSet<Index1D>::const_iterator                                   const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator             const_coeff1d_it;

//Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>                                    Dij_Basis1D;
typedef Basis<T,Orthogonal,Interval,Multi>                                  MW_Basis1D;
typedef Basis<T,Primal,Interval,SparseMulti>                                SparseMW_Basis1D;

//Operator definitions
typedef CompressionPDE1D<T, Dij_Basis1D>                                    Dij_CompressionPDE1D;

typedef IdentityOperator1D<T, Dij_Basis1D>                                  Dij_IdentityOp1D;
typedef DiagonalMatrixPreconditioner1D<T, Dij_Basis1D, Dij_IdentityOp1D >   Dij_IdentityPrec1D;
typedef MapMatrix<T, Index1D, Dij_IdentityOp1D, Dij_CompressionPDE1D,
                  Dij_IdentityPrec1D>                                       Dij_MA_L2;
typedef AdaptiveIdentityOperator1D<T,Orthogonal,Interval,Multi>             MW_MA_L2;
typedef AdaptiveIdentityOperator1D<T,Primal,Interval,SparseMulti>           SparseMW_MA_L2;


typedef HelmholtzOperator1D<T, Dij_Basis1D>                                 Dij_HelmholtzOp1D;
typedef DiagonalMatrixPreconditioner1D<T, Dij_Basis1D, Dij_HelmholtzOp1D >  Dij_HelmholtzPrec1D;
typedef MapMatrix<T, Index1D, Dij_HelmholtzOp1D, Dij_CompressionPDE1D,
                  Dij_HelmholtzPrec1D>                                      Dij_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Interval,Multi>   MW_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,Interval,SparseMulti> SparseMW_MA;


template <typename Basis>

IndexSet<Index1D>
LambdaForEigenvalues(const Basis &basis, int jmin, int jmax);

void
computelargestEV(const SparseMatrixT &A, T &CB);

void
computesmallestEV(const SparseMatrixT &A, T &cB);

void
computeEV(DenseMatrixT &A, T &cB, T &CB);

void
computeSV(DenseMatrixT &A, T &cB, T &CB);

int main (int argc, char *argv[]) {
    if (argc!=7) {
        cout << "usage " << argv[0] << " basistype d d_ jmin s c" << endl; exit(1);
    }
    cout.precision(10);

    int d =atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0=atoi(argv[4]);
    int s =atoi(argv[5]);
    T   c =atof(argv[6]);

    int max_level  = 10;

    std::stringstream filename;
    if (s==0) {
        filename << "eigenvalues_L2_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                                      << argv[4] << ".dat";
    }
    else if (s==1) {
        filename << "eigenvalues_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                                             << argv[4] << "_" << argv[6] << ".dat";
    }
    else {
        std::cerr << "Not implemented for s=" << s << std::endl;
        return 0;
    }
    std::ofstream file_eigenvalues(filename.str().c_str());
    file_eigenvalues.precision(16);

    if (strcmp(argv[1],"Dijkema")==0) {
        Dij_Basis1D             Dij_basis(d,d_,j0);
        Dij_basis.enforceBoundaryCondition<DirichletBC>();
        Dij_CompressionPDE1D    Dij_compression1D(Dij_basis);

        for (int jmax=j0; jmax<=max_level; jmax+=1) {
            IndexSet<Index1D> Lambda = LambdaForEigenvalues(Dij_basis, j0, jmax);
            int N = Lambda.size();
            cout << "Size of Lambda: " << N << endl;
            SparseMatrixT A(N,N);
            if (s==0) {
                Dij_IdentityOp1D        Dij_identityOp1D(Dij_basis);
                Dij_IdentityPrec1D      Dij_prec1D(Dij_identityOp1D);   // wavelets are not L2-normalized by default!!
                Dij_MA_L2               Dij_A(Dij_identityOp1D, Dij_prec1D, Dij_compression1D);
                Dij_A.toFlensSparseMatrix(Lambda, Lambda, A);
            }
            else if (s==1) {
                Dij_HelmholtzOp1D       Dij_helmholtzOp1D(Dij_basis,c);
                Dij_HelmholtzPrec1D     Dij_prec1D(Dij_helmholtzOp1D);
                Dij_MA                  Dij_A(Dij_helmholtzOp1D, Dij_prec1D, Dij_compression1D);
                Dij_A.toFlensSparseMatrix(Lambda, Lambda, A);
            }

            T cB=0., CB=0.;
            //DenseMatrixT A_dense;
            //densify(cxxblas::NoTrans,A,A_dense);
            //computeEV(A_dense, cB, CB);

            T cB2, CB2;
            computesmallestEV(A,cB2);
            computelargestEV(A,CB2);

            file_eigenvalues << " " << jmax << " "
                             << " " << cB << " " << cB2 << " " << CB  << " " << CB2 << endl;
            cout             << " " << jmax << " "
                             << " " << cB << " " << cB2 << " " << CB  << " " << CB2 << endl;
        }
        file_eigenvalues << endl;
    }
    else if (strcmp(argv[1],"MW")==0) {
        MW_Basis1D             MW_basis(d,j0);
        MW_basis.enforceBoundaryCondition<DirichletBC>();
        MW_MA                  MW_A(MW_basis,c);
        MW_MA_L2               MW_A_L2(MW_basis);

        for (int jmax=std::max(0,MW_basis.j0); jmax<=max_level; jmax+=1) {
            IndexSet<Index1D> Lambda = LambdaForEigenvalues(MW_basis, j0, jmax);
            int N = Lambda.size();
            cout << "Size of Lambda: " << N << endl;
            SparseMatrixT A(N,N);
            if (s==0) {
                MW_A_L2.toFlensSparseMatrix(Lambda, Lambda, A);
            }
            else if (s==1) {
                MW_A.toFlensSparseMatrix(Lambda, Lambda, A);
            }

            T cB=0., CB=0.;
            //DenseMatrixT A_dense;
            //densify(cxxblas::NoTrans,A,A_dense);
            //computeEV(A_dense, cB, CB);

            T cB2, CB2;
            computesmallestEV(A,cB2);
            computelargestEV(A,CB2);

            file_eigenvalues << " " << jmax << " "
                             << " " << cB << " " << cB2 << " " << CB  << " " << CB2 << endl;
            cout             << " " << jmax << " "
                             << " " << cB << " " << cB2 << " " << CB  << " " << CB2 << endl;
        }
        file_eigenvalues << endl;
    }
    else if (strcmp(argv[1],"SparseMW")==0) {
        cout << "   SparseMW realline" << endl;
        SparseMW_Basis1D       SparseMW_basis(d,j0);
        SparseMW_basis.enforceBoundaryCondition<DirichletBC>();
        SparseMW_MA            SparseMW_A(SparseMW_basis,c);
        SparseMW_MA_L2         SparseMW_A_L2(SparseMW_basis);
        for (int jmax=std::max(0,SparseMW_basis.j0); jmax<=max_level; jmax+=1) {

            IndexSet<Index1D> Lambda = LambdaForEigenvalues(SparseMW_basis, j0, jmax);

            if (jmax>=12) continue;

            int N = Lambda.size();
            //cout << "Lambda = " << Lambda << endl;
            cout << "Size of Lambda: " << N << endl;
            SparseMatrixT A(N,N);
            if (s==0) {
                SparseMW_A_L2.toFlensSparseMatrix(Lambda, Lambda, A);
            }
            else if (s==1) {
                SparseMW_A.toFlensSparseMatrix(Lambda, Lambda, A);
            }

            T cB=0., CB=0.;
            //DenseMatrixT A_dense;
            //densify(cxxblas::NoTrans,A,A_dense);
            //cout << "A = " << A_dense << endl;
            //computeEV(A_dense, cB, CB);

            T cB2, CB2;
            computesmallestEV(A,cB2);
            computelargestEV(A,CB2);

            file_eigenvalues << " " << jmax << " "
                             << " " << cB << " " << cB2 << " " << CB  << " " << CB2 << endl;
            cout             << " " << jmax << " "
                             << " " << cB << " " << cB2 << " " << CB  << " " << CB2 << endl;

        }
        file_eigenvalues << endl;
    }

    return 0;

}



template <typename Basis>
IndexSet<Index1D>
LambdaForEigenvalues(const Basis &basis, int jmin, int jmax)
{
    IndexSet<Index1D> Lambda;
    for (int j=jmin; j<=jmax; ++j) {
        for (long k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            Lambda.insert(Index1D(j,k,XWavelet));
        }
    }
    for (long k=basis.mra.rangeI(jmin).firstIndex(); k<=basis.mra.rangeI(jmin).lastIndex(); ++k) {
        Lambda.insert(Index1D(jmin,k,XBSpline));
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
