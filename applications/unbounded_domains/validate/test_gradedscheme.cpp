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

typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >        SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >      DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                            DenseVectorT;

//Iterator definitions
typedef IndexSet<Index1D>::const_iterator                               const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator         const_coeff1d_it;
typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator           const_coeff1d_abs_it;

//Basis definitions
typedef Basis<T,Primal,R,CDF>                                           CDF_Basis1D;

//Operator definitions
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>            CDF_MA;

typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, CDF_MA>       CDF_Prec;

//Righthandsides definitions
typedef RHSWithPeaks1D<T, CDF_Basis1D>                                  CDF_RhsIntegral1D;

typedef RHS1D<T, CDF_RhsIntegral1D, CDF_Prec>                           CDF_Rhs;

void
getGradedLambda(IndexSet<Index1D> &Lambda, T gamma, int j0, int J);

int main (int argc, char *argv[]) {
    if (argc!=3) {
        cout << "usage " << argv[0] << " d d_ " << endl; exit(1);
    }
    cout.precision(8);


    int d=atoi(argv[1]);
    int d_=atoi(argv[2]);
    int j0=0;
    int J = 20;
    T gamma = 0.5;

    T c = 1.;
    int example=1;

    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,1.,0.,c);
    Function<T>        rhsFct(refsol.rhs,refsol.sing_pts);

    CDF_Basis1D             CDF_basis(d,d_,j0);
    CDF_MA                  CDF_A(CDF_basis,true,c);
    CDF_Prec                CDF_prec(CDF_A);

    CDF_RhsIntegral1D       CDF_rhsintegral1d(CDF_basis, rhsFct, refsol.deltas, 70);
    CDF_Rhs                 CDF_F(CDF_rhsintegral1d,CDF_prec);

    ofstream file("gradedscheme_conv.txt");
    for (int j=j0; j<=J; ++j) {
        cerr << "j = " << j << endl;
        IndexSet<Index1D> Lambda;
        getGradedLambda(Lambda,gamma,j0,j);
        int N = Lambda.size();
        cerr << "   GradedLambda has size " << N << endl;

        SparseMatrixT A(N,N);
        DenseVectorT  F(N), u(N);
        int row_count = 1;
        for (const_set1d_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
            F(row_count) = CDF_F(*row);
            u(row_count) = 0.;
            int col_count = 1;
            for (const_set1d_it col=Lambda.begin(); col!=Lambda.end(); ++col) {
                T val = CDF_A((*row),(*col));
                if (fabs(val)>0) A(row_count,col_count) = val;
                ++col_count;
            }
            ++row_count;
        }
        A.finalize();
        cerr << "   cg-solve started..." << endl;
        int numOfIterations = cg(A, u, F, 1e-18);
        cerr << "   ... finished." << endl;

        DenseVectorT  Au(N);
        Au = A*u;
        T uAu;
        uAu = Au*u;
        T energyerror = std::sqrt(fabs(refsol.H1norm()*refsol.H1norm() - uAu));
        cout << "N = " << N << ", " << energyerror << ", " << refsol.H1norm()*refsol.H1norm() << ", " << uAu << endl;
        file << j << " " <<  N << " " << energyerror << endl;
    }

    return 0;
}

void
getGradedLambda(IndexSet<Index1D> &Lambda, T gamma, int j0, int J)
{
    for (int j=j0; j<=J; ++j) {
        cerr << "   j = " << j << " " << 5*std::pow(2.,(T)j) << " " << 5* std::pow(2.,(T)(J-j)) << endl;
        for (int k=0; k<=3 + 3*std::pow(2.,(T)j) + 5*std::pow(2.,(T)(J-j)); ++k) {

            if (j==j0) {
                Lambda.insert(Index1D(j,k,XBSpline));
                Lambda.insert(Index1D(j,-k,XBSpline));
            }
            Lambda.insert(Index1D(j,k,XWavelet));
            Lambda.insert(Index1D(j,-k,XWavelet));
        }
    }
}
