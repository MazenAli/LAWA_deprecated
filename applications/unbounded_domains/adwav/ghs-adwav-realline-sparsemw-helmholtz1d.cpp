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
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMWBasis1D;

//Operator definitions
//typedef AdaptiveHelmholtzOperatorSparseMW1D<T, SparseMWBasis1D>         MW_MA;
typedef AdaptiveHelmholtzOperator1D<T,Primal,R,SparseMulti>             MW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, MW_MA>        MW_Prec;

//Righthandsides definitions (tensor)
typedef RHSWithPeaks1D<T, SparseMWBasis1D>                              RhsIntegral1D;
typedef RHS<T,Index1D, RhsIntegral1D, MW_Prec>                          Rhs;

typedef GHS_ADWAV<T, Index1D, MW_MA, Rhs>                               SparseMW_GHS_ADWAV_SOLVER;


IndexSet<Index1D>
computeRHSLambda_SmoothPart(const SparseMWBasis1D &basis, T a, T b, int J_plus);

template <typename T>
Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const SparseMWBasis1D &basis, const RhsIntegral1D &rhsintegral1d_singular,
                    const RhsIntegral1D &rhsintegral1d_smooth, MW_Prec &P,
                    RefSols_PDE_Realline1D<T> &refsol);


int main (int argc, char *argv[]) {
    if (argc!=5) {
        cout << "usage " << argv[0] << " d max_its example jmin" << endl; exit(1);
    }
    cout.precision(3);

    int d=atoi(argv[1]);
    int NumOfIterations=atoi(argv[2]);
    int example=atoi(argv[3]);
    int j0=atoi(argv[4]);

    SparseMWBasis1D     sparsemw_basis(d,j0);
    MW_MA               sparsemw_A(sparsemw_basis,1.);
    MW_Prec             sparsemw_prec(sparsemw_A);

    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,1.,0.,1.);

    Function<T>        rhs_func(refsol.rhs,refsol.sing_pts);
    RhsIntegral1D      rhsintegral1d_singular(sparsemw_basis, rhs_func, refsol.deltas, 20, true, false);
    RhsIntegral1D      rhsintegral1d_smooth(sparsemw_basis, rhs_func, refsol.deltas, 20, false, true);
    RhsIntegral1D      rhsintegral1d(sparsemw_basis, rhs_func, refsol.deltas, 20, true, true);

    Coefficients<Lexicographical,T,Index1D> f;
    f = initializeRHSVector(sparsemw_basis, rhsintegral1d_singular, rhsintegral1d_smooth, sparsemw_prec, refsol);
    Rhs F(rhsintegral1d,sparsemw_prec,f);

    Coefficients<Lexicographical,T,Index1D> u;

    SparseMW_GHS_ADWAV_SOLVER sparsemw_ghs_adwav_solver(sparsemw_A, F);
    u = sparsemw_ghs_adwav_solver.SOLVE(f.norm(2.), 1e-10, NumOfIterations, refsol.H1norm());

/*
    IndexSet<Index1D> Lambda = computeRHSLambda_SmoothPart(sparsemw_basis, -10., 10., sparsemw_basis.j0+2);
    cout << "Size of Lambda: " << Lambda.size() << endl;

    DenseVectorT rhs(Lambda.size()), x(Lambda.size());
    int count=1;
    for (const_set1d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
       rhs(count) = F(*it);
       ++count;
    }

    cout << "Assembling stiffness matrix started..." << endl;
    SparseMatrixT A(Lambda.size(),Lambda.size());
    sparsemw_A.toFlensSparseMatrix(Lambda,Lambda,A);
    DenseMatrixT A_dense;
    densify(cxxblas::NoTrans,A,A_dense);
    DenseMatrixT U(A.numRows(),A.numRows()), V(A.numCols(),A.numCols());
    int N = A_dense.numRows();
    DenseVectorT wr(N), wi(N);
    DenseMatrixT vl,vr;
    ev(false, false, A_dense, wr, wi, vl, vr);
    T cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
    for (int i=1; i<=wr.lastIndex(); ++i) {
        cB = std::min(cB,wr(i));
        CB = std::max(CB,wr(i));
    }
    cout << "... cB = " << cB << ", CB = " << CB << endl;
    cout << "... finished." << endl;

    cout << "CG solver started..." << endl;
    int iters = lawa::cg(A,x,rhs,1e-8);
    cout << "... finished after cg-iters: " << iters << endl;

    count=1;
    for (const_set1d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
       u[*it] = x(count);
       ++count;
    }



*/
    stringstream plot_filename;
    plot_filename << "ghs-adwav-realline-sparsemw-helmholtz1d-plot_" << example
                  << "_" << d << "_" << j0 << ".dat";
    cout << "Plot of solution started." << endl;
    plot<T, SparseMWBasis1D, MW_Prec>(sparsemw_basis, u, sparsemw_prec, refsol.u,
                                      refsol.d_u, -10., 10., pow2i<T>(-5),
                                      plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;

    return 0;
}


IndexSet<Index1D>
computeRHSLambda_SingularPart(const SparseMWBasis1D &basis, const DenseVectorT &_f_singularPoints,
                              int J_plus)
{
    IndexSet<Index1D> ret;

    T l1, l2;
    l1 = basis.mra.phi.max_support().l1, l2 = basis.mra.phi.max_support().l2;
    int numScaling =basis.mra.phi._numSplines;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        int k_left =  std::floor(float(pow2i<T>(basis.j0)*x-l2))-2;
        int k_right = std::ceil(float(pow2i<T>(basis.j0)*x-l1))+2;
        for (int k_help=k_left; k_help<=k_right; ++k_help) {
            for (int k=(k_help-1)*numScaling+1; k<=(k_help)*numScaling; ++k) {
                //cout << "Singular: Insert Scaling function -> " << basis.mra.phi.support(basis.j0,k) << endl;
                ret.insert(Index1D(basis.j0,k,XBSpline));
            }
        }
    }

    l1 = basis.psi.max_support().l1, l2 = basis.psi.max_support().l2;
    int numWavelets =basis.psi._numSplines;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        for (int j=basis.j0; j<=J_plus; ++j) {
            int k_left =  std::floor(float(pow2i<T>(j)*x-l2))-2;
            int k_right = std::ceil(float(pow2i<T>(j)*x-l1))+2;
            for (int k_help=k_left; k_help<=k_right; ++k_help) {
                for (int k=(k_help-1)*numWavelets+1; k<=(k_help)*numWavelets; ++k) {
                    //cout << "Singular: Insert Wavelet function -> " << basis.psi.support(j,k) << endl;
                    ret.insert(Index1D(j,k,XWavelet));
                }
            }
        }
    }

    return ret;
}

IndexSet<Index1D>
computeRHSLambda_SmoothPart(const SparseMWBasis1D &basis, T a, T b, int J_plus)
{
    IndexSet<Index1D> ret;

    cout << "computeRHSLambda_SmoothPart: j0 = " << basis.j0 << endl;
    T l1, l2;
    l1 = basis.mra.phi.max_support().l1, l2 = basis.mra.phi.max_support().l2;
    int numScaling =basis.mra.phi._numSplines;
    int k_left =  std::floor(float(pow2i<T>(basis.j0)*a-l2));
    int k_right = std::ceil(float(pow2i<T>(basis.j0)*b-l1));
    for (int k_help=k_left; k_help<=k_right; ++k_help) {
        for (int k=(k_help-1)*numScaling+1; k<=(k_help)*numScaling; ++k) {
            //cout << "Smooth: Insert Scaling function -> " << basis.mra.phi.support(basis.j0,k) << endl;
            ret.insert(Index1D(basis.j0,k,XBSpline));
        }
    }

    l1 = basis.psi.support(0,0).l1, l2 = basis.psi.support(0,0).l2;
    int numWavelets =basis.psi._numSplines;
    for (int j=basis.j0; j<=J_plus; ++j) {
        int k_left =  std::floor(float(pow2i<T>(j)*a-l2));
        int k_right = std::ceil(float(pow2i<T>(j)*b-l1));
        for (int k_help=k_left; k_help<=k_right; ++k_help) {
            for (int k=(k_help-1)*numWavelets+1; k<=(k_help)*numWavelets; ++k) {
                //cout << "Smooth: Insert Wavelet function -> " << basis.psi.support(j,k) << endl;
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }

    return ret;
}


template <typename T>
Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const SparseMWBasis1D &basis, const RhsIntegral1D &rhsintegral1d_singular,
                    const RhsIntegral1D &rhsintegral1d_smooth, MW_Prec &P,
                    RefSols_PDE_Realline1D<T> &refsol)
{

    T left_bound = -120., right_bound = 120.;
    int J_plus_smooth = 5, J_plus_singular = 30;
    bool singular_integral=true;
//    refsol.getRHS_W_XBSplineParameters(basis.d, basis.d_, left_bound, right_bound, J_plus_smooth,
//                                       J_plus_singular, singular_integral);


    Coefficients<Lexicographical,T,Index1D> f_singular;
    IndexSet<Index1D> LambdaRHS_singular = computeRHSLambda_SingularPart(basis,refsol.sing_pts,
                                                                         J_plus_singular);
    cout << "Initializing singular part of rhs, size of indexset: "
         << LambdaRHS_singular.size() << endl;
    for (const_set1d_it it=LambdaRHS_singular.begin(); it!=LambdaRHS_singular.end(); ++it) {
        f_singular[*it] = P(*it)*rhsintegral1d_singular(*it);
    }


    Coefficients<Lexicographical,T,Index1D> f_smooth;
    IndexSet<Index1D> LambdaRHS_smooth =computeRHSLambda_SmoothPart(basis,left_bound,
                                                                    right_bound,J_plus_smooth);
    cout << "Initializing smooth part of rhs, size of indexset: "
         << LambdaRHS_smooth.size() << endl;
    for (const_set1d_it it=LambdaRHS_smooth.begin(); it!=LambdaRHS_smooth.end(); ++it) {
        f_smooth[*it] = P(*it)*rhsintegral1d_smooth(*it);
    }
    if (singular_integral) {
        for (const_set1d_it it=LambdaRHS_singular.begin(); it!=LambdaRHS_singular.end(); ++it) {
            if (LambdaRHS_smooth.count(*it)>0) continue;
            f_smooth[*it] = P(*it)*rhsintegral1d_smooth(*it);
        }
    }

    return f_smooth + f_singular;
}

