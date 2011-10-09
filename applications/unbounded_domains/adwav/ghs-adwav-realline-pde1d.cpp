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
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;

//Operator definitions
typedef AdaptivePDEOperatorOptimized1D<T,Primal,R,SparseMulti>          SparseMW_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,SparseMulti>    SparseMW_PP_MA;

typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, SparseMW_MA>  SparseMW_Prec;

//Righthandsides definitions (tensor)
typedef RHSWithPeaks1D<T, SparseMW_Basis1D>                             SparseMW_RhsIntegral1D;

typedef RHS1D<T, SparseMW_RhsIntegral1D, SparseMW_Prec>                 SparseMW_Rhs;

typedef GHS_NONSYM_ADWAV<T, Index1D, SparseMW_MA, SparseMW_Rhs,
                         SparseMW_PP_MA, SparseMW_Rhs>               SparseMW_GHS_NONSYM_ADWAV_SOLVER;


int main (int argc, char *argv[]) {
    if (argc!=7) {
        cout << "usage " << argv[0] << " basistype d d_ jmin example max_its" << endl; exit(1);
    }
    cout.precision(3);

    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0=atoi(argv[4]);

    T reaction = 1., convection = 10., diffusion = 1.;
    int example=atoi(argv[5]);
    int NumOfIterations=atoi(argv[6]);

    stringstream rhsfilename;
    rhsfilename << "rhs/rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_"
                << argv[3] << "_" << argv[4] << "_" << 1. << "_" << argv[5] << ".dat";
    stringstream convfilename;
    convfilename << "ghs_adwav_conv_realline_pde1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << 1. << "_" << argv[5] << ".dat";
    stringstream plotfilename;
    plotfilename << "ghs_adwav_plot_realline_pde1d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << 1. << "_" << argv[5] << ".dat";

    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,reaction,convection,diffusion);
    Function<T>                 rhsFct(refsol.rhs,refsol.sing_pts);
    Function<T>                 PP_rhsFct(refsol.H1_rhs,refsol.sing_pts);

    SparseMW_Basis1D            SparseMW_basis(d,j0);
    SparseMW_MA                 SparseMW_A(SparseMW_basis,reaction,convection,diffusion);
    SparseMW_PP_MA              SparseMW_PP_A(SparseMW_basis,1.);
    SparseMW_Prec               SparseMW_prec(SparseMW_A);
    SparseMW_RhsIntegral1D      SparseMW_rhsintegral1d(SparseMW_basis, rhsFct, refsol.deltas, 17);
    SparseMW_RhsIntegral1D      SparseMW_PP_rhsintegral1d(SparseMW_basis, PP_rhsFct, refsol.H1_deltas, 17);
    SparseMW_Rhs                SparseMW_F(SparseMW_rhsintegral1d,SparseMW_prec);
    SparseMW_Rhs                SparseMW_PP_F(SparseMW_PP_rhsintegral1d,SparseMW_prec);
    SparseMW_GHS_NONSYM_ADWAV_SOLVER   SparseMW_ghs_adwav_solver(SparseMW_A, SparseMW_F,
                                                                 SparseMW_PP_A, SparseMW_PP_F,true);


    if (SparseMW_F.readIndexSets(rhsfilename.str().c_str()) ) {
        cout << "Index sets for rhs read... Ready to start."  << endl;
    }
    else {
        cout << "RHS: Could not open file." << endl;
        return 0;
    }

    Coefficients<Lexicographical,T,Index1D> w;
    w = SparseMW_ghs_adwav_solver.SOLVE(SparseMW_F.norm_estimate, 1e-10, convfilename.str().c_str(),
                                        NumOfIterations, refsol.H1norm());
    cout << "Plot of solution started." << endl;
    plot<T, SparseMW_Basis1D, SparseMW_Prec>(SparseMW_basis, w, SparseMW_prec, refsol.u,
                                             refsol.d_u, -80., 80., pow2i<T>(-5),
                                             plotfilename.str().c_str());
    cout << "Plot of solution finished." << endl;





    IndexSet<Index1D> Lambda;
    int max_k=15;
    for (int k=-max_k; k<=max_k; ++k) {
        Lambda.insert(Index1D(j0,k,XBSpline));
    }
    for (int j=j0; j<=j0+3; ++j) {
        for (int k=-pow2i<int>(j-j0)*max_k; k<=pow2i<int>(j-j0)*max_k; ++k) {
            Lambda.insert(Index1D(j,k,XWavelet));
        }
    }
    cout << "Size of Lambda: " << Lambda.size() << endl;

    Coefficients<Lexicographical,T,Index1D> f,u;
    f = SparseMW_F(Lambda);
    T res=0.;
    int numOfIterations = GMRES_Solve(Lambda, SparseMW_A, u, f, res, 1e-12);

    plot<T,SparseMW_Basis1D,SparseMW_Prec>(SparseMW_basis, u, SparseMW_prec, refsol.u, refsol.d_u, -80., 80.,
                                           pow2i<T>(-5), "convection2.dat");

    cout << "Number of iterations: " << numOfIterations << endl;

    T alpha, beta, gammaPrev, gamma, b_norm;
    Coefficients<Lexicographical,T,Index1D> b,x;
    Coefficients<Lexicographical,T,Index1D> r, q, s, p;
    b = f;
    x = f;
    b_norm = b.norm(2.);
    SparseMW_A.apply(x,0.,r,NoTrans);
    r -= b;
    r *= -1.;
    SparseMW_A.apply(r,0.,Lambda,s,Trans);
    p = s;
    gammaPrev = s*s;

    for (int k=1; k<=30; k++) {
        SparseMW_A.apply(p,0.,q,NoTrans);   //q = A*p;
        alpha = gammaPrev/(q*q);
        x +=   alpha *p;
        r -=   alpha*q;
        SparseMW_A.apply(r,0.,Lambda,s,Trans);  // flens::blas::mv(cxxblas::Trans, typename _cg<VB>::T(1), A, r, typename _cg<VB>::T(0), s);

        gamma = s*s;
        cout << "Iteration " << k << ": gamma = " << gamma << endl;
        /*
        if (sqrt(gamma)<=1e-12) {
            std::cerr << "    cgls: gamma = " << gamma << std::endl;
            return k-1;
        }
        */
        beta  = gamma/gammaPrev;
        p *= beta;
        p += s;
        gammaPrev = gamma;
    }

    plot<T,SparseMW_Basis1D,SparseMW_Prec>(SparseMW_basis, x, SparseMW_prec, refsol.u, refsol.d_u, -80., 80.,
                                           pow2i<T>(-5), "convection2.dat");

    return 0;
}
