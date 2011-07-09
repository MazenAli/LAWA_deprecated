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
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator         const_coeff2d_it;
typedef Coefficients<AbsoluteValue,T,Index1D>::const_iterator           const_coeff1d_abs_it;
typedef Coefficients<AbsoluteValue,T,Index2D>::const_iterator           const_coeff2d_abs_it;

//Basis definitions
typedef Basis<T,Orthogonal,R,Multi>                                     MWBasis1D;
typedef TensorBasis2D<Adaptive, MWBasis1D, MWBasis1D>                   MWBasis2D;

//Operator definitions
typedef AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>                     MW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D, MW_MA>        MW_Prec;

//Righthandsides definitions (tensor)
typedef RHSWithPeaks1D<T, MWBasis1D>                                    RHS1D;
typedef SeparableRHS2D<T,MWBasis2D>                                     SeparableRhsIntegral2D;
typedef SumOfTwoRHSIntegrals<T,Index2D,SeparableRhsIntegral2D,
                             SeparableRhsIntegral2D>                    SumOfSeparableRhsIntegral2D;
typedef RHS<T,Index2D, SumOfSeparableRhsIntegral2D,
            MW_Prec>                                                    Rhs;

typedef GHS_ADWAV<T, Index2D, MW_MA, Rhs>                               MW_GHS_ADWAV_SOLVER;

template<typename T>
IndexSet<Index1D>
ReferenceLambda(int jmin, int jmax, const MWBasis1D &basis, T radius);

template<typename T, typename Basis2D>
Coefficients<Lexicographical,T,Index2D>
precoomputeRHS(int example, const Basis2D &mwbasis2d, MW_Prec &mw_prec, int j0, int J, int r, T thresh);

int main (int argc, char *argv[]) {
    if (argc!=4) {
        cout << "usage " << argv[0] << " d max_its example [jmin_x jmin_y]" << endl; exit(1);
    }
    cout.precision(16);

    int d=atoi(argv[1]);
    int NumOfIterations=atoi(argv[2]);
    int example=atoi(argv[3]);
    int j0=-2;
    int jmin_x=j0, jmin_y=j0;
    int order=20;


    MWBasis1D mw_basis_x(d,jmin_x);
    MWBasis1D mw_basis_y(d,jmin_y);

    MWBasis2D           mw_basis2d(mw_basis_x,mw_basis_y);
    MW_MA               mw_A(mw_basis2d,1.);
    MW_Prec             mw_prec(mw_A);

    TensorRefSols_PDE_Realline2D<T> refsol;
    refsol.setExample(example, 1.);
    SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                    refsol.exact_y, refsol.sing_pts_y);

    SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                    refsol.rhs_y, refsol.sing_pts_y);
    GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;
    SeparableRhsIntegral2D rhsintegral_x(mw_basis2d, SepFunc1, refsol.deltas_x, no_deltas, order);
    SeparableRhsIntegral2D rhsintegral_y(mw_basis2d, SepFunc2, no_deltas, refsol.deltas_y, order);
    SumOfSeparableRhsIntegral2D rhsintegral2d(rhsintegral_x,rhsintegral_y);

    Coefficients<Lexicographical,T,Index2D> f;
    f = precoomputeRHS(example, mw_basis2d, mw_prec, j0, 8, 13., 1e-8);

    Rhs mw_F_test(rhsintegral2d,mw_prec);
    Rhs mw_F(rhsintegral2d,mw_prec,f);


    Index2D index(Index1D(0,0,XWavelet),Index1D(0,1,XWavelet));
    cout << mw_F_test(index) << " " << mw_F(index) << endl;

    MW_GHS_ADWAV_SOLVER mw_ghs_adwav_solver(mw_A,mw_F);

    mw_ghs_adwav_solver.SOLVE(f.norm(2.), 1e-5, NumOfIterations, refsol.H1norm());

    return 0;
}

template<typename T>
IndexSet<Index1D>
ReferenceLambda(int jmin, int jmax, const MWBasis1D &basis, T radius)
{
    IndexSet<Index1D> Lambda;
    int k_left, k_right;
    int numWavelets = (int)basis.psi._numSplines;
    int numScaling = (int)basis.mra.phi._numSplines;

    int wavelet_count = 0;

    for (int j=jmin; j<=jmax; ++j) {
        k_left = std::floor(-pow2i<T>(j)*radius-basis.psi.max_support().l2);
        k_right = std::ceil(pow2i<T>(j)*radius-basis.psi.max_support().l1);
        for (int k_help=k_left; k_help<=k_right; ++k_help) {
            for (int k=(k_help-1)*numWavelets+1; k<=k_help*numWavelets; ++k) {
                Lambda.insert(Index1D(j,k,XWavelet));
                ++wavelet_count;
            }
        }
    }

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

    return Lambda;
}

template<typename T, typename Basis2D>
Coefficients<Lexicographical,T,Index2D>
precoomputeRHS(int example, const Basis2D &mw_basis2d, MW_Prec &mw_prec, int j0, int J, int r, T thresh)
{
    int order = 20;
    DenseMatrixT nodeltas;
    TensorRefSols_PDE_Realline2D<T> refsol;
    refsol.setExample(example, 1.);
    Function<T> u1Fct(refsol.exact_x, refsol.sing_pts_x);
    Function<T> u2Fct(refsol.exact_y, refsol.sing_pts_y);
    RHS1D       u1_integral(mw_basis2d.first, u1Fct, nodeltas, order);
    RHS1D       u2_integral(mw_basis2d.second, u2Fct, nodeltas, order);
    Function<T> f1Fct(refsol.rhs_x, refsol.sing_pts_x);
    Function<T> f2Fct(refsol.rhs_y, refsol.sing_pts_y);
    RHS1D       f1_integral(mw_basis2d.first, f1Fct, nodeltas, order);
    RHS1D       f2_integral(mw_basis2d.second, f2Fct, nodeltas, order);

    IndexSet<Index1D> Lambda1d = ReferenceLambda(j0,J,mw_basis2d.first,r);
    cout << "Size of reference index set: " << Lambda1d.size() << endl;
    cout << "Computation of 1d vectors for u and f started." << endl;
    Coefficients<Lexicographical,T,Index1D> u1_coeff, u2_coeff, f1_coeff, f2_coeff;
    for (const_set1d_it it=Lambda1d.begin(); it!=Lambda1d.end(); ++it) {
        u1_coeff[*it] = u1_integral(*it);
        u2_coeff[*it] = u2_integral(*it);
        f1_coeff[*it] = f1_integral(*it);
        f2_coeff[*it] = f2_integral(*it);
    }
    cout << "Computation of 1d vectors for u and f finished." << endl;
    cout << "Computation of 2d vector u and f started." << endl;

    Coefficients<AbsoluteValue,T,Index1D> u1_coeff_abs, u2_coeff_abs, f1_coeff_abs, f2_coeff_abs;
    u1_coeff_abs = u1_coeff;
    u2_coeff_abs = u2_coeff;
    u1_coeff = THRESH(u1_coeff_abs,1e-8);
    u2_coeff = THRESH(u2_coeff_abs,1e-10);
    f1_coeff_abs = f1_coeff;
    f2_coeff_abs = f2_coeff;
    f1_coeff = THRESH(f1_coeff_abs,1e-10);
    f2_coeff = THRESH(f2_coeff_abs,1e-10);

    cout << "Sizes of 1d u vectors after thresholding: " << u1_coeff.size() << ", " << u2_coeff.size() << endl;
    cout << "Sizes of 1d f vectors after thresholding: " << f1_coeff.size() << ", " << f2_coeff.size() << endl;

    Coefficients<Lexicographical,T,Index2D> f_coeff;

    for (const_coeff1d_it it_x=u1_coeff.begin(); it_x!=u1_coeff.end(); ++it_x) {
       for (const_coeff1d_it it_y=f2_coeff.begin(); it_y!=f2_coeff.end(); ++it_y) {
           Index2D index((*it_x).first,(*it_y).first);
           T f_val = (*it_x).second * (*it_y).second * mw_prec(index);
           if (fabs(f_val)>1e-12) {
               f_coeff[index] += f_val;
           }
       }
    }
    for (const_coeff1d_it it_x=f1_coeff.begin(); it_x!=f1_coeff.end(); ++it_x) {
       for (const_coeff1d_it it_y=u2_coeff.begin(); it_y!=u2_coeff.end(); ++it_y) {
           Index2D index((*it_x).first,(*it_y).first);
           T f_val = (*it_x).second * (*it_y).second * mw_prec(index);
           if (fabs(f_val)>1e-12) {
               f_coeff[index] += f_val;
           }
       }
    }
    cout << "Size of 2d vector f: " << f_coeff.size() << endl;

    return f_coeff;
}
