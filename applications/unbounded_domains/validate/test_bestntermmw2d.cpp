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
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,R,Multi,
                                               Orthogonal,R,Multi>      MW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D, MW_MA>        MW_Prec;

typedef RHSWithPeaks1D<T, MWBasis1D>                                    RhsIntegral1D;


IndexSet<Index1D>
computeRHSLambda_SingularPart(const MWBasis1D &basis, const DenseVectorT &_f_singularPoints,
                              int J_plus);

IndexSet<Index1D>
computeRHSLambda_SmoothPart(const MWBasis1D &basis, T a, T b, int J_plus);

Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const RhsIntegral1D &rhsintegral,
                    const IndexSet<Index1D> LambdaRHS_singular,
                    const IndexSet<Index1D> LambdaRHS_smooth);

int main(int argc, char *argv[]) {
    int d=2;
    int j0_x=-2;
    int j0_y=0;
    int J=atoi(argv[1]);
    int example = 1;
    T   r_x = 10.;
    T   r_y = 10.;
    T   c = 1.;

    MWBasis1D           mwbasis1d_x(d,j0_x);
    MWBasis1D           mwbasis1d_y(d,j0_y);
    MWBasis2D           mwbasis2d(mwbasis1d_x,mwbasis1d_y);

    MW_MA               mw_A(mwbasis2d,1.);
    MW_Prec             mw_prec(mw_A);

    // Assemble f vector and u vector
    DenseMatrixT nodeltas;
    TensorRefSols_PDE_Realline2D<T> refsol;
    refsol.setExample(example, c);
    Function<T> u1Fct(refsol.exact_x, refsol.sing_pts_x);
    Function<T> u2Fct(refsol.exact_y, refsol.sing_pts_y);
    RhsIntegral1D       u1_integral(mwbasis2d.first, u1Fct, nodeltas, 35);
    RhsIntegral1D       u2_integral(mwbasis2d.second, u2Fct, nodeltas, 35);

    Function<T> f1Fct(refsol.rhs_x, refsol.sing_pts_x);
    Function<T> f2Fct(refsol.rhs_y, refsol.sing_pts_y);
    RhsIntegral1D       f1_integral(mwbasis2d.first, f1Fct, refsol.deltas_x, 35);
    RhsIntegral1D       f2_integral(mwbasis2d.second, f2Fct, refsol.deltas_y, 35);

    IndexSet<Index1D> Lambda_smooth, Lambda_singular;
    Coefficients<Lexicographical,T,Index1D> u1_coeff, u2_coeff, f1_coeff, f2_coeff;

    Lambda_smooth = computeRHSLambda_SmoothPart(mwbasis1d_x, -r_x, r_x, J);
    Lambda_singular = computeRHSLambda_SingularPart(mwbasis1d_x,refsol.sing_pts_x ,30);
    u1_coeff = initializeRHSVector(u1_integral, Lambda_singular, Lambda_smooth );

    Lambda_smooth = computeRHSLambda_SmoothPart(mwbasis1d_y, -r_y, r_y, J);
    Lambda_singular = computeRHSLambda_SingularPart(mwbasis1d_y,refsol.sing_pts_y ,30);
    u2_coeff = initializeRHSVector(u2_integral, Lambda_singular, Lambda_smooth );

    Lambda_smooth = computeRHSLambda_SmoothPart(mwbasis1d_x, -r_x, r_x, J);
    Lambda_singular = computeRHSLambda_SingularPart(mwbasis1d_x,refsol.sing_pts_x ,30);
    f1_coeff = initializeRHSVector(f1_integral, Lambda_singular, Lambda_smooth );

    Lambda_smooth = computeRHSLambda_SmoothPart(mwbasis1d_y, -r_y, r_y, J);
    Lambda_singular = computeRHSLambda_SingularPart(mwbasis1d_y,refsol.sing_pts_y ,30);
    f2_coeff = initializeRHSVector(f2_integral, Lambda_singular, Lambda_smooth );



    cout << "Computation of 1d vectors for u and f finished." << endl;

    u1_coeff = ABSOLUTE_THRESH(u1_coeff,1e-10);
    u2_coeff = ABSOLUTE_THRESH(u2_coeff,1e-10);
    f1_coeff = ABSOLUTE_THRESH(f1_coeff,1e-10);
    f2_coeff = ABSOLUTE_THRESH(f2_coeff,1e-10);

    cout << "Sizes of 1d u vectors after thresholding: " << u1_coeff.size() << ", " << u2_coeff.size() << endl;
    cout << "Sizes of 1d f vectors after thresholding: " << f1_coeff.size() << ", " << f2_coeff.size() << endl;

    cout << "Computation of 2d vector u and f started." << endl;

    Coefficients<Lexicographical,T,Index2D> u_coeff, f_coeff;

    for (const_coeff1d_it it_x=u1_coeff.begin(); it_x!=u1_coeff.end(); ++it_x) {
        for (const_coeff1d_it it_y=u2_coeff.begin(); it_y!=u2_coeff.end(); ++it_y) {
            Index2D index((*it_x).first,(*it_y).first);
            T u_val = (*it_x).second * (*it_y).second * (1./mw_prec(index));
            if (fabs(u_val)>1e-6) {
                u_coeff[index] = u_val;
                T f_val = (  f1_coeff[(*it_x).first] * (*it_y).second
                           + f2_coeff[(*it_y).first] * (*it_x).second ) * mw_prec(index);
                f_coeff[index] = f_val;
            }
        }
    }
    cout << "Computation of 2d vectors u and f finished." << endl;
    cout << "Size of 2d u vector: " << u_coeff.size() << endl;
    cout << "Size of 2d f vector: " << f_coeff.size() << endl;
    //u_coeff=THRESH(u_coeff,1e-10);



    Coefficients<AbsoluteValue,T,Index2D> u_coeff_abs;
    u_coeff_abs = u_coeff;


    for(int n=1; n<=21000; n+=1000) {
        int count=1;
        Coefficients<Lexicographical,T,Index2D> u, Au;
        for (const_coeff2d_abs_it it=u_coeff_abs.begin(); it!=u_coeff_abs.end(); ++it,++count) {
            u[(*it).second] = (*it).first;
            if (count>n) break;
        }
        T fu = u*f_coeff;
        Au = mw_A.mv(supp(u), u);
        T uAu = u*Au;
        cout << n << " " << sqrt(fabs(std::pow(refsol.H1norm(),2.)- 2*fu + uAu)) << endl;
    }



    return 0.;
}

IndexSet<Index1D>
computeRHSLambda_SingularPart(const MWBasis1D &basis, const DenseVectorT &_f_singularPoints,
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
computeRHSLambda_SmoothPart(const MWBasis1D &basis, T a, T b, int J_plus)
{
    IndexSet<Index1D> ret;

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

Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const RhsIntegral1D &rhsintegral,
                    const IndexSet<Index1D> LambdaRHS_singular,
                    const IndexSet<Index1D> LambdaRHS_smooth)
{

    Coefficients<Lexicographical,T,Index1D> f;

    for (const_set1d_it it=LambdaRHS_singular.begin(); it!=LambdaRHS_singular.end(); ++it) {
        T tmp = rhsintegral(*it);
        if (fabs(tmp)>0) f[*it] = tmp;
    }

    for (const_set1d_it it=LambdaRHS_smooth.begin(); it!=LambdaRHS_smooth.end(); ++it) {
        if (LambdaRHS_singular.count(*it)>0) continue;
        f[*it] = rhsintegral(*it);
    }

    return f;

}
