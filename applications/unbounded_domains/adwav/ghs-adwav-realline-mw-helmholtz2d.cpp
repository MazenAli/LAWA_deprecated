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
typedef IndexSet<Index2D>::const_iterator                               const_set2d_it;
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


IndexSet<Index1D>
computeRHSLambda_SmoothPart(const MWBasis1D &basis, T a, T b, int J_plus);

template<typename T, typename Basis2D>
Coefficients<Lexicographical,T,Index2D>
precomputeRHS(int example, const Basis2D &mwbasis2d, MW_Prec &mw_prec, int j0_x, int j0_y,
               int J, T r_x, T r_y);

int main (int argc, char *argv[]) {
    if (argc!=4) {
        cout << "usage " << argv[0] << " d max_its example [jmin_x jmin_y]" << endl; exit(1);
    }
    cout.precision(16);

    int d=atoi(argv[1]);
    int NumOfIterations=atoi(argv[2]);
    int example=atoi(argv[3]);
    int j0_x=-1;
    int j0_y=-2;
    int order=35;
    T r=5.;

    MWBasis1D mw_basis_x(d,j0_x);
    MWBasis1D mw_basis_y(d,j0_y);

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

    Coefficients<Lexicographical,T,Index2D> f, u;
    f = precomputeRHS(example, mw_basis2d, mw_prec, j0_x, j0_y, 7, 15., 15.);

    for (const_coeff2d_it it=f.begin(); it!=f.end(); ++it) {
        if (((*it).first.index2.xtype==XBSpline) && ((*it).first.index2.j!=j0_y)) {
            cout << "Error in RHS: " << (*it).first << endl;
        }
    }

    Rhs mw_F_test(rhsintegral2d,mw_prec);
    Rhs mw_F(rhsintegral2d,mw_prec,f);

    MW_GHS_ADWAV_SOLVER mw_ghs_adwav_solver(mw_A,mw_F);
    u = mw_ghs_adwav_solver.SOLVE(f.norm(2.), 1e-5, NumOfIterations, refsol.H1norm());

/*
    IndexSet<Index1D> Lambda_x, Lambda_y;
    Lambda_x = computeRHSLambda_SmoothPart(mw_basis_x, -r, r, 1);
    Lambda_y = computeRHSLambda_SmoothPart(mw_basis_y, -r, r, 1);

    IndexSet<Index2D> Lambda;
    /*
    for (const_set1d_it it_x=Lambda_x.begin(); it_x!=Lambda_x.end(); ++it_x) {
        for (const_set1d_it it_y=Lambda_y.begin(); it_y!=Lambda_y.end(); ++it_y) {
            Lambda.insert(Index2D(*it_x,*it_y));
        }
    }

    cout << Lambda << endl;
    cout << "Size of Lambda: " << Lambda.size() << endl;

    DenseVectorT rhs(Lambda.size()), x(Lambda.size());
    int count=1;
    for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        rhs(count) = mw_F_test(*it);
        ++count;
    }

    SparseMatrixT A(Lambda.size(),Lambda.size());
    mw_A.toFlensSparseMatrix(Lambda,Lambda,A);

    int iters = lawa::cg(A,x,rhs,1e-8);
    cout << "cg-iters: " << iters << endl;


    count=1;
    for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        u[*it] = x(count);
        ++count;
    }
*/




    stringstream plot_filename;
    plot_filename << "ghs-adwav-helmholtz2d-realline-mw-plot.dat";
    ofstream plotfile(plot_filename.str().c_str());
    cout << "Plot of solution started." << endl;
    T a1=-5., b1=5, a2=-5., b2=5.;
    T h1=0.125, h2=0.125;
    for (T x=a1; x<=b1; x+=h1) {
        for (T y=a2; y<=b2; y+=h2) {
            T appr = 0.0;
            T exact= refsol.exact(x,y);
            for (const_coeff2d_it it = u.begin(); it != u.end(); ++it) {
                XType xtype_x = (*it).first.index1.xtype;
                XType xtype_y = (*it).first.index2.xtype;
                int j_x = (*it).first.index1.j, k_x = (*it).first.index1.k;
                int j_y = (*it).first.index2.j, k_y = (*it).first.index2.k;

                T coeff = (*it).second, prec = mw_A.prec((*it).first);

                appr    += prec * coeff * mw_basis2d.first.generator(xtype_x)(x,j_x,k_x,0)
                                        * mw_basis2d.second.generator(xtype_y)(y,j_y,k_y,0);
            }
            plotfile << x << " " << y << " " << exact << " " << appr  << " " << endl;
        }
        plotfile << std::endl;
    }
    cout << "Plot of solution finished." << endl;
    return 0;
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

Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const RHS1D &rhsintegral,
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

template<typename T, typename Basis2D>
Coefficients<Lexicographical,T,Index2D>
precomputeRHS(int example, const Basis2D &mw_basis2d, MW_Prec &mw_prec, int j0_x, int j0_y,
               int J, T r_x, T r_y)
{
    // Assemble f vector and u vector
    DenseMatrixT nodeltas;
    TensorRefSols_PDE_Realline2D<T> refsol;
    refsol.setExample(example, 1.);
    Function<T> u1Fct(refsol.exact_x, refsol.sing_pts_x);
    Function<T> u2Fct(refsol.exact_y, refsol.sing_pts_y);
    RHS1D       u1_integral(mw_basis2d.first, u1Fct, nodeltas, 35);
    RHS1D       u2_integral(mw_basis2d.second, u2Fct, nodeltas, 35);

    Function<T> f1Fct(refsol.rhs_x, refsol.sing_pts_x);
    Function<T> f2Fct(refsol.rhs_y, refsol.sing_pts_y);
    RHS1D       f1_integral(mw_basis2d.first, f1Fct, refsol.deltas_x, 35);
    RHS1D       f2_integral(mw_basis2d.second, f2Fct, refsol.deltas_y, 35);

    IndexSet<Index1D> Lambda_smooth, Lambda_singular;
    Coefficients<Lexicographical,T,Index1D> u1_coeff, u2_coeff, f1_coeff, f2_coeff;

    Lambda_smooth = computeRHSLambda_SmoothPart(mw_basis2d.first, -r_x, r_x, J);
    Lambda_singular = computeRHSLambda_SingularPart(mw_basis2d.first,refsol.sing_pts_x ,30);
    u1_coeff = initializeRHSVector(u1_integral, Lambda_singular, Lambda_smooth );

    Lambda_smooth = computeRHSLambda_SmoothPart(mw_basis2d.second, -r_y, r_y, J);
    Lambda_singular = computeRHSLambda_SingularPart(mw_basis2d.second,refsol.sing_pts_y ,30);
    u2_coeff = initializeRHSVector(u2_integral, Lambda_singular, Lambda_smooth );

    Lambda_smooth = computeRHSLambda_SmoothPart(mw_basis2d.first, -r_x, r_x, J);
    Lambda_singular = computeRHSLambda_SingularPart(mw_basis2d.first,refsol.sing_pts_x ,30);
    f1_coeff = initializeRHSVector(f1_integral, Lambda_singular, Lambda_smooth );

    Lambda_smooth = computeRHSLambda_SmoothPart(mw_basis2d.second, -r_y, r_y, J);
    Lambda_singular = computeRHSLambda_SingularPart(mw_basis2d.second,refsol.sing_pts_y ,30);
    f2_coeff = initializeRHSVector(f2_integral, Lambda_singular, Lambda_smooth );



    cout << "Computation of 1d vectors for u and f finished." << endl;

    u1_coeff = ABSOLUTE_THRESH(u1_coeff,1e-10);
    u2_coeff = ABSOLUTE_THRESH(u2_coeff,1e-10);
    f1_coeff = ABSOLUTE_THRESH(f1_coeff,1e-10);
    f2_coeff = ABSOLUTE_THRESH(f2_coeff,1e-10);

    cout << "Sizes of 1d u vectors after thresholding: " << u1_coeff.size() << ", " << u2_coeff.size() << endl;
    cout << "Sizes of 1d f vectors after thresholding: " << f1_coeff.size() << ", " << f2_coeff.size() << endl;

    cout << "Computation of 2d vector u and f started." << endl;

    Coefficients<Lexicographical,T,Index2D> f_coeff;


    for (const_coeff1d_it it_x=u1_coeff.begin(); it_x!=u1_coeff.end(); ++it_x) {
       for (const_coeff1d_it it_y=f2_coeff.begin(); it_y!=f2_coeff.end(); ++it_y) {
           Index2D index((*it_x).first,(*it_y).first);
           T f_val = (*it_x).second * (*it_y).second * mw_prec(index);
           if (fabs(f_val)>1e-10) {
               f_coeff[index] += f_val;
           }
       }
    }
    for (const_coeff1d_it it_x=f1_coeff.begin(); it_x!=f1_coeff.end(); ++it_x) {
       for (const_coeff1d_it it_y=u2_coeff.begin(); it_y!=u2_coeff.end(); ++it_y) {
           Index2D index((*it_x).first,(*it_y).first);
           T f_val = (*it_x).second * (*it_y).second * mw_prec(index);
           if (fabs(f_val)>1e-10) {
               f_coeff[index] += f_val;
           }
       }
    }

    cout << "Size of 2d vector f: " << f_coeff.size() << endl;

    return f_coeff;
}
