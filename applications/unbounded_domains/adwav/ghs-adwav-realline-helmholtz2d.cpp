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
typedef Basis<T,Orthogonal,R,Multi>                                     MW_Basis1D;
typedef TensorBasis2D<Adaptive, MW_Basis1D, MW_Basis1D>                 MW_Basis2D;
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;
typedef TensorBasis2D<Adaptive, SparseMW_Basis1D,SparseMW_Basis1D>      SparseMW_Basis2D;

//Operator definitions
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,R,Multi,
                                               Orthogonal,R,Multi>      MW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D, MW_MA>        MW_Prec;
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Primal,R,SparseMulti,
                                               Primal,R,SparseMulti>    SparseMW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D,SparseMW_MA>   SparseMW_Prec;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,MW_Basis2D>                                    MW_SeparableRhsIntegral2D;
typedef SeparableRHS2D<T,SparseMW_Basis2D>                              SparseMW_SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,MW_SeparableRhsIntegral2D,
                             MW_SeparableRhsIntegral2D>                 MW_SumOfSeparableRhsIntegral2D;
typedef SumOfTwoRHSIntegrals<T,Index2D,SparseMW_SeparableRhsIntegral2D,
                             SparseMW_SeparableRhsIntegral2D>           SparseMW_SumOfSeparableRhsIntegral2D;

typedef RHS2D<T,MW_SumOfSeparableRhsIntegral2D,MW_Prec>                 MW_SumOfSeparableRhs;
typedef RHS2D<T,SparseMW_SumOfSeparableRhsIntegral2D,SparseMW_Prec>     SparseMW_SumOfSeparableRhs;

//Righthandsides definitions (non-separable)
typedef SmoothRHSWithAlignedSing2D<T, MW_Basis2D, SparseGridGP>         MW_NonSeparableRhsIntegralSG2D;
typedef SmoothRHSWithAlignedSing2D<T, SparseMW_Basis2D, SparseGridGP>   SparseMW_NonSeparableRhsIntegralSG2D;

typedef SmoothRHSWithAlignedSing2D<T, MW_Basis2D, FullGridGL>           MW_NonSeparableRhsIntegralFG2D;
typedef SmoothRHSWithAlignedSing2D<T, SparseMW_Basis2D, FullGridGL>     SparseMW_NonSeparableRhsIntegralFG2D;

typedef SumOfThreeRHSIntegrals<T, Index2D,
                               MW_NonSeparableRhsIntegralFG2D>          MW_SumOfNonSeparableRhsIntegral2D;
typedef SumOfThreeRHSIntegrals<T, Index2D,
                               SparseMW_NonSeparableRhsIntegralFG2D>    SparseMW_SumOfNonSeparableRhsIntegral2D;

typedef RHS2D<T,MW_NonSeparableRhsIntegralSG2D, MW_Prec>                MW_NonSeparableRhs2D;
typedef RHS2D<T,SparseMW_NonSeparableRhsIntegralSG2D, SparseMW_Prec>    SparseMW_NonSeparableRhs2D;

typedef RHS2D<T,MW_SumOfNonSeparableRhsIntegral2D,MW_Prec>              MW_SumOfNonSeparableRhs2D;
typedef RHS2D<T,SparseMW_SumOfNonSeparableRhsIntegral2D,SparseMW_Prec>  SparseMW_SumOfNonSeparableRhs2D;


//Algorithm definition
typedef GHS_ADWAV<T, Index2D, MW_MA, MW_SumOfSeparableRhs>              MW_GHS_ADWAV_SOLVER_SumofSeparableRhs;
typedef GHS_ADWAV<T, Index2D, SparseMW_MA, SparseMW_SumOfSeparableRhs>  SparseMW_GHS_ADWAV_SOLVER_SumofSeparableRhs;

typedef GHS_ADWAV<T, Index2D, MW_MA, MW_NonSeparableRhs2D>              MW_GHS_ADWAV_SOLVER_NonSeparableRhs;
typedef GHS_ADWAV<T, Index2D, SparseMW_MA, SparseMW_NonSeparableRhs2D>  SparseMW_GHS_ADWAV_SOLVER_NonSeparableRhs;

typedef GHS_ADWAV<T, Index2D, MW_MA, MW_SumOfNonSeparableRhs2D>         MW_GHS_ADWAV_SOLVER_SumNonSeparableRhs;
typedef GHS_ADWAV<T, Index2D, SparseMW_MA,
                  SparseMW_SumOfNonSeparableRhs2D>                      SparseMW_GHS_ADWAV_SOLVER_SumNonSeparableRhs;


IndexSet<Index1D>
computeRHSLambda_SmoothPart(const MW_Basis1D &basis, T a, T b, int J_plus);

template<typename T, typename Basis2D>
Coefficients<Lexicographical,T,Index2D>
precomputeRHS(int example, const Basis2D &mwbasis2d, MW_Prec &mw_prec, int j0_x, int j0_y,
               int J, T r_x, T r_y);

int main (int argc, char *argv[]) {
    if (argc!=8) {
        cout << "usage " << argv[0] << " basistype d d_ j0_x j0_y example NumOfIterations" << endl; exit(1);
    }
    cout.precision(20);

    int d   =atoi(argv[2]);
    int d_  =atoi(argv[3]);
    int j0_x=atoi(argv[4]);
    int j0_y=atoi(argv[5]);
    T c = 1.;
    int example=atoi(argv[6]);
    int NumOfIterations=atoi(argv[7]);

    stringstream rhsfilename;
    rhsfilename << "rhs/rhs_realline_pde2d_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                << argv[4] << "_" << argv[5] << "_cx_0_cy_0_" << argv[6] << ".dat";

    /*
    stringstream rhsfilename;
    rhsfilename << "rhs/rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";
    */
    stringstream convfilename;
    convfilename << "ghs_adwav_conv_realline_helmholtz2d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";

    int order=20;

    if (strcmp(argv[1],"MW")==0) {
        MW_Basis1D MW_basis_x(d,j0_x);
        MW_Basis1D MW_basis_y(d,j0_y);
        MW_Basis2D MW_basis2d(MW_basis_x,MW_basis_y);
        MW_MA      MW_A(MW_basis2d,1.);
        MW_Prec    MW_P(MW_A);
        int assemble_matrix = 2;    //assemble matrix, use operator internal method

        if (example==1 || example==2 || example==3) {

            TensorRefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(example, 1.,0., 0., 1.);
            SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                            refsol.exact_y, refsol.sing_pts_y);

            SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                            refsol.rhs_y, refsol.sing_pts_y);
            GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;
            MW_SeparableRhsIntegral2D MW_rhsintegral_x(MW_basis2d, SepFunc1, refsol.deltas_x, no_deltas, order);
            MW_SeparableRhsIntegral2D MW_rhsintegral_y(MW_basis2d, SepFunc2, no_deltas, refsol.deltas_y, order);
            MW_SumOfSeparableRhsIntegral2D MW_rhsintegral2d(MW_rhsintegral_x,MW_rhsintegral_y);
            MW_SumOfSeparableRhs MW_F(MW_rhsintegral2d,MW_P);
            MW_F.readIndexSets(rhsfilename.str().c_str());

            MW_GHS_ADWAV_SOLVER_SumofSeparableRhs MW_ghs_adwav_solver(MW_A,MW_F,true,assemble_matrix);


            Coefficients<Lexicographical,T,Index2D> u;
            u = MW_ghs_adwav_solver.SOLVE(MW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                          NumOfIterations, refsol.H1norm());
        }
        else if (example==4) {
            int order = 31;
            RefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(1, 1.);
            Function2D<T> Func2d(refsol.rhs, refsol.sing_pts_x, refsol.sing_pts_y);
            MW_NonSeparableRhsIntegralSG2D MW_rhsintegral2d(MW_basis2d, Func2d, order);
            MW_NonSeparableRhs2D           MW_F(MW_rhsintegral2d,MW_P);
            MW_F.readIndexSets(rhsfilename.str().c_str());

            MW_GHS_ADWAV_SOLVER_NonSeparableRhs MW_ghs_adwav_solver(MW_A,MW_F,true,assemble_matrix);

            Coefficients<Lexicographical,T,Index2D> u;
            u = MW_ghs_adwav_solver.SOLVE(MW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                          NumOfIterations, refsol.H1norm());
        }
        else if (example==5) {
            RefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(2, 1.);
            Function2D<T> Func2d(refsol.exact, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_x(refsol.exact_dx, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_y(refsol.exact_dy, refsol.sing_pts_x, refsol.sing_pts_y);

            if (d==2) {
                int order = 20;
                MW_NonSeparableRhsIntegralFG2D    MW_rhsintegral_reaction(MW_basis2d, Func2d, order);
                MW_NonSeparableRhsIntegralFG2D    MW_rhsintegral_diffusion_x(MW_basis2d, Func2d_x, order, 1, 0);
                MW_NonSeparableRhsIntegralFG2D    MW_rhsintegral_diffusion_y(MW_basis2d, Func2d_y, order, 0, 1);
                MW_SumOfNonSeparableRhsIntegral2D MW_rhsintegral2d(MW_rhsintegral_diffusion_x,
                                                                   MW_rhsintegral_diffusion_y,
                                                                   MW_rhsintegral_reaction);
                MW_SumOfNonSeparableRhs2D         MW_F(MW_rhsintegral2d,MW_P);
                MW_F.readIndexSets(rhsfilename.str().c_str());

                MW_GHS_ADWAV_SOLVER_SumNonSeparableRhs MW_ghs_adwav_solver(MW_A,MW_F,true,assemble_matrix);

                Coefficients<Lexicographical,T,Index2D> u;
                u = MW_ghs_adwav_solver.SOLVE(MW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                              NumOfIterations, refsol.H1norm());
            }
        }
    }
    else if (strcmp(argv[1],"SparseMW")==0) {
        SparseMW_Basis1D SparseMW_basis_x(d,j0_x);
        SparseMW_Basis1D SparseMW_basis_y(d,j0_y);
        SparseMW_Basis2D SparseMW_basis2d(SparseMW_basis_x,SparseMW_basis_y);
        SparseMW_MA      SparseMW_A(SparseMW_basis2d,1.);
        SparseMW_Prec    SparseMW_P(SparseMW_A);
        int assemble_matrix = 0;    //compute matrix vector product solely on index set
        if (example==1 || example==2 || example==3) {

            TensorRefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(example, 1., 0., 0., 1.);
            SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                            refsol.exact_y, refsol.sing_pts_y);

            SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                            refsol.rhs_y, refsol.sing_pts_y);
            GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;
            SparseMW_SeparableRhsIntegral2D      SparseMW_rhsintegral_x(SparseMW_basis2d, SepFunc1, refsol.deltas_x, no_deltas, order);
            SparseMW_SeparableRhsIntegral2D      SparseMW_rhsintegral_y(SparseMW_basis2d, SepFunc2, no_deltas, refsol.deltas_y, order);
            SparseMW_SumOfSeparableRhsIntegral2D SparseMW_rhsintegral2d(SparseMW_rhsintegral_x,SparseMW_rhsintegral_y);
            SparseMW_SumOfSeparableRhs           SparseMW_F(SparseMW_rhsintegral2d,SparseMW_P);
            SparseMW_F.readIndexSets(rhsfilename.str().c_str());

            SparseMW_GHS_ADWAV_SOLVER_SumofSeparableRhs SparseMW_ghs_adwav_solver(SparseMW_A,SparseMW_F,true,assemble_matrix);
            T alpha = 0.6, omega = 0.2, gamma = 0.15, theta = 2*omega/(1+omega);
            SparseMW_ghs_adwav_solver.setParameters(alpha, omega, gamma, theta);

            Coefficients<Lexicographical,T,Index2D> u;
            u = SparseMW_ghs_adwav_solver.SOLVE(SparseMW_F.norm_estimate, 1e-16, convfilename.str().c_str(),
                                                NumOfIterations, refsol.H1norm());
        }
    }
    return 0;
}

/*
IndexSet<Index1D>
computeRHSLambda_SingularPart(const MW_Basis1D &basis, const DenseVectorT &_f_singularPoints,
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
computeRHSLambda_SmoothPart(const MW_Basis1D &basis, T a, T b, int J_plus)
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
    RhsIntegral1D       u1_integral(mw_basis2d.first, u1Fct, nodeltas, 35);
    RhsIntegral1D       u2_integral(mw_basis2d.second, u2Fct, nodeltas, 35);

    Function<T> f1Fct(refsol.rhs_x, refsol.sing_pts_x);
    Function<T> f2Fct(refsol.rhs_y, refsol.sing_pts_y);
    RhsIntegral1D       f1_integral(mw_basis2d.first, f1Fct, refsol.deltas_x, 35);
    RhsIntegral1D       f2_integral(mw_basis2d.second, f2Fct, refsol.deltas_y, 35);

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
*/
