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
typedef IndexSet<Index2D>::const_iterator                               const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator         const_coeff2d_it;
typedef Coefficients<AbsoluteValue,T,Index2D>::const_iterator           const_coeff2d_abs_it;

//Basis definitions
typedef Basis<T,Primal,R,CDF>                                           CDF_Basis1D;
typedef TensorBasis2D<Adaptive, CDF_Basis1D,CDF_Basis1D>                CDF_Basis2D;
typedef Basis<T,Orthogonal,R,Multi>                                     MW_Basis1D;
typedef TensorBasis2D<Adaptive, MW_Basis1D,MW_Basis1D>                  MW_Basis2D;
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;
typedef TensorBasis2D<Adaptive, SparseMW_Basis1D,SparseMW_Basis1D>      SparseMW_Basis2D;


//Operator definitions
typedef HelmholtzOperator2D<T, CDF_Basis2D>                             CDF_HelmholtzOp2D;
typedef DiagonalMatrixPreconditioner2D<T,CDF_Basis2D,
                                       CDF_HelmholtzOp2D >              CDF_Prec;
typedef AdaptiveHelmholtzOperator2D<T,CDF_Basis2D,CDF_Prec>             CDF_MA;

typedef AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,R,Multi,
                                               Orthogonal,R,Multi>      MW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D,MW_MA>         MW_Prec;

typedef AdaptiveHelmholtzOperatorOptimized2D<T,Primal,R,SparseMulti,
                                               Primal,R,SparseMulti>    SparseMW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D,SparseMW_MA>   SparseMW_Prec;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,CDF_Basis2D >                                  CDF_SeparableRhsIntegral2D;
typedef SeparableRHS2D<T,MW_Basis2D >                                   MW_SeparableRhsIntegral2D;
typedef SeparableRHS2D<T,SparseMW_Basis2D >                             SparseMW_SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,CDF_SeparableRhsIntegral2D,
                             CDF_SeparableRhsIntegral2D>                CDF_SumOfSeparableRhsIntegral2D;
typedef SumOfTwoRHSIntegrals<T,Index2D,MW_SeparableRhsIntegral2D,
                             MW_SeparableRhsIntegral2D>                 MW_SumOfSeparableRhsIntegral2D;
typedef SumOfTwoRHSIntegrals<T,Index2D,SparseMW_SeparableRhsIntegral2D,
                             SparseMW_SeparableRhsIntegral2D>           SparseMW_SumOfSeparableRhsIntegral2D;

typedef RHS<T,Index2D, CDF_SumOfSeparableRhsIntegral2D,
            CDF_Prec>                                                   CDF_SumOfSeparableRhs;
typedef RHS<T,Index2D, MW_SumOfSeparableRhsIntegral2D,
            MW_Prec>                                                    MW_SumOfSeparableRhs;
typedef RHS<T,Index2D,SparseMW_SumOfSeparableRhsIntegral2D,
            SparseMW_Prec>                                              SparseMW_SumOfSeparableRhs;

//Righthandsides definitions (non-separable)
typedef SmoothRHSWithAlignedSing2D<T, CDF_Basis2D, SparseGridGP>        CDF_NonSeparableRhsIntegralSG2D;
typedef SmoothRHSWithAlignedSing2D<T, MW_Basis2D, SparseGridGP>         MW_NonSeparableRhsIntegralSG2D;
typedef SmoothRHSWithAlignedSing2D<T, SparseMW_Basis2D, SparseGridGP>   SparseMW_NonSeparableRhsIntegralSG2D;

typedef SmoothRHSWithAlignedSing2D<T, CDF_Basis2D, FullGridGL>          CDF_NonSeparableRhsIntegralFG2D;
typedef SmoothRHSWithAlignedSing2D<T, MW_Basis2D, FullGridGL>           MW_NonSeparableRhsIntegralFG2D;
typedef SmoothRHSWithAlignedSing2D<T, SparseMW_Basis2D, FullGridGL>     SparseMW_NonSeparableRhsIntegralFG2D;

typedef SumOfThreeRHSIntegrals<T, Index2D,
                               CDF_NonSeparableRhsIntegralFG2D>         CDF_SumOfNonSeparableRhsIntegral2D;
typedef SumOfThreeRHSIntegrals<T, Index2D,
                               MW_NonSeparableRhsIntegralFG2D>          MW_SumOfNonSeparableRhsIntegral2D;
typedef SumOfThreeRHSIntegrals<T, Index2D,
                               SparseMW_NonSeparableRhsIntegralFG2D>    SparseMW_SumOfNonSeparableRhsIntegral2D;

typedef RHS<T,Index2D, CDF_NonSeparableRhsIntegralSG2D,
            CDF_Prec>                                                   CDF_NonSeparableRhs2D;
typedef RHS<T,Index2D, MW_NonSeparableRhsIntegralSG2D,
            MW_Prec>                                                    MW_NonSeparableRhs2D;
typedef RHS<T,Index2D, SparseMW_NonSeparableRhsIntegralSG2D,
            SparseMW_Prec>                                              SparseMW_NonSeparableRhs2D;

typedef RHS<T,Index2D, CDF_SumOfNonSeparableRhsIntegral2D,
            CDF_Prec>                                                   CDF_SumOfNonSeparableRhs2D;
typedef RHS<T,Index2D, MW_SumOfNonSeparableRhsIntegral2D,
            MW_Prec>                                                    MW_SumOfNonSeparableRhs2D;
typedef RHS<T,Index2D, SparseMW_SumOfNonSeparableRhsIntegral2D,
            SparseMW_Prec>                                              SparseMW_SumOfNonSeparableRhs2D;

//Algorithm definition
typedef S_ADWAV<T,Index2D,CDF_Basis2D,CDF_MA,CDF_SumOfSeparableRhs>     CDF_S_ADWAV_SOLVER_SeparableRhs;
typedef S_ADWAV<T,Index2D,MW_Basis2D, MW_MA,MW_SumOfSeparableRhs>       MW_S_ADWAV_SOLVER_SeparableRhs;
typedef S_ADWAV<T,Index2D,SparseMW_Basis2D,SparseMW_MA,
                SparseMW_SumOfSeparableRhs>                             SparseMW_S_ADWAV_SOLVER_SeparableRhs;

typedef S_ADWAV<T,Index2D,CDF_Basis2D,CDF_MA,
                CDF_NonSeparableRhs2D>                                  CDF_S_ADWAV_SOLVER_NonSeparableRhs;
typedef S_ADWAV<T,Index2D,MW_Basis2D,MW_MA,
                MW_NonSeparableRhs2D>                                   MW_S_ADWAV_SOLVER_NonSeparableRhs;
typedef S_ADWAV<T,Index2D,SparseMW_Basis2D,SparseMW_MA,
                SparseMW_NonSeparableRhs2D>                             SparseMW_S_ADWAV_SOLVER_NonSeparableRhs;

typedef S_ADWAV<T,Index2D,CDF_Basis2D,CDF_MA,
                CDF_SumOfNonSeparableRhs2D>                             CDF_S_ADWAV_SOLVER_SumNonSeparableRhs;
typedef S_ADWAV<T,Index2D,MW_Basis2D,MW_MA,
                MW_SumOfNonSeparableRhs2D>                              MW_S_ADWAV_SOLVER_SumNonSeparableRhs;
typedef S_ADWAV<T,Index2D,MW_Basis2D,MW_MA,
                MW_SumOfNonSeparableRhs2D>                              SparseMW_S_ADWAV_SOLVER_SumNonSeparableRhs;

int main (int argc, char *argv[]) {
    if (argc!=8) {
        cout << "usage " << argv[0] << " basistype d d_ jmin_x jmin_y example max_its" << endl; exit(1);
    }
    cout.precision(3);

    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0_x=atoi(argv[4]);
    int j0_y=atoi(argv[5]);
    int example=atoi(argv[6]);
    int NumOfIterations=atoi(argv[7]);

    T c=1.; //for other values of c, on the fly error computation does not work!!
    T contraction = 0.125;
    T threshTol = 0.4;
    T cgTol = 0.1*threshTol;//1e-12;
    T resTol=1e-4;

    IndexSet<Index2D> InitialLambda;
    Index1D index_x(j0_x,0,XBSpline);
    Index1D index_y(j0_y,0,XBSpline);
    InitialLambda.insert(Index2D(index_x,index_y));

    stringstream convfilename;
    convfilename << "s_adwav_conv_realline_helmholtz2d_" << argv[1] << "_" << argv[2] << "_"
                 << argv[3] << "_" << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";

    //Righthand side construction for tensor solution
    if (strcmp(argv[1],"CDF")==0) {
        CDF_Basis1D       CDF_basis_x(d,d_,j0_x);
        CDF_Basis1D       CDF_basis_y(d,d_,j0_y);
        CDF_Basis2D       CDF_basis2d(CDF_basis_x,CDF_basis_y);
        CDF_HelmholtzOp2D CDF_Bil(CDF_basis2d, c);
        CDF_Prec          CDF_P(CDF_Bil);
        CDF_MA            CDF_A(CDF_basis2d, c, CDF_P, 0., 8092, 8092);


        if (example==1 || example==2 || example==3) {
            int order = 35;
            TensorRefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(example, c, 0. ,0. ,1.);
            SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                            refsol.exact_y, refsol.sing_pts_y);

            SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                            refsol.rhs_y, refsol.sing_pts_y);
            GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;

            CDF_SeparableRhsIntegral2D CDF_rhsintegral_x(CDF_basis2d, SepFunc1, refsol.deltas_x,
                                                                     no_deltas, order);
            CDF_SeparableRhsIntegral2D CDF_rhsintegral_y(CDF_basis2d, SepFunc2, no_deltas,
                                                         refsol.deltas_y, order);
            CDF_SumOfSeparableRhsIntegral2D CDF_rhsintegral2d(CDF_rhsintegral_x,CDF_rhsintegral_y);
            CDF_SumOfSeparableRhs CDF_F(CDF_rhsintegral2d,CDF_P);

            CDF_S_ADWAV_SOLVER_SeparableRhs CDF_s_adwav_solver(CDF_basis2d, CDF_A, CDF_F, contraction,
                                                         threshTol, cgTol, resTol, NumOfIterations,
                                                         2, 1e-2,1000000);
            CDF_s_adwav_solver.solve(InitialLambda, "cg", convfilename.str().c_str(), 2,
                                     refsol.H1norm());

            plot2D<T,CDF_Basis2D,CDF_Prec>(CDF_basis2d, CDF_s_adwav_solver.solutions[NumOfIterations-1], CDF_P, refsol.exact,
                   -80., 80, -80., 80., pow2i<T>(-4), "example2");
        }
        else if (example==4) {
            int order = 7;
            RefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(1, 1.);
            Function2D<T> Func2d(refsol.rhs, refsol.sing_pts_x, refsol.sing_pts_y);
            CDF_NonSeparableRhsIntegralSG2D CDF_rhsintegral2d(CDF_basis2d, Func2d, order);
            CDF_NonSeparableRhs2D           CDF_F(CDF_rhsintegral2d,CDF_P);

            CDF_S_ADWAV_SOLVER_NonSeparableRhs CDF_s_adwav_solver(CDF_basis2d, CDF_A, CDF_F, contraction,
                                                                  threshTol, cgTol, resTol, NumOfIterations,
                                                                  2, 1e-2,1000000);
            CDF_s_adwav_solver.solve(InitialLambda, "cg", convfilename.str().c_str(), 2,
                                                 refsol.H1norm());

            plot2D<T,CDF_Basis2D,CDF_Prec>(CDF_basis2d, CDF_s_adwav_solver.solutions[NumOfIterations-1], CDF_P, refsol.exact,
                   -2.5, 2.5, -2.5, 2.5, pow2i<T>(-4), "example2");
        }
        else if (example==5) {
            RefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(2, 1.);
            Function2D<T> Func2d(refsol.exact, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_x(refsol.exact_dx, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_y(refsol.exact_dy, refsol.sing_pts_x, refsol.sing_pts_y);

            if (d==2) {
                int order = 20;
                CDF_NonSeparableRhsIntegralFG2D CDF_rhsintegral_reaction(CDF_basis2d, Func2d, order);
                CDF_NonSeparableRhsIntegralFG2D CDF_rhsintegral_diffusion_x(CDF_basis2d, Func2d_x, order, 1, 0);
                CDF_NonSeparableRhsIntegralFG2D CDF_rhsintegral_diffusion_y(CDF_basis2d, Func2d_y, order, 0, 1);
                CDF_SumOfNonSeparableRhsIntegral2D CDF_rhsintegral2d(CDF_rhsintegral_diffusion_x,
                                                                       CDF_rhsintegral_diffusion_y,
                                                                       CDF_rhsintegral_reaction);
                CDF_SumOfNonSeparableRhs2D CDF_F(CDF_rhsintegral2d,CDF_P);
                CDF_S_ADWAV_SOLVER_SumNonSeparableRhs CDF_s_adwav_solver(CDF_basis2d, CDF_A, CDF_F, contraction,
                                                                threshTol, cgTol, resTol, NumOfIterations,
                                                                2, 1e-2,1000000);
                CDF_s_adwav_solver.solve(InitialLambda, "cg", convfilename.str().c_str(), 2,
                                         refsol.H1norm());

                plot2D<T,CDF_Basis2D,CDF_Prec>(CDF_basis2d, CDF_s_adwav_solver.solutions[NumOfIterations-1],
                                               CDF_P, refsol.exact, -4., 4, -4., 4., pow2i<T>(-4), "example5");
            }
        }
    }
    else if (strcmp(argv[1],"MW")==0) {
        MW_Basis1D       MW_basis_x(d,j0_x);
        MW_Basis1D       MW_basis_y(d,j0_y);
        MW_Basis2D       MW_basis2d(MW_basis_x,MW_basis_y);
        MW_MA            MW_A(MW_basis2d, c);
        MW_Prec          MW_P(MW_A);

        if (example==1 || example==2 || example==3) {
            int order = 35;
            TensorRefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(example, 1., 0., 0., 1.);
            SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                            refsol.exact_y, refsol.sing_pts_y);

            SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                            refsol.rhs_y, refsol.sing_pts_y);
            GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;

            MW_SeparableRhsIntegral2D MW_rhsintegral_x(MW_basis2d, SepFunc1, refsol.deltas_x,
                                                         no_deltas, order);
            MW_SeparableRhsIntegral2D MW_rhsintegral_y(MW_basis2d, SepFunc2, no_deltas,
                                                         refsol.deltas_y, order);
            MW_SumOfSeparableRhsIntegral2D MW_rhsintegral2d(MW_rhsintegral_x,MW_rhsintegral_y);
            MW_SumOfSeparableRhs MW_F(MW_rhsintegral2d,MW_P);

            MW_S_ADWAV_SOLVER_SeparableRhs MW_s_adwav_solver(MW_basis2d, MW_A, MW_F, contraction,
                                                         threshTol, cgTol, resTol, NumOfIterations,
                                                         2, 1e-2,1000000);
            MW_s_adwav_solver.solve(InitialLambda, "cg", convfilename.str().c_str(), 2,
                                    refsol.H1norm());

            stringstream rhsfilename;
            rhsfilename << "rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                        << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";

            Coefficients<Lexicographical,T,Index2D> f;
            IndexSet<Index2D> Lambda;
            Lambda = supp(MW_s_adwav_solver.solutions[NumOfIterations-1]);

            IndexSet<Index2D> Extension;
            Extension = C(Lambda,1.,MW_basis2d);
            Lambda = Lambda + Extension;
            f = MW_F(Lambda);
            Coefficients<AbsoluteValue,T,Index2D> f_abs;
            f_abs = f;
            cout << f.norm(2.) << " " << f_abs.norm(2.) << endl;

            ofstream rhsfile(rhsfilename.str().c_str());
            rhsfile << f.norm(2.) << endl;
            for (int k=0; k<=30; ++k) {
                T eta=pow(2.,(T)-k);
                f = MW_F(eta);
                cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

                IndexSet<Index2D> supp_f;
                supp_f = supp(f);
                rhsfile << "#," << eta << endl;
                for (const_set2d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
                    if (Lambda.count(*it)>0) {
                        Lambda.erase(*it);
                        rhsfile << *it << endl;
                    }
                }
                rhsfile << endl;
            }
            rhsfile.close();

            stringstream plot_filename;
            plot_filename << "s-adwav-realline-helmholtz2d-plot_" << example << "_" << d << "_" << d_
                          << "_" << j0_x << "_" << j0_y;
            cout << "Plot of solution started." << endl;
            plot2D(MW_basis2d, MW_s_adwav_solver.solutions[NumOfIterations-1], MW_P, refsol.exact,
                   -40., 40., -40., 40., pow2i<T>(-3), plot_filename.str().c_str());
            cout << "Plot of solution finished." << endl;
        }
        else if (example==4) {
            int order = 27;
            RefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(1, 1.);
            Function2D<T> Func2d(refsol.rhs, refsol.sing_pts_x, refsol.sing_pts_y);
            MW_NonSeparableRhsIntegralSG2D MW_rhsintegral2d(MW_basis2d, Func2d, order);
            MW_NonSeparableRhs2D           MW_F(MW_rhsintegral2d,MW_P);

            MW_S_ADWAV_SOLVER_NonSeparableRhs MW_s_adwav_solver(MW_basis2d, MW_A, MW_F, contraction,
                                                                threshTol, cgTol, resTol, NumOfIterations,
                                                                2, 1e-2,1000000);
            MW_s_adwav_solver.solve(InitialLambda, "cg", convfilename.str().c_str(), 2,
                                                 refsol.H1norm());

            stringstream rhsfilename;
            rhsfilename << "rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                        << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";

            Coefficients<Lexicographical,T,Index2D> f;
            IndexSet<Index2D> Lambda;
            Lambda = supp(MW_s_adwav_solver.solutions[NumOfIterations-1]);

            IndexSet<Index2D> Extension;
            Extension = C(Lambda,4.,MW_basis2d);
            Lambda = Lambda + Extension;
            f = MW_F(Lambda);
            Coefficients<AbsoluteValue,T,Index2D> f_abs;
            f_abs = f;
            cout << f.norm(2.) << " " << f_abs.norm(2.) << endl;

            ofstream rhsfile(rhsfilename.str().c_str());
            rhsfile << f.norm(2.) << endl;
            for (int k=0; k<=30; ++k) {
                T eta=pow(2.,(T)-k);
                f = MW_F(eta);
                cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

                IndexSet<Index2D> supp_f;
                supp_f = supp(f);
                rhsfile << "#," << eta << endl;
                for (const_set2d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
                    if (Lambda.count(*it)>0) {
                        Lambda.erase(*it);
                        rhsfile << *it << endl;
                    }
                }
                rhsfile << endl;
            }
            rhsfile.close();

            plot2D<T,MW_Basis2D,MW_Prec>(MW_basis2d, MW_s_adwav_solver.solutions[NumOfIterations-1], MW_P, refsol.exact,
                   -2.5, 2.5, -2.5, 2.5, pow2i<T>(-4), "example4");
        }
        else if (example==5) {
            RefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(2, 1.);
            Function2D<T> Func2d(refsol.exact, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_x(refsol.exact_dx, refsol.sing_pts_x, refsol.sing_pts_y);
            Function2D<T> Func2d_y(refsol.exact_dy, refsol.sing_pts_x, refsol.sing_pts_y);

            if (d==2) {
                int order = 20;
                MW_NonSeparableRhsIntegralFG2D       MW_rhsintegral_reaction(MW_basis2d, Func2d, order);
                MW_NonSeparableRhsIntegralFG2D       MW_rhsintegral_diffusion_x(MW_basis2d, Func2d_x, order, 1, 0);
                MW_NonSeparableRhsIntegralFG2D       MW_rhsintegral_diffusion_y(MW_basis2d, Func2d_y, order, 0, 1);
                MW_SumOfNonSeparableRhsIntegral2D    MW_rhsintegral2d(MW_rhsintegral_diffusion_x,
                                                                      MW_rhsintegral_diffusion_y,
                                                                      MW_rhsintegral_reaction);
                MW_SumOfNonSeparableRhs2D            MW_F(MW_rhsintegral2d,MW_P);
                MW_S_ADWAV_SOLVER_SumNonSeparableRhs MW_s_adwav_solver(MW_basis2d, MW_A, MW_F, contraction,
                                                              threshTol, cgTol, resTol, NumOfIterations,
                                                              3, 1e-2,1000000);
                MW_s_adwav_solver.solve(InitialLambda, "cg", convfilename.str().c_str(), 2,
                                         refsol.H1norm());

                plot2D<T,MW_Basis2D,MW_Prec>(MW_basis2d, MW_s_adwav_solver.solutions[NumOfIterations-1],
                                             MW_P, refsol.exact, -4., 4, -4., 4., pow2i<T>(-4), "example5");

                stringstream rhsfilename;
               rhsfilename << "rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                           << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";

               Coefficients<Lexicographical,T,Index2D> f;
               IndexSet<Index2D> Lambda;
               Lambda = supp(MW_s_adwav_solver.solutions[NumOfIterations-1]);

               IndexSet<Index2D> Extension;
               Extension = C(Lambda,4.,MW_basis2d);
               Lambda = Lambda + Extension;
               Extension = C(Lambda,4.,MW_basis2d);
               Lambda = Lambda + Extension;
               f = MW_F(Lambda);
               Coefficients<AbsoluteValue,T,Index2D> f_abs;
               f_abs = f;
               cout << f.norm(2.) << " " << f_abs.norm(2.) << endl;

               ofstream rhsfile(rhsfilename.str().c_str());
               rhsfile << f.norm(2.) << endl;
               for (int k=0; k<=30; ++k) {
                   T eta=pow(2.,(T)-k);
                   f = MW_F(eta);
                   cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

                   IndexSet<Index2D> supp_f;
                   supp_f = supp(f);
                   rhsfile << "#," << eta << endl;
                   for (const_set2d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
                       if (Lambda.count(*it)>0) {
                           Lambda.erase(*it);
                           rhsfile << *it << endl;
                       }
                   }
                   rhsfile << endl;
               }
               rhsfile.close();
            }


        }
    }
    else if (strcmp(argv[1],"SparseMW")==0) {
        SparseMW_Basis1D       SparseMW_basis_x(d,j0_x);
        SparseMW_Basis1D       SparseMW_basis_y(d,j0_y);
        SparseMW_Basis2D       SparseMW_basis2d(SparseMW_basis_x,SparseMW_basis_y);
        SparseMW_MA            SparseMW_A(SparseMW_basis2d, c);
        SparseMW_Prec          SparseMW_P(SparseMW_A);

        if (example==1 || example==2 || example==3) {
            int order = 50;
            TensorRefSols_PDE_Realline2D<T> refsol;
            refsol.setExample(example, 1., 0., 0., 1.);
            SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                            refsol.exact_y, refsol.sing_pts_y);

            SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                            refsol.rhs_y, refsol.sing_pts_y);
            GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;

            SparseMW_SeparableRhsIntegral2D SparseMW_rhsintegral_x(SparseMW_basis2d, SepFunc1, refsol.deltas_x,
                                                         no_deltas, order);
            SparseMW_SeparableRhsIntegral2D SparseMW_rhsintegral_y(SparseMW_basis2d, SepFunc2, no_deltas,
                                                         refsol.deltas_y, order);
            SparseMW_SumOfSeparableRhsIntegral2D SparseMW_rhsintegral2d(SparseMW_rhsintegral_x,SparseMW_rhsintegral_y);
            SparseMW_SumOfSeparableRhs SparseMW_F(SparseMW_rhsintegral2d,SparseMW_P);
            SparseMW_S_ADWAV_SOLVER_SeparableRhs SparseMW_s_adwav_solver(SparseMW_basis2d, SparseMW_A, SparseMW_F, contraction,
                                                         threshTol, cgTol, resTol, NumOfIterations,
                                                         1, 1e-2,1000000);
            SparseMW_s_adwav_solver.solve(InitialLambda, "cg", convfilename.str().c_str(), 0,
                                          refsol.H1norm());

            stringstream rhsfilename;
                        rhsfilename << "rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                                    << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";

            Coefficients<Lexicographical,T,Index2D> f;
            IndexSet<Index2D> Lambda;
            Lambda = supp(SparseMW_s_adwav_solver.solutions[NumOfIterations-1]);

            IndexSet<Index2D> Extension;
            Extension = C(Lambda,0.25,SparseMW_basis2d);
            Lambda = Lambda + Extension;
            std::cerr << "  Size of enlarged Lambda: " << Lambda.size() << std::endl;
            Extension = C(Lambda,0.25,SparseMW_basis2d);
            Lambda = Lambda + Extension;
            std::cerr << "  Size of enlarged Lambda: " << Lambda.size() << std::endl;
//            Extension = C(Lambda,4.,SparseMW_basis2d);
//            Lambda = Lambda + Extension;
            f = SparseMW_F(Lambda);
            Coefficients<AbsoluteValue,T,Index2D> f_abs;
            f_abs = f;
            cout << f.norm(2.) << " " << f_abs.norm(2.) << endl;

            ofstream rhsfile(rhsfilename.str().c_str());
            rhsfile << f.norm(2.) << endl;
            for (int k=0; k<=50; ++k) {
                T eta=pow(2.,(T)-k);
                f = SparseMW_F(eta);
                cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

                IndexSet<Index2D> supp_f;
                supp_f = supp(f);
                rhsfile << "#," << eta << endl;
                for (const_set2d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
                    if (Lambda.count(*it)>0) {
                        Lambda.erase(*it);
                        rhsfile << *it << endl;
                    }
                }
                rhsfile << endl;
            }
            rhsfile.close();

            stringstream plot_filename;
            plot_filename << "s-adwav-realline-helmholtz2d-plot_" << example << "_" << d << "_" << d_
                          << "_" << j0_x << "_" << j0_y;
            cout << "Plot of solution started." << endl;
            plot2D(SparseMW_basis2d, SparseMW_s_adwav_solver.solutions[NumOfIterations-1], SparseMW_P, refsol.exact, -10., 10., -10., 10.,
                   pow2i<T>(-3), plot_filename.str().c_str());
            cout << "Plot of solution finished." << endl;
        }
    }
    else {
        std::cerr << "Not yet implemented for " << argv[1] << std::endl;
        return 0;
    }

    return 0;
}

