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
typedef TensorBasis2D<Adaptive,MW_Basis1D,MW_Basis1D>                   MW_Basis2D;

//Operator definitions
typedef NoPreconditioner<T,Index1D>                                     NoPrec1D;
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,R,Multi,
                                               Orthogonal,R,Multi>      MW_MA;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index2D, MW_MA>        MW_Prec;

//Righthandsides definitions (tensor)
typedef RHSWithPeaks1D<T,MW_Basis1D>                                    MW_RhsIntegral1D;
typedef RHS<T,Index1D,MW_RhsIntegral1D,NoPrec1D>                        MW_Rhs_Ref1D;
typedef SeparableRHS2D<T,MW_Basis2D>                                    MW_SeparableRhsIntegral2D;
typedef SumOfTwoRHSIntegrals<T,Index2D,MW_SeparableRhsIntegral2D,
                             MW_SeparableRhsIntegral2D>                 MW_SumOfSeparableRhsIntegral2D;
typedef RHS<T,Index2D,MW_SumOfSeparableRhsIntegral2D,
            MW_Prec>                                                    MW_Rhs_Ref;
typedef RHS2D<T,MW_SumOfSeparableRhsIntegral2D,MW_Prec>                 MW_Rhs;

IndexSet<Index1D>
computeRHSLambda_SingularPart(const MW_Basis1D &basis,
                              const DenseVectorT &_f_singularPoints, int J_plus);

IndexSet<Index1D>
computeRHSLambda_SmoothPart(const MW_Basis1D &basis, T a, T b, int J_plus);


int main (int argc, char *argv[]) {

    if (argc!=7) {
        cout << "usage " << argv[0] << " basistype d d_ j0_x j0_y example" << endl; exit(1);
    }
    cout.precision(20);

    int d   =atoi(argv[2]);
    int d_  =atoi(argv[3]);
    int j0_x=atoi(argv[4]);
    int j0_y=atoi(argv[5]);
    T c = 1.;
    int example=atoi(argv[6]);

    stringstream rhsfilename;
    rhsfilename << "rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_"
                << argv[4] << "_" << argv[5] << "_" << c << "_" << argv[6] << ".dat";

    NoPrec1D      NoP;

    //Righthand side construction for tensor solution
    if (example==1 || example==2 || example==3) {
        TensorRefSols_PDE_Realline2D<T> refsol;
        refsol.setExample(example, c);
        SeparableFunction2D<T> SepFunc1(refsol.rhs_x, refsol.sing_pts_x,
                                        refsol.exact_y, refsol.sing_pts_y);

        SeparableFunction2D<T> SepFunc2(refsol.exact_x, refsol.sing_pts_x,
                                        refsol.rhs_y, refsol.sing_pts_y);
        GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > no_deltas;

        int order=35;
        T r_x =  40.;
        T r_y =  40.;
        int J_smooth = 6, J_singular=40;

        if (strcmp(argv[1],"MW")==0) {
            MW_Basis1D       MW_basis_x(d,j0_x);
            MW_Basis1D       MW_basis_y(d,j0_y);
            MW_Basis2D       MW_basis2d(MW_basis_x,MW_basis_y);
            MW_MA            MW_A(MW_basis2d, c);
            MW_Prec          MW_P(MW_A);

            Function<T> u1Fct(refsol.exact_x, refsol.sing_pts_x);
            Function<T> u2Fct(refsol.exact_y, refsol.sing_pts_y);
            MW_RhsIntegral1D       u1_integral(MW_basis2d.first, u1Fct, no_deltas, order);
            MW_Rhs_Ref1D           u1_rhs(u1_integral,NoP);
            MW_RhsIntegral1D       u2_integral(MW_basis2d.second, u2Fct, no_deltas, order);
            MW_Rhs_Ref1D           u2_rhs(u2_integral,NoP);

            Function<T> f1Fct(refsol.rhs_x, refsol.sing_pts_x);
            Function<T> f2Fct(refsol.rhs_y, refsol.sing_pts_y);
            MW_RhsIntegral1D       f1_integral(MW_basis2d.first, f1Fct, refsol.deltas_x, order);
            MW_Rhs_Ref1D           f1_rhs(f1_integral,NoP);
            MW_RhsIntegral1D       f2_integral(MW_basis2d.second, f2Fct, refsol.deltas_y, order);
            MW_Rhs_Ref1D           f2_rhs(f2_integral,NoP);

            IndexSet<Index1D> Lambda1d, Lambda_smooth1d, Lambda_singular1d;
            Coefficients<Lexicographical,T,Index1D> u1_coeff, u2_coeff, f1_coeff, f2_coeff;

            Lambda_smooth1d = computeRHSLambda_SmoothPart(MW_basis2d.first, -r_x, r_x, J_smooth);
            Lambda_singular1d = computeRHSLambda_SingularPart(MW_basis2d.first,refsol.sing_pts_x, J_singular);
            Lambda1d = Lambda_smooth1d + Lambda_singular1d;
            u1_coeff = u1_rhs(Lambda1d);
            f1_coeff = f1_rhs(Lambda1d);

            Lambda_smooth1d = computeRHSLambda_SmoothPart(MW_basis2d.second, -r_y, r_y, J_smooth);
            Lambda_singular1d = computeRHSLambda_SingularPart(MW_basis2d.second,refsol.sing_pts_y, J_singular);
            Lambda1d = Lambda_smooth1d + Lambda_singular1d;
            u2_coeff = u2_rhs(Lambda1d);
            f2_coeff = f2_rhs(Lambda1d);

            cout << "Sizes of 1d u vectors before thresholding: " << u1_coeff.size() << ", " << u2_coeff.size() << endl;
            cout << "Sizes of 1d f vectors before thresholding: " << f1_coeff.size() << ", " << f2_coeff.size() << endl;

            u1_coeff = ABSOLUTE_THRESH(u1_coeff,1e-15);
            u2_coeff = ABSOLUTE_THRESH(u2_coeff,1e-15);
            f1_coeff = ABSOLUTE_THRESH(f1_coeff,1e-15);
            f2_coeff = ABSOLUTE_THRESH(f2_coeff,1e-15);

            cout << "Sizes of 1d u vectors after thresholding: " << u1_coeff.size() << ", " << u2_coeff.size() << endl;
            cout << "Sizes of 1d f vectors after thresholding: " << f1_coeff.size() << ", " << f2_coeff.size() << endl;

            cout << "Computation of 2d vector u and f started." << endl;

            IndexSet<Index2D> Lambda;
            Coefficients<Lexicographical,T,Index2D> f;


            for (const_coeff1d_it it_x=u1_coeff.begin(); it_x!=u1_coeff.end(); ++it_x) {
               for (const_coeff1d_it it_y=f2_coeff.begin(); it_y!=f2_coeff.end(); ++it_y) {
                   Index2D index((*it_x).first,(*it_y).first);
                   T f_val = (*it_x).second * (*it_y).second * MW_P(index);
                   if (fabs(f_val)>1e-12) {
                       f[index] += f_val;
                   }
               }
            }
            for (const_coeff1d_it it_x=f1_coeff.begin(); it_x!=f1_coeff.end(); ++it_x) {
               for (const_coeff1d_it it_y=u2_coeff.begin(); it_y!=u2_coeff.end(); ++it_y) {
                   Index2D index((*it_x).first,(*it_y).first);
                   T f_val = (*it_x).second * (*it_y).second * MW_P(index);
                   if (fabs(f_val)>1e-12) {
                       f[index] += f_val;
                   }
               }
            }

            Lambda = supp(f);
            cout << "Size of 2d vector f: " << f.size() << endl;


            MW_SeparableRhsIntegral2D      MW_rhsintegral_x(MW_basis2d, SepFunc1, refsol.deltas_x,
                                                            no_deltas, order);
            MW_SeparableRhsIntegral2D      MW_rhsintegral_y(MW_basis2d, SepFunc2, no_deltas,
                                                            refsol.deltas_y, order);
            MW_SumOfSeparableRhsIntegral2D MW_rhsintegral2d(MW_rhsintegral_x,MW_rhsintegral_y);
            MW_Rhs_Ref                     MW_F_Ref(MW_rhsintegral2d,MW_P,f);
            MW_Rhs_Ref                     MW_F_Ref2(MW_rhsintegral2d,MW_P);
            MW_Rhs                         MW_F(MW_rhsintegral2d,MW_P);

            Coefficients<AbsoluteValue,T,Index2D> f_abs;
            f_abs = f;
            cout << f.norm(2.) << " " << f_abs.norm(2.) << endl;

            ofstream rhsfile(rhsfilename.str().c_str());
            rhsfile << f.norm(2.) << endl;
            for (int k=0; k<=30; ++k) {
                T eta=pow(2.,(T)-k);
                f = MW_F_Ref(eta);
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

            if (MW_F.readIndexSets(rhsfilename.str().c_str()) ) {
                cout << "Index sets for rhs read... Ready to start."  << endl;
            }
            else {
                cout << "RHS: Could not open file." << endl;
                return 0;
            }

            cout.precision(16);
            Coefficients<Lexicographical,T,Index2D> f_eta, diff;
            Coefficients<AbsoluteValue,T,Index2D> diff_abs;

            for (T eta=1.1; eta>=0.01; eta*=0.5) {
                f_eta = MW_F(eta);
                diff = f-f_eta;
                diff_abs = diff;
                cout << "|| f - f_thresh ||_2 = " << diff_abs.norm(2.) << " (should be " << eta << ")" << endl;
                cout << "f norm: " << MW_F.rhs_data.norm(2.) << endl;

                f_eta = MW_F_Ref(eta);
                diff = f-f_eta;
                cout << "|| f - f_thresh_Ref ||_2 = " << diff.norm(2.) << " (should be " << eta << ")" << endl;
                cout << "f norm: " << MW_F_Ref.rhs_data.norm(2.) << endl << endl;;
            }
        }

    }




    return 0;
}


IndexSet<Index1D>
computeRHSLambda_SingularPart(const MW_Basis1D &basis, const DenseVector<Array<T> > &_f_singularPoints,
                              int J_plus)
{
    IndexSet<Index1D> ret;

    T l1, l2;
    l1 = basis.mra.phi.max_support().l1, l2 = basis.mra.phi.max_support().l2;
    int numScaling =basis.mra.phi._numSplines;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        int k_left =  std::floor(float(pow2i<T>(basis.j0)*x-l2));
        int k_right = std::ceil(float(pow2i<T>(basis.j0)*x-l1));
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
            int k_left =  std::floor(float(pow2i<T>(j)*x-l2));
            int k_right = std::ceil(float(pow2i<T>(j)*x-l1));
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



