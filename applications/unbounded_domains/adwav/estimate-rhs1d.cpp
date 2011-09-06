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
#include <vector>
#include <lawa/lawa.h>
#include <applications/unbounded_domains/referencesolutions/refsols_pde_realline1d.h>
#include <applications/unbounded_domains/parameters/parameters.h>

typedef double T;
using namespace lawa;
using namespace std;

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
typedef Basis<T,Primal,R,CDF>                                           CDF_Basis1D;
typedef Basis<T,Orthogonal,R,Multi>                                     MW_Basis1D;
typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;

//Operator definitions
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,CDF>            CDF_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,R,Multi>      MW_MA;
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Primal,R,SparseMulti>    SparseMW_MA;

typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, CDF_MA>       CDF_Prec;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, MW_MA>        MW_Prec;
typedef DiagonalPreconditionerAdaptiveOperator<T,Index1D, SparseMW_MA>  SparseMW_Prec;

typedef RHSWithPeaks1D<T, CDF_Basis1D>                                  CDF_RhsIntegral1D;
typedef RHSWithPeaks1D_WO_XBSpline<T>                                   CDF_RhsIntegral1D_WO_XBSpline;
typedef RHSWithPeaks1D<T, MW_Basis1D>                                   MW_RhsIntegral1D;
typedef RHSWithPeaks1D<T, SparseMW_Basis1D>                             SparseMW_RhsIntegral1D;

typedef RHS<T,Index1D, CDF_RhsIntegral1D, CDF_Prec>                     CDF_Rhs_Ref;
typedef RHS<T,Index1D, CDF_RhsIntegral1D_WO_XBSpline, CDF_Prec>         CDF_Rhs_Ref_WO_XBSpline;
typedef RHS<T,Index1D, MW_RhsIntegral1D, MW_Prec>                       MW_Rhs_Ref;
typedef RHS<T,Index1D, SparseMW_RhsIntegral1D, SparseMW_Prec>           SparseMW_Rhs_Ref;

typedef RHS1D<T, CDF_RhsIntegral1D, CDF_Prec>                           CDF_Rhs1D;
typedef RHS1D<T, CDF_RhsIntegral1D_WO_XBSpline, CDF_Prec>               CDF_Rhs1D_WO_XBSpline;
typedef RHS1D<T, MW_RhsIntegral1D, MW_Prec>                             MW_Rhs1D;
typedef RHS1D<T, SparseMW_RhsIntegral1D, SparseMW_Prec>                 SparseMW_Rhs1D;

IndexSet<Index1D>
computeRHSLambda_SingularPart(const CDF_Basis1D &basis, const DenseVectorT &_f_singularPoints,
                              int J_plus);

IndexSet<Index1D>
computeRHSLambda_SmoothPart(const CDF_Basis1D &basis, T a, T b, int J_plus);

IndexSet<Index1D>
computeRHSLambda_SingularPart_WO_XBSpline(const CDF_Basis1D &basis,
                                          const DenseVectorT &_f_singularPoints,
                                          int J_plus, int J_minus);

IndexSet<Index1D>
computeRHSLambda_SmoothPart_WO_XBSpline(const CDF_Basis1D &basis, T a, T b, int J_plus,
                                        int J_minus);

IndexSet<Index1D>
computeRHSLambda_SingularPart(const MW_Basis1D &basis,
                              const DenseVectorT &_f_singularPoints, int J_plus);

IndexSet<Index1D>
computeRHSLambda_SmoothPart(const MW_Basis1D &basis, T a, T b, int J_plus);

IndexSet<Index1D>
computeRHSLambda_SingularPart(const SparseMW_Basis1D &basis, const DenseVectorT &_f_singularPoints,
                              int J_plus);

IndexSet<Index1D>
computeRHSLambda_SmoothPart(const SparseMW_Basis1D &basis, T a, T b, int J_plus);

int main (int argc, char *argv[]) {

    if (argc!=6) {
        cout << "usage " << argv[0] << " basistype d d_ jmin example" << endl; exit(1);
    }
    cout.precision(20);

    int d=atoi(argv[2]);
    int d_=atoi(argv[3]);
    int j0; bool w_XBSpline;
    T c = 1.;
    int example=atoi(argv[5]);

    RefSols_PDE_Realline1D<T> refsol;
    refsol.setExample(example,1.,0.,1.);
    Function<T>        rhsFct(refsol.rhs,refsol.sing_pts);

    stringstream rhsfilename;
    if      (strcmp(argv[4],"-inf")==0) {
        j0=0;             w_XBSpline=false;
        rhsfilename << "rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_"
                    << argv[3] << "_" << argv[4] << "_" << c << "_" << argv[5] << ".dat";
    }
    else if (strcmp(argv[4],"best")==0) {
        j0=refsol.getMinimalLevel(d,d_); w_XBSpline=true;
        rhsfilename << "rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_"
                    << argv[3] << "_" << j0 << "_" << c << "_" << argv[5] << ".dat";
    }
    else {
        j0=atoi(argv[4]); w_XBSpline=true;
        rhsfilename << "rhs_realline_helmholtz_" << argv[1] << "_" << argv[2] << "_"
                    << argv[3] << "_" << argv[4] << "_" << c << "_" << argv[5] << ".dat";
    }

    T left_bound = 0., right_bound = 0.;
    int J_plus_smooth = 0, J_plus_singular = 0, J_minus_smooth=0, J_minus_singular=0;
    bool singular_integral=false;


    Coefficients<Lexicographical,T,Index1D> f;

    if (strcmp(argv[1],"CDF")==0) {
        CDF_Basis1D             CDF_basis(d,d_,j0);
        CDF_MA                  CDF_A(CDF_basis,w_XBSpline,c);
        CDF_Prec                CDF_prec(CDF_A);

        if (w_XBSpline) {
            refsol.getRHS_W_XBSplineParameters(CDF_basis.d, CDF_basis.d_, left_bound, right_bound,
                                               J_plus_smooth, J_plus_singular, singular_integral);

            CDF_RhsIntegral1D       CDF_rhsintegral1d(CDF_basis, rhsFct, refsol.deltas, 120);
            CDF_Rhs1D               CDF_F(CDF_rhsintegral1d,CDF_prec);
            CDF_Rhs_Ref             CDF_F_Ref(CDF_rhsintegral1d,CDF_prec);

            IndexSet<Index1D> Lambda, Lambda_singular, Lambda_smooth;
            Lambda_singular =  computeRHSLambda_SingularPart
                                (CDF_basis,refsol.sing_pts, J_plus_singular);
            Lambda_smooth   =  computeRHSLambda_SmoothPart
                                (CDF_basis, left_bound, right_bound, J_plus_smooth);
            Lambda = Lambda_smooth + Lambda_singular;

            f = CDF_F_Ref(Lambda);
            Coefficients<AbsoluteValue,T,Index1D> f_abs;
            f_abs = f;
            cout << f.norm(2.) << " " << f_abs.norm(2.) << endl;


            ofstream rhsfile(rhsfilename.str().c_str());
            rhsfile << f.norm(2.) << endl;
            for (int k=0; k<=50; ++k) {
                T eta=pow(2.,(T)-k);
                f = CDF_F_Ref(eta);
                cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

                IndexSet<Index1D> supp_f;
                supp_f = supp(f);
                rhsfile << "#," << eta << endl;
                for (const_set1d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
                    if (Lambda.count(*it)>0) {
                        Lambda.erase(*it);
                        rhsfile << *it << endl;
                    }
                }
                rhsfile << endl;
            }
            rhsfile.close();

            if (CDF_F.readIndexSets(rhsfilename.str().c_str()) ) {
                cout << "Index sets for rhs read... Ready to start."  << endl;
            }
            else {
                cout << "RHS: Could not open file." << endl;
                return 0;
            }

            Coefficients<Lexicographical,T,Index1D> f_eta, diff;

            for (T eta=0.1; eta>1e-8; eta*=0.5) {
                f_eta = CDF_F(eta);
                diff = f-f_eta;
                cout << "|| f - f_thresh ||_2 = " << diff.norm(2.) << " (should be " << eta << ")" << endl;
                cout << "f norm: " << CDF_F.rhs_data.norm(2.) << endl;
                f_eta = CDF_F_Ref(eta);
                diff = f-f_eta;
                cout << "|| f - f_thresh_Ref ||_2 = " << diff.norm(2.) << " (should be " << eta << ")" << endl;
                cout << "f norm: " << CDF_F_Ref.rhs_data.norm(2.) << endl << endl;;
            }
        }
        else {
            refsol.getRHS_WO_XBSplineParameters(CDF_basis.d, CDF_basis.d_, left_bound, right_bound,
                                                J_plus_smooth, J_minus_smooth,
                                                J_plus_singular, J_minus_singular, singular_integral);

            CDF_RhsIntegral1D_WO_XBSpline    CDF_rhsintegral1d_WO_XBSpline
                                               (CDF_basis.psi, refsol.rhs, refsol.sing_pts,
                                                refsol.deltas, left_bound, right_bound,
                                                1., 20);
            CDF_Rhs1D_WO_XBSpline            CDF_F_WO_XBSpline(CDF_rhsintegral1d_WO_XBSpline,
                                                               CDF_prec);
            CDF_Rhs_Ref_WO_XBSpline          CDF_F_Ref_WO_XBSpline(CDF_rhsintegral1d_WO_XBSpline,
                                                                   CDF_prec);

            IndexSet<Index1D> Lambda, Lambda_singular, Lambda_smooth;
            Lambda_singular =  computeRHSLambda_SingularPart_WO_XBSpline
                                (CDF_basis,refsol.sing_pts, J_plus_singular, J_minus_singular);
            Lambda_smooth   =  computeRHSLambda_SmoothPart_WO_XBSpline
                                (CDF_basis, left_bound, right_bound, J_plus_smooth, J_minus_smooth);
            Lambda = Lambda_smooth + Lambda_singular;

            f = CDF_F_Ref_WO_XBSpline(Lambda);
            ofstream rhsfile(rhsfilename.str().c_str());
            rhsfile << f.norm(2.) << endl;
            for (int k=0; k<=50; ++k) {
                T eta=pow(2.,(T)-k);
                f = CDF_F_Ref_WO_XBSpline(eta);
                cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

                IndexSet<Index1D> supp_f;
                supp_f = supp(f);
                rhsfile << "#," << eta << endl;
                for (const_set1d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
                    if (Lambda.count(*it)>0) {
                        Lambda.erase(*it);
                        rhsfile << *it << endl;
                    }
                }
                rhsfile << endl;
            }
            rhsfile.close();


            if (CDF_F_WO_XBSpline.readIndexSets(rhsfilename.str().c_str()) ) {
                cout << "Index sets for rhs read... Ready to start."  << endl;
            }
            else {
                cout << "RHS: Could not open file." << endl;
                return 0;
            }
        }
    }
    else if (strcmp(argv[1],"MW")==0) {
        MW_Basis1D              MW_basis(d,j0);
        MW_MA                   MW_A(MW_basis,c);
        MW_Prec                 MW_prec(MW_A);
        MW_RhsIntegral1D        MW_rhsintegral1d(MW_basis, rhsFct, refsol.deltas, 20);
        MW_Rhs1D                MW_F(MW_rhsintegral1d,MW_prec);
        MW_Rhs_Ref              MW_F_Ref(MW_rhsintegral1d,MW_prec);

        refsol.getRHS_W_XBSplineParameters(MW_basis.d, MW_basis.d, left_bound, right_bound,
                                           J_plus_smooth, J_plus_singular, singular_integral);

        IndexSet<Index1D> Lambda, Lambda_singular, Lambda_smooth;
        Lambda_singular =  computeRHSLambda_SingularPart
                            (MW_basis,refsol.sing_pts, J_plus_singular);
        Lambda_smooth   =  computeRHSLambda_SmoothPart
                            (MW_basis, left_bound, right_bound, J_plus_smooth);
        Lambda = Lambda_smooth + Lambda_singular;

        cout << "MW: #Lambda = " << Lambda.size() << endl;

        f = MW_F_Ref(Lambda);

        ofstream rhsfile(rhsfilename.str().c_str());
        rhsfile << f.norm(2.) << endl;
        for (int k=0; k<=50; ++k) {
            T eta=pow(2.,(T)-k);
            f = MW_F_Ref(eta);
            cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

            IndexSet<Index1D> supp_f;
            supp_f = supp(f);
            rhsfile << "#," << eta << endl;
            for (const_set1d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
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
        Coefficients<Lexicographical,T,Index1D> f_eta, diff;
        Coefficients<AbsoluteValue,T,Index1D> diff_abs;

        for (T eta=0.1; eta>1e-8; eta*=0.5) {
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
    else if (strcmp(argv[1],"SparseMW")==0) {
        SparseMW_Basis1D            SparseMW_basis(d,j0);
        SparseMW_MA                 SparseMW_A(SparseMW_basis,c);
        SparseMW_Prec               SparseMW_prec(SparseMW_A);
        SparseMW_RhsIntegral1D      SparseMW_rhsintegral1d(SparseMW_basis, rhsFct, refsol.deltas, 40);
        SparseMW_Rhs1D              SparseMW_F(SparseMW_rhsintegral1d,SparseMW_prec);
        SparseMW_Rhs_Ref            SparseMW_F_Ref(SparseMW_rhsintegral1d,SparseMW_prec);

        refsol.getRHS_W_XBSplineParameters(SparseMW_basis.d, -1, left_bound, right_bound,
                                           J_plus_smooth, J_plus_singular, singular_integral);

        IndexSet<Index1D> Lambda, Lambda_singular, Lambda_smooth;
        Lambda_singular =  computeRHSLambda_SingularPart
                            (SparseMW_basis,refsol.sing_pts, J_plus_singular);
        Lambda_smooth   =  computeRHSLambda_SmoothPart
                            (SparseMW_basis, left_bound, right_bound, J_plus_smooth);
        Lambda = Lambda_smooth + Lambda_singular;

        cout << "SparseMW: #Lambda = " << Lambda.size() << endl;

        f = SparseMW_F_Ref(Lambda);


        ofstream rhsfile(rhsfilename.str().c_str());
        rhsfile << f.norm(2.) << endl;
        for (int k=0; k<=50; ++k) {
            T eta=pow(2.,(T)-k);
            f = SparseMW_F_Ref(eta);
            cout << "Size of index set for eta = " << eta  << ": " << f.size() << endl;

            IndexSet<Index1D> supp_f;
            supp_f = supp(f);
            rhsfile << "#," << eta << endl;
            for (const_set1d_it it=supp_f.begin(); it!=supp_f.end(); ++it) {
                if (Lambda.count(*it)>0) {
                    Lambda.erase(*it);
                    rhsfile << *it << endl;
                }
            }
            rhsfile << endl;
        }
        rhsfile.close();

        if (SparseMW_F.readIndexSets(rhsfilename.str().c_str()) ) {
            cout << "Index sets for rhs read... Ready to start."  << endl;
        }
        else {
            cout << "RHS: Could not open file." << endl;
            return 0;
        }

        cout.precision(16);
        Coefficients<Lexicographical,T,Index1D> f_eta, diff;

        for (T eta=0.1; eta>1e-8; eta*=0.5) {
            f_eta = SparseMW_F(eta);
            diff = f-f_eta;
            cout << "|| f - f_thresh ||_2 = " << diff.norm(2.) << " (should be " << eta << ")" << endl;
            cout << "f norm: " << SparseMW_F.rhs_data.norm(2.) << endl;
            f_eta = SparseMW_F_Ref(eta);
            diff = f-f_eta;
            cout << "|| f - f_thresh_Ref ||_2 = " << diff.norm(2.) << " (should be " << eta << ")" << endl;
            cout << "f norm: " << SparseMW_F_Ref.rhs_data.norm(2.) << endl << endl;;
        }
    }



    return 0;
}


IndexSet<Index1D>
computeRHSLambda_SingularPart(const CDF_Basis1D &basis, const DenseVectorT &_f_singularPoints,
                              int J_plus)
{
    IndexSet<Index1D> ret;
    T l1, l2;
    l1 = basis.mra.phi.support(0,0).l1, l2 = basis.mra.phi.support(0,0).l2;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        int k_left =  std::floor(float(pow2i<T>(basis.j0)*x-l2));
        int k_right = std::ceil(float(pow2i<T>(basis.j0)*x-l1));
        for (int k=k_left; k<=k_right; ++k) {
            ret.insert(Index1D(basis.j0,k,XBSpline));
        }
    }

    l1 = basis.psi.support(0,0).l1, l2 = basis.psi.support(0,0).l2;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        for (int j=basis.j0; j<=J_plus; ++j) {
            int k_left =  std::floor(float(pow2i<T>(j)*x-l2));
            int k_right = std::ceil(float(pow2i<T>(j)*x-l1));
            for (int k=k_left; k<=k_right; ++k) {
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }
    return ret;
}

IndexSet<Index1D>
computeRHSLambda_SmoothPart(const CDF_Basis1D &basis, T a, T b, int J_plus)
{
    IndexSet<Index1D> ret;
    T l1, l2;
    l1 = basis.mra.phi.support(0,0).l1, l2 = basis.mra.phi.support(0,0).l2;
    int k_left =  std::floor(float(pow2i<T>(basis.j0)*a-l2));
    int k_right = std::ceil(float(pow2i<T>(basis.j0)*b-l1));
    for (int k=k_left; k<=k_right; ++k) {
        ret.insert(Index1D(basis.j0,k,XBSpline));
    }

    l1 = basis.psi.support(0,0).l1, l2 = basis.psi.support(0,0).l2;
    for (int j=basis.j0; j<=J_plus; ++j) {
        int k_left =  std::floor(float(pow2i<T>(j)*a-l2));
        int k_right = std::ceil(float(pow2i<T>(j)*b-l1));
        for (int k=k_left; k<=k_right; ++k) {
            ret.insert(Index1D(j,k,XWavelet));
        }
    }
    return ret;
}


IndexSet<Index1D>
computeRHSLambda_SingularPart_WO_XBSpline(const CDF_Basis1D &basis,
                                          const DenseVectorT &_f_singularPoints,
                                          int J_plus, int J_minus)
{
    T l1 = basis.psi.support(0,0).l1, l2 = basis.psi.support(0,0).l2;
    IndexSet<Index1D> ret;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        for (int j=0; j<=J_plus; ++j) {
            int k_left =  std::floor(float(pow2i<T>(j)*x-l2));
            int k_right = std::ceil(float(pow2i<T>(j)*x-l1));
            for (int k=k_left; k<=k_right; ++k) {
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
        for (int j=-1; j>=J_minus; --j) {
            int k_left =  std::floor(float(pow2i<T>(j)*x-l2));
            int k_right = std::ceil(float(pow2i<T>(j)*x-l1));
            for (int k=k_left; k<=k_right; ++k) {
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }
    return ret;
}

IndexSet<Index1D>
computeRHSLambda_SmoothPart_WO_XBSpline(const CDF_Basis1D &basis, T a, T b, int J_plus, int J_minus)
{
    T l1 = basis.psi.support(0,0).l1, l2 = basis.psi.support(0,0).l2;
    IndexSet<Index1D> ret;
    for (int j=0; j<=J_plus; ++j) {
        int k_left =  std::floor(float(pow2i<T>(j)*a-l2));
        int k_right = std::ceil(float(pow2i<T>(j)*b-l1));
        for (int k=k_left; k<=k_right; ++k) {
            ret.insert(Index1D(j,k,XWavelet));
        }
    }
    for (int j=-1; j>=J_minus; --j) {
        int k_left =  std::floor(float(pow2i<T>(j)*a-l2));
        int k_right = std::ceil(float(pow2i<T>(j)*b-l1));
        for (int k=k_left; k<=k_right; ++k) {
            ret.insert(Index1D(j,k,XWavelet));
        }
    }
    return ret;
}


IndexSet<Index1D>
computeRHSLambda_SingularPart(const SparseMW_Basis1D &basis, const DenseVectorT &_f_singularPoints,
                              int J_plus)
{
    IndexSet<Index1D> ret;

    T l1, l2;
    l1 = basis.mra.phi.max_support().l1, l2 = basis.mra.phi.max_support().l2;
    int numSplines =basis.mra.phi._numSplines;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        int kMin = floor( pow2i<T>(basis.j0)*x - l2) - 3;
        int kMax =  ceil( pow2i<T>(basis.j0)*x - l1) + 3;

        kMin = kMin*numSplines;
        kMax = kMax*numSplines;

        for (int k=kMin*numSplines; k<=kMax*numSplines; ++k) {
            ret.insert(Index1D(basis.j0,k,XBSpline));
        }
    }

    l1 = basis.psi.max_support().l1, l2 = basis.psi.max_support().l2;
    int numWavelets =basis.psi._numSplines;
    for (int i=1; i<=_f_singularPoints.length(); ++i) {
        T x=_f_singularPoints(i);
        for (int j=basis.j0; j<=J_plus; ++j) {
            int kMin = floor( pow2i<T>(j)*x - l2) / 2 - 3;
            int kMax =  ceil( pow2i<T>(j)*x - l1) / 2 + 3;

            kMin = kMin*numWavelets;
            kMax = kMax*numWavelets;

            for (int k=kMin; k<=kMax; ++k) {
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }

    return ret;
}

IndexSet<Index1D>
computeRHSLambda_SmoothPart(const SparseMW_Basis1D &basis, T a, T b, int J_plus)
{
    IndexSet<Index1D> ret;

    T l1, l2;
    l1 = basis.mra.phi.max_support().l1, l2 = basis.mra.phi.max_support().l2;
    int numSplines =basis.mra.phi._numSplines;

    int kMin = floor( pow2i<T>(basis.j0)*a - l2) - 3;
    int kMax =  ceil( pow2i<T>(basis.j0)*b - l1) + 3;

    kMin = kMin*numSplines;
    kMax = kMax*numSplines;

    for (int k=kMin*numSplines; k<=kMax*numSplines; ++k) {
        ret.insert(Index1D(basis.j0,k,XBSpline));
    }

    l1 = basis.psi.support(0,0).l1, l2 = basis.psi.support(0,0).l2;
    int numWavelets =basis.psi._numSplines;
    for (int j=basis.j0; j<=J_plus; ++j) {
        int kMin = floor( pow2i<T>(j)*a - l2) / 2 - 3;
        int kMax =  ceil( pow2i<T>(j)*b - l1) / 2 + 3;

        kMin = kMin*numWavelets;
        kMax = kMax*numWavelets;

        for (int k=kMin; k<=kMax; ++k) {
            ret.insert(Index1D(j,k,XWavelet));
        }
    }
    return ret;
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


