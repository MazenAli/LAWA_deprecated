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

template<typename T>
IndexSet<Index1D>
ReferenceLambda(int jmin, int jmax, const MWBasis1D &basis, T radius);


int main(int argc, char *argv[]) {
    int d=2;
    int j0=-2;
    int J=atoi(argv[1]);
    T   r = 30.;
    T   c = 1.;

    MWBasis1D           mwbasis1d(d,j0);
    MWBasis2D           mwbasis2d(mwbasis1d,mwbasis1d);

    MW_MA               mw_A(mwbasis2d,1.);
    MW_Prec             mw_prec(mw_A);

    // Assemble f vector and u vector
    DenseMatrixT nodeltas;
    TensorRefSols_PDE_Realline2D<T> refsol;
    refsol.setExample(1, c);
    Function<T> u1Fct(refsol.exact_x, refsol.sing_pts_x);
    Function<T> u2Fct(refsol.exact_y, refsol.sing_pts_y);
    RhsIntegral1D       u1_integral(mwbasis2d.first, u1Fct, nodeltas, 35);
    RhsIntegral1D       u2_integral(mwbasis2d.second, u2Fct, nodeltas, 35);

    Function<T> f1Fct(refsol.rhs_x, refsol.sing_pts_x);
    Function<T> f2Fct(refsol.rhs_y, refsol.sing_pts_y);
    RhsIntegral1D       f1_integral(mwbasis2d.first, f1Fct, nodeltas, 35);
    RhsIntegral1D       f2_integral(mwbasis2d.second, f2Fct, nodeltas, 35);

    IndexSet<Index1D> Lambda1d = ReferenceLambda(j0,J,mwbasis1d,r);
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

    cout << "Computation of 2d vector u and f started." << endl;

    Coefficients<Lexicographical,T,Index2D> u_coeff, f_coeff;

    for (const_coeff1d_it it_x=u1_coeff.begin(); it_x!=u1_coeff.end(); ++it_x) {
        for (const_coeff1d_it it_y=u2_coeff.begin(); it_y!=u2_coeff.end(); ++it_y) {
            //if ((*it_x).first.j-j0+(*it_y).first.j-j0 <= J) {
                Index2D index((*it_x).first,(*it_y).first);
                T u_val = (*it_x).second * (*it_y).second * (1./mw_prec(index));
                if (fabs(u_val)>1e-12) {
                    u_coeff[index] = u_val;
                }
                T f_val = (  f1_coeff[(*it_x).first] * (*it_y).second
                           + f2_coeff[(*it_y).first] * (*it_x).second ) * mw_prec(index);
                if (fabs(f_val)>1e-12) {
                    f_coeff[index] = f_val;
                }
            //}
        }
    }
    //u_coeff=THRESH(u_coeff,1e-10);
    cout << "Computation of 2d vectors u and f finished." << endl;
    cout << "Size of 2d u vector: " << u_coeff.size() << endl;
    cout << "Size of 2d f vector: " << f_coeff.size() << endl;


    Coefficients<AbsoluteValue,T,Index2D> u_coeff_abs;
    u_coeff_abs = u_coeff;


    for(int n=1; n<=20000; n+=1000) {
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


    f_coeff =  mw_A.apply(u_coeff, 21);
    cout << "|| f ||_{ell2}: " << f_coeff.norm(2.) << endl;

    Coefficients<Lexicographical,T,Index2D> Apply_u_coeff, Au_mv_coeff, diff;
    //Au_mv_coeff = mw_A.mv(supp(u_coeff), u_coeff);
    //cout << "|| Au_mv ||_{ell2}: " << Au_mv_coeff.norm(2.) << endl;

    for (int s=0; s<=20; ++s) {
        Timer time;
        time.start();
        Apply_u_coeff = mw_A.apply(u_coeff, s);
        time.stop();
        diff = Apply_u_coeff - f_coeff;
        //for (const_coeff2d_it it=Au_mv_coeff.begin(); it!=Au_mv_coeff.end(); ++it) {
        //    diff[(*it).first] = (*it).second - Apply_u_coeff[(*it).first];
        //}
        T ell2error = diff.norm(2.);
        cout << s << " " << ell2error << " " << time.elapsed() << " "
             << mw_A.findK(u_coeff_abs,ell2error) << endl;
    }


    return 0.;
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


/*
    Coefficients<AbsoluteValue,T,Index2D> u_coeff_abs;
    u_coeff_abs = u_coeff;
    ofstream file_bestnterm("bestnterm_u.dat");
    for (int n=1; n<=10000; n+=100) {
        file_bestnterm << n << " " << u_coeff_abs.l2bestnterm(n) << endl;
    }
    file_bestnterm.close();
*/

/*
    stringstream plot_filename;
    plot_filename << "sol_applyhelmholtzmw2d.dat";
    cout << "Plot of solution started." << endl;
    Coefficients<Lexicographical,T,Index2D> u_coeff_thresh;
    u_coeff_thresh = THRESH(u_coeff,1e-1);
    cout << "Size of threshold u: " << u_coeff_thresh.size() << endl;
    plot2D(mwbasis2d, u_coeff_thresh, prec, refsol.exact, -2., 2., -2., 2.,
           pow2i<T>(-3), plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;
*/


/*
    T thresh=1e-10;
    T ell2normf1_sq = 0.;
    for (const_coeff1d_it it_x=f1_coeff.begin(); it_x!=f1_coeff.end(); ++it_x) {
        for (const_coeff1d_it it_y=u2_coeff.begin(); it_y!=u2_coeff.end(); ++it_y) {
            Index2D index((*it_x).first,(*it_y).first);
            ell2normf1_sq +=   std::pow((*it_x).second * (*it_y).second * mw_prec(index) ,2.);
        }
    }
    cout << "ell2normf1_sq = " << ell2normf1_sq << endl;

    T ell2normf2_sq = 0.;
    for (const_coeff1d_it it_x=u1_coeff.begin(); it_x!=u1_coeff.end(); ++it_x) {
        for (const_coeff1d_it it_y=f2_coeff.begin(); it_y!=f2_coeff.end(); ++it_y) {
            Index2D index((*it_x).first,(*it_y).first);
            ell2normf2_sq +=   std::pow((*it_x).second * (*it_y).second * mw_prec(index) ,2.);
        }
    }
    cout << "ell2normf2_sq = " << ell2normf2_sq << endl;


    Coefficients<Lexicographical,T,Index2D> f_coeff;
    for (const_coeff1d_abs_it it_x=f1_coeff_abs.begin(); it_x!=f1_coeff_abs.end(); ++it_x) {
        bool stop = false;
        for (const_coeff1d_abs_it it_y=u2_coeff_abs.begin(); it_y!=u2_coeff_abs.end(); ++it_y) {
            Index2D index((*it_x).second,(*it_y).second);
            T tmp = (*it_x).first * (*it_y).first * mw_prec(index);
            if (fabs(tmp)<1e-15) break;
            ell2normf1_sq -= tmp*tmp;
            f_coeff[index] += tmp;
            if (ell2normf1_sq < thresh*thresh) {
                cout << "ell2normf1_sq finished" << endl;
                stop = true;
                break;
            }
        }
        if (stop) break;
    }

    for (const_coeff1d_abs_it it_x=u1_coeff_abs.begin(); it_x!=u1_coeff_abs.end(); ++it_x) {
        bool stop = false;
        for (const_coeff1d_abs_it it_y=f2_coeff_abs.begin(); it_y!=f2_coeff_abs.end(); ++it_y) {
            Index2D index((*it_x).second,(*it_y).second);
            T tmp = (*it_x).first * (*it_y).first * mw_prec(index);
            if (fabs(tmp)<1e-15) break;
            ell2normf2_sq -= tmp*tmp;
            f_coeff[index] += tmp;
            if (ell2normf2_sq < thresh*thresh) {
                cout << "ell2normf2_sq finished" << endl;
                stop = true;
                break;
            }
        }
        if (stop) break;
    }
    //f_coeff=THRESH(f_coeff,thresh);

    cout << "Size of 2d vector f: " << f_coeff.size() << endl;

    Coefficients<Lexicographical,T,Index2D> u_coeff;
    for (const_coeff2d_it it=f_coeff.begin(); it!=f_coeff.end(); ++it) {
        Index2D index((*it).first.index1,(*it).first.index2);
        u_coeff[(*it).first] = u1_coeff[(*it).first.index1]*u2_coeff[(*it).first.index2]/mw_prec(index);
    }

    cout << "Size of 2d vector u: " << u_coeff.size() << endl;
*/

