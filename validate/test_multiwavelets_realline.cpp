/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <lawa/methods/adaptive/adaptive.h>

using namespace std;
using namespace lawa;

typedef double T;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<long double, cxxblas::ColMajor> >  DenseMatrixLongT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Typedefs for problem components:
///     Multiwavelet Basis over an interval
typedef Basis<T,Orthogonal,R,Multi>                                 MWBasis;

///     IdentityOperator in 1D, i.e. for $a(v,u) = \cdot \int(v \cdot u)$
typedef IdentityOperator1D<T, MWBasis>                              IdentityOp;

///     HelmholtzOperator in 1D, i.e. for $a(v,u) = \int(v_x \cdot u_x) + c \cdot \int(v \cdot u)$
typedef HelmholtzOperator1D<T, MWBasis>                             HelmholtzOp;

T
p0(T x) { return 1.; }

T
p1(T x) { return x; }

T
p2(T x) { return x*x; }

T
p3(T x) { return x*x*x; }


int main(int argc, char *argv[]) {
    if (argc != 4) {
        cout << "usage " << argv[0] << " d j0 J" << endl; exit(1);
    }
    /// wavelet basis parameters:
    int d  = atoi(argv[1]);          // d-wavelets
    int j0 = atoi(argv[2]);         // minimal level
    int J  = atoi(argv[3]);         // maximal level
    cout.precision(16);

    /// Basis initialization
    MWBasis basis(d,j0);
    BSpline<T,Orthogonal,R,Multi> phi(d);
    Wavelet<T,Orthogonal,R,Multi> psi(d);


    /// Plot multi-scaling function and multi-wavelets
    int k_left_bspline  = -basis.mra.phi._numSplines;
    int k_right_bspline =  basis.mra.phi._numSplines;
    int k_left_wavelet  = -basis.psi._numSplines;
    int k_right_wavelet =  basis.psi._numSplines;

    T left_boundary_bspline = 0., right_boundary_bspline = 1.;
    for (int k=k_left_bspline; k<=k_right_bspline; ++k) {
        left_boundary_bspline  = std::min(left_boundary_bspline, basis.mra.phi.support(j0,k).l1);
        right_boundary_bspline = std::max(right_boundary_bspline, basis.mra.phi.support(j0,k).l2);
    }
    T left_boundary_wavelet = 0., right_boundary_wavelet = 1.;
    for (int k=k_left_wavelet; k<=k_right_wavelet; ++k) {
        left_boundary_wavelet  = std::min(left_boundary_wavelet, basis.psi.support(j0,k).l1);
        right_boundary_wavelet = std::max(right_boundary_wavelet, basis.psi.support(j0,k).l2);
    }

    ofstream plotfile_scaling("realline_multiscaling.txt");
    for (T x=left_boundary_bspline; x<=right_boundary_bspline; x+=pow2i<T>(-6)) {
        plotfile_scaling << x;
        for (int k=k_left_bspline; k<=k_right_bspline; ++k) {
            plotfile_scaling << " " << basis.generator(XBSpline)(x,j0,k,0);
        }
        plotfile_scaling << endl;
    }
    ofstream plotfile_wavelet("realline_multiwavelet.txt");
    for (T x=left_boundary_wavelet; x<=right_boundary_wavelet; x+=pow2i<T>(-6)) {
        plotfile_wavelet << x;
        for (int j=j0; j<=J-1; ++j) {
            for (int k=k_left_wavelet; k<=k_right_wavelet; ++k) {
                plotfile_wavelet << " " << basis.generator(XWavelet)(x,j,k,0);
            }
        }
        plotfile_wavelet << endl;
    }


    /// Operator and integral initialization
    IdentityOp   identity_op(basis);

    DenseVectorT singPts;
    Function<T> p0_Fct(p0,singPts);
    Function<T> p1_Fct(p1,singPts);
    Function<T> p2_Fct(p2,singPts);
    Function<T> p3_Fct(p3,singPts);
    IntegralF<Gauss,MWBasis> mw_p0(p0_Fct, basis);
    IntegralF<Gauss,MWBasis> mw_p1(p1_Fct, basis);
    IntegralF<Gauss,MWBasis> mw_p2(p2_Fct, basis);
    IntegralF<Gauss,MWBasis> mw_p3(p3_Fct, basis);


    /// Check for orthogonality and vanishing moments
    int N = (k_right_bspline-k_left_bspline+1)+(J-j0)*(k_right_wavelet-k_left_wavelet+1);
    DenseMatrixT identity_A_dense(N,N);
    DenseVectorT p0(N), p1(N), p2(N), p3(N);

    int rowcount = 1;
    for (int k_row=k_left_bspline; k_row<=k_right_bspline; ++k_row) {
        int colcount=1;
        for (int k_col=k_left_bspline; k_col<=k_right_bspline; ++k_col) {
            identity_A_dense(rowcount,colcount) = identity_op(XBSpline,j0,k_row,XBSpline,j0,k_col);
            ++colcount;
        }
        for (int j_col=j0; j_col<=J-1; ++j_col) {
            for (int k_col=k_left_wavelet; k_col<=k_right_wavelet; ++k_col) {
                identity_A_dense(rowcount,colcount) = identity_op(XBSpline,j0,k_row,XWavelet,j_col,k_col);
                ++colcount;
            }
        }
        // By construction, scaling functions do not have vanishing moments (in general)
        p0(rowcount) = 0.;//mw_p0(j0,k_row,XBSpline,0);
        p1(rowcount) = 0.;//mw_p1(j0,k_row,XBSpline,0);
        p2(rowcount) = 0.;//mw_p2(j0,k_row,XBSpline,0);
        p3(rowcount) = 0.;//mw_p3(j0,k_row,XBSpline,0);
        ++rowcount;
    }
    for (int j_row=j0; j_row<=J-1; ++j_row) {
        for (int k_row=k_left_wavelet; k_row<=k_right_wavelet; ++k_row) {
            int colcount=1;
            for (int k_col=k_left_bspline; k_col<=k_right_bspline; ++k_col) {
                identity_A_dense(rowcount,colcount) = identity_op(XWavelet,j_row,k_row,XBSpline,j0,k_col);
                ++colcount;
            }
            for (int j_col=j0; j_col<=J-1; ++j_col) {
                for (int k_col=k_left_wavelet; k_col<=k_right_wavelet; ++k_col) {
                    identity_A_dense(rowcount,colcount) = identity_op(XWavelet,j_row,k_row,XWavelet,j_col,k_col);
                    ++colcount;
                }
            }
            p0(rowcount) = mw_p0(j_row,k_row,XWavelet,0);
            p1(rowcount) = mw_p1(j_row,k_row,XWavelet,0);
            p2(rowcount) = mw_p2(j_row,k_row,XWavelet,0);
            p3(rowcount) = mw_p3(j_row,k_row,XWavelet,0);
            ++rowcount;
        }
    }

    for (int i=1; i<=identity_A_dense.numRows(); ++i) {
        for (int j=1; j<=identity_A_dense.numCols(); ++j) {
            if (i==j) {
                if (fabs(identity_A_dense(i,j)-1.)>1e-14) {
                    cout << "(" << i << ", " << j << "): " << fabs(identity_A_dense(i,j)) << endl;
                }
            }
            else {
                if (fabs(identity_A_dense(i,j))>1e-14) {
                    cout << "(" << i << ", " << j << "): " << fabs(identity_A_dense(i,j)) << endl;
                }
            }
        }
    }
    cout << endl << endl;
    for (int i=1; i<=p0.length(); ++i) {
        if (fabs(p0(i)>1e-13) && d>=1) {
            cout << "p0: (" << i << "): " << p0(i) << endl;
        }
        if (fabs(p1(i)>1e-13) && d>=2) {
            cout << "p1: (" << i << "): " << p1(i) << endl;
        }
        if (fabs(p2(i)>1e-13) && d>=3) {
            cout << "p2: (" << i << "): " << p2(i) << endl;
        }
        if (fabs(p3(i)>1e-13) && d>=4) {
            cout << "p3: (" << i << "): " << p3(i) << endl;
        }
    }

/*
    /// Check for orthogonality
    Integral<Gauss,MWBasis,MWBasis> integral(basis,basis);
    cout << "Phi_1 vs. Phi1: " << integral(j0,1,XBSpline,0, j0,1,XBSpline,0) << endl;
    cout << "Phi_1 vs. Phi2: " << integral(j0,1,XBSpline,0, j0,2,XBSpline,0) << endl;
    cout << "Phi_1 vs. Phi3: " << integral(j0,1,XBSpline,0, j0,3,XBSpline,0) << endl;
    cout << "Phi_1 vs. Psi1: " << integral(j0,1,XBSpline,0, j0,1,XWavelet,0) << endl;
    cout << "Phi_1 vs. Psi2: " << integral(j0,1,XBSpline,0, j0,2,XWavelet,0) << endl;
    cout << "Phi_1 vs. Psi3: " << integral(j0,1,XBSpline,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Phi_2 vs. Phi1: " << integral(j0,2,XBSpline,0, j0,1,XBSpline,0) << endl;
    cout << "Phi_2 vs. Phi2: " << integral(j0,2,XBSpline,0, j0,2,XBSpline,0) << endl;
    cout << "Phi_2 vs. Phi3: " << integral(j0,2,XBSpline,0, j0,3,XBSpline,0) << endl;
    cout << "Phi_2 vs. Psi1: " << integral(j0,2,XBSpline,0, j0,1,XWavelet,0) << endl;
    cout << "Phi_2 vs. Psi2: " << integral(j0,2,XBSpline,0, j0,2,XWavelet,0) << endl;
    cout << "Phi_2 vs. Psi3: " << integral(j0,2,XBSpline,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Phi_3 vs. Phi1: " << integral(j0,3,XBSpline,0, j0,1,XBSpline,0) << endl;
    cout << "Phi_3 vs. Phi2: " << integral(j0,3,XBSpline,0, j0,2,XBSpline,0) << endl;
    cout << "Phi_3 vs. Phi3: " << integral(j0,3,XBSpline,0, j0,3,XBSpline,0) << endl;
    cout << "Phi_3 vs. Psi1: " << integral(j0,3,XBSpline,0, j0,1,XWavelet,0) << endl;
    cout << "Phi_3 vs. Psi2: " << integral(j0,3,XBSpline,0, j0,2,XWavelet,0) << endl;
    cout << "Phi_3 vs. Psi3: " << integral(j0,3,XBSpline,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Psi_1 vs. Phi1: " << integral(j0,1,XWavelet,0, j0,1,XBSpline,0) << endl;
    cout << "Psi_1 vs. Phi2: " << integral(j0,1,XWavelet,0, j0,2,XBSpline,0) << endl;
    cout << "Psi_1 vs. Phi3: " << integral(j0,1,XWavelet,0, j0,3,XBSpline,0) << endl;
    cout << "Psi_1 vs. Psi1: " << integral(j0,1,XWavelet,0, j0,1,XWavelet,0) << endl;
    cout << "Psi_1 vs. Psi2: " << integral(j0,1,XWavelet,0, j0,2,XWavelet,0) << endl;
    cout << "Psi_1 vs. Psi3: " << integral(j0,1,XWavelet,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Psi_2 vs. Phi1: " << integral(j0,2,XWavelet,0, j0,1,XBSpline,0) << endl;
    cout << "Psi_2 vs. Phi2: " << integral(j0,2,XWavelet,0, j0,2,XBSpline,0) << endl;
    cout << "Psi_2 vs. Phi3: " << integral(j0,2,XWavelet,0, j0,3,XBSpline,0) << endl;
    cout << "Psi_2 vs. Psi1: " << integral(j0,2,XWavelet,0, j0,1,XWavelet,0) << endl;
    cout << "Psi_2 vs. Psi2: " << integral(j0,2,XWavelet,0, j0,2,XWavelet,0) << endl;
    cout << "Psi_2 vs. Psi3: " << integral(j0,2,XWavelet,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Psi_3 vs. Phi1: " << integral(j0,3,XWavelet,0, j0,1,XBSpline,0) << endl;
    cout << "Psi_3 vs. Phi2: " << integral(j0,3,XWavelet,0, j0,2,XBSpline,0) << endl;
    cout << "Psi_3 vs. Phi3: " << integral(j0,3,XWavelet,0, j0,3,XBSpline,0) << endl;
    cout << "Psi_3 vs. Psi1: " << integral(j0,3,XWavelet,0, j0,1,XWavelet,0) << endl;
    cout << "Psi_3 vs. Psi2: " << integral(j0,3,XWavelet,0, j0,2,XWavelet,0) << endl;
    cout << "Psi_3 vs. Psi3: " << integral(j0,3,XWavelet,0, j0,3,XWavelet,0) << endl << endl;

    /// Test for vanishing moments
    DenseVectorT singPts;
    Function<T> p0_Fct(p0,singPts);
    Function<T> p1_Fct(p1,singPts);
    Function<T> p2_Fct(p2,singPts);

    IntegralF<Gauss,MWBasis> mw_p0(p0_Fct, basis);
    IntegralF<Gauss,MWBasis> mw_p1(p1_Fct, basis);
    IntegralF<Gauss,MWBasis> mw_p2(p2_Fct, basis);


    for (int k=1; k<=(int)phi._numSplines; ++k) {
        cout << "<Phi_" << k << ",const>=" << mw_p0(j0,k,XBSpline,0) << endl;
        cout << "<Phi_" << k << ",x>    =" << mw_p1(j0,k,XBSpline,0) << endl;
        cout << "<Phi_" << k << ",x*x>  =" << mw_p2(j0,k,XBSpline,0) << endl <<endl;

        cout << "<d_Phi_" << k << ",const>=" << mw_p0(j0,k,XBSpline,1) << endl;
        cout << "<d_Phi_" << k << ",x>    =" << mw_p1(j0,k,XBSpline,1) << endl;
        cout << "<d_Phi_" << k << ",x*x>  =" << mw_p2(j0,k,XBSpline,1) << endl<<endl;
    }
    for (int k=1; k<=(int)psi._numSplines; ++k) {
        cout << "<Psi_" << k << ",const>=" << mw_p0(j0,k,XWavelet,0) << endl;
        cout << "<Psi_" << k << ",x>    =" << mw_p1(j0,k,XWavelet,0) << endl;
        cout << "<Psi_" << k << ",x*x>  =" << mw_p2(j0,k,XWavelet,0) << endl <<endl;

        cout << "<d_Psi_" << k << ",const>=" << mw_p0(j0,k,XWavelet,1) << endl;
        cout << "<d_Psi_" << k << ",x>    =" << mw_p1(j0,k,XWavelet,1) << endl;
        cout << "<d_Psi_" << k << ",x*x>  =" << mw_p2(j0,k,XWavelet,1) << endl<<endl;
    }
*/

    return 0;
}
