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
///     Special Multiwavelet Basis over an interval
typedef Basis<T,Primal,RPlus,SparseMulti>                             SparseMWBasis;

///     IdentityOperator in 1D, i.e. for $a(v,u) = \cdot \int(v \cdot u)$
typedef IdentityOperator1D<T,SparseMWBasis>                           IdentityOp;

///     LaplaceOperator in 1D, i.e. for $a(v,u) = \int(v_x \cdot u_x)$
typedef LaplaceOperator1D<T,SparseMWBasis>                            LaplaceOp;

///     ConvectionOperator in 1D, i.e. for $a(v,u) = \int(v \cdot u_x)$
typedef ConvectionOperator1D<T,SparseMWBasis>                         ConvectionOp;

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
    SparseMWBasis basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();


    /// Plot multi-scaling function and multi-wavelets
    int k_left_bspline  = basis.mra.rangeIL(j0).firstIndex();
    int k_right_bspline = 2*basis.mra.phi._numSplines;
    int k_left_wavelet  = basis.rangeJL(j0).firstIndex();
    int k_right_wavelet = 2*basis.psi._numSplines;

    T right_boundary_bspline = 1.;
    for (int k=k_left_bspline; k<=k_right_bspline; ++k) {
        right_boundary_bspline = std::max(right_boundary_bspline, basis.mra.phi.support(j0,k).l2);
    }
    T left_boundary_wavelet = 0., right_boundary_wavelet = 1.;
    for (int k=k_left_wavelet; k<=k_right_wavelet; ++k) {
        right_boundary_wavelet = std::max(right_boundary_wavelet, basis.psi.support(j0,k).l2);
    }

    ofstream plotfile_scaling("rplus_multiscaling.txt");
    for (T x=0.; x<=right_boundary_bspline; x+=pow2i<T>(-6)) {
        plotfile_scaling << x;
        for (int k=k_left_bspline; k<=k_right_bspline; ++k) {
            plotfile_scaling << " " << basis.generator(XBSpline)(x,j0,k,0);
        }
        plotfile_scaling << endl;
    }
    ofstream plotfile_wavelet("rplus_multiwavelet.txt");
    for (T x=0.; x<=right_boundary_wavelet; x+=pow2i<T>(-6)) {
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
    LaplaceOp    laplace_op(basis);
    ConvectionOp convection_op(basis);

    DenseVectorT singPts;
    Function<T> p0_Fct(p0,singPts);
    Function<T> p1_Fct(p1,singPts);
    Function<T> p2_Fct(p2,singPts);
    Function<T> p3_Fct(p3,singPts);
    IntegralF<Gauss,SparseMWBasis> sparsemw_p0(p0_Fct, basis);
    IntegralF<Gauss,SparseMWBasis> sparsemw_p1(p1_Fct, basis);
    IntegralF<Gauss,SparseMWBasis> sparsemw_p2(p2_Fct, basis);
    IntegralF<Gauss,SparseMWBasis> sparsemw_p3(p3_Fct, basis);


    /// Check for orthogonality and vanishing moments
    for (int k_row=k_left_bspline; k_row<=k_right_bspline; ++k_row) {
        for (int j_col=j0; j_col<=J-1; ++j_col) {
            for (int k_col=k_left_wavelet; k_col<=k_right_wavelet; ++k_col) {
                T identity_val = identity_op(XBSpline,j0,k_row,XWavelet,j_col,k_col);
                if (fabs(identity_val)>1e-15 && j_col-j0>=1) {
                    cout << "identity: (" << j0 << ", " << k_row << "), ("
                         << j_col << ", " << k_col << "): " << identity_val << endl;
                }
                T laplace_val = laplace_op(XBSpline,j0,k_row,XWavelet,j_col,k_col);
                if (fabs(laplace_val)>1e-15 && j_col-j0>=1) {
                    cout << "laplace: (" << j0 << ", " << k_row << "), ("
                         << j_col << ", " << k_col << "): " << laplace_val << endl;
                }
                T convection_val = convection_op(XBSpline,j0,k_row,XWavelet,j_col,k_col);
                if (fabs(convection_val)>1e-15 && j_col-j0>=1) {
                    cout << "convection: (" << j0 << ", " << k_row << "), ("
                         << j_col << ", " << k_col << "): " << convection_val << endl;
                }
            }
        }
    }
    for (int j_row=j0; j_row<=J-1; ++j_row) {
        for (int k_row=k_left_wavelet; k_row<=k_right_wavelet; ++k_row) {

            T val_p0 = sparsemw_p0(j_row,k_row,XWavelet,0);
            T val_p1 = sparsemw_p1(j_row,k_row,XWavelet,0);
            T val_p2 = sparsemw_p2(j_row,k_row,XWavelet,0);
            T val_p3 = sparsemw_p3(j_row,k_row,XWavelet,0);

            if (d>=1 && fabs(val_p0)>1e-15) {
                cout << "p0: (" << j_row << ", " << k_row << "): " << val_p0 << endl;
            }
            if (d>=2 && fabs(val_p1)>1e-15) {
                cout << "p1: (" << j_row << ", " << k_row << "): " << val_p1 << endl;
            }
            if (d>=3 && fabs(val_p2)>1e-15) {
                cout << "p2: (" << j_row << ", " << k_row << "): " << val_p2 << endl;
            }
            if (d>=4 && fabs(val_p3)>1e-15) {
                cout << "p3: (" << j_row << ", " << k_row << "): " << val_p3 << endl;
            }

            for (int k_col=k_left_bspline; k_col<=k_right_bspline; ++k_col) {
                T identity_val = identity_op(XWavelet,j_row,k_row,XBSpline,j0,k_col);
                if (fabs(identity_val)>1e-15 && j_row-j0>=1) {
                    cout << "identity: (" << j_row << ", " << k_row << "), ("
                         << j0 << ", " << k_col << "): " << identity_val << endl;
                }
            }
            for (int j_col=j0; j_col<=J-1; ++j_col) {
                for (int k_col=k_left_wavelet; k_col<=k_right_wavelet; ++k_col) {
                    T identity_val = identity_op(XWavelet,j_row,k_row,XWavelet,j_col,k_col);
                    if (fabs(identity_val)>1e-15 && fabs(j_row-j_col)>1) {
                        cout << "identity: (" << j_row << ", " << k_row << "), ("
                             << j_col << ", " << k_col << "): " << identity_val << endl;
                    }
                    T laplace_val = laplace_op(XWavelet,j_row,k_row,XWavelet,j_col,k_col);
                    if (fabs(laplace_val)>1e-15 && fabs(j_row-j_col)>1) {
                        cout << "laplace: (" << j_row << ", " << k_row << "), ("
                             << j_col << ", " << k_col << "): " << laplace_val << endl;
                    }
                    T convection_val = convection_op(XWavelet,j_row,k_row,XWavelet,j_col,k_col);
                    if (fabs(convection_val)>1e-15 && fabs(j_row-j_col)>1) {
                        cout << "convection: (" << j_row << ", " << k_row << "), ("
                             << j_col << ", " << k_col << "): " << convection_val << endl;
                    }
                }
            }
        }
    }

    return 0;
}
