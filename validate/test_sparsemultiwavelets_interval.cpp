/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Typedefs for problem components:

///     Special Multiwavelet Basis over an interval
typedef Basis<T,Primal,Interval,SparseMulti> SparseMWBasis;

///     IdentityOperator in 1D, i.e. for $a(v,u) = \cdot \int(v \cdot u)$
typedef IdentityOperator1D<T,SparseMWBasis>                           IdentityOp;

///     LaplaceOperator in 1D, i.e. for $a(v,u) = \int(v_x \cdot u_x)$
typedef LaplaceOperator1D<T,SparseMWBasis>                            LaplaceOp;

///     ConvectionOperator in 1D, i.e. for $a(v,u) = \int(v \cdot u_x)$
typedef ConvectionOperator1D<T,SparseMWBasis>                         ConvectionOp;

T
p0(T x) {   return 1.; }

T
p1(T x) {   return x; }

T
p2(T x) {   return x*x; }

T
p3(T x) {   return x*x*x; }

int main (int argc, char *argv[]) {

    cout.precision(16);
    int d=4;
    int j0=1;
    int j=1;
    int J=5;


    /// Basis initialization, using Dirichlet boundary conditions
    SparseMWBasis basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    /// Check for correct dimension of single and multi scale spaces
    int dim_Vj = basis.mra.cardI(j0);
    for (int j=j0; j<=J-1; ++j) {
        cout << "Range of W_" << j << ": " << basis.rangeJ(j) << endl;
        dim_Vj += basis.cardJ(j);
        if (dim_Vj != basis.mra.cardI(j+1)) {
            cout << "Incorrect dimensions: dim(multi scale space) = " << dim_Vj
                 << ", dim(single scale space) ) " << basis.mra.cardI(j+1) << endl;
        }
    }

    /// Plot multi-scaling function and multi-wavelets
    ofstream plotfile_scaling("interval_sparsemultiscaling.txt");
    for (T x=0.; x<=2.; x+=pow2i<T>(-12)) {
        plotfile_scaling << x;
        for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            plotfile_scaling << " " << basis.generator(XBSpline)(x,j0,k,0);
        }
        plotfile_scaling << endl;
    }
    ofstream plotfile_wavelet("interval_sparsemultiwavelet.txt");
    for (T x=0.; x<=1.; x+=pow2i<T>(-7)) {
        plotfile_wavelet << x;
        for (int j=j0; j<=J-1; ++j) {
            for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
                plotfile_wavelet << " " << basis.generator(XWavelet)(x,j,k,0);
            }
        }
        plotfile_wavelet << endl;
    }

    /// Operator initialization
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
    int N = basis.mra.cardI(J);

    for (int k_row=basis.mra.rangeI(j0).firstIndex(); k_row<=basis.mra.rangeI(j0).lastIndex(); ++k_row) {
        for (int j_col=j0; j_col<=J-1; ++j_col) {
            for (int k_col=basis.rangeJ(j_col).firstIndex(); k_col<=basis.rangeJ(j_col).lastIndex(); ++k_col) {
                T identity_val = identity_op(XBSpline,j0,k_row,XWavelet,j_col,k_col);
                if (fabs(identity_val)>1e-13 && j_col-j0>=1) {
                    cout << "identity: (" << j0 << ", " << k_row << "), ("
                         << j_col << ", " << k_col << "): " << identity_val << endl;
                }
                T laplace_val = laplace_op(XBSpline,j0,k_row,XWavelet,j_col,k_col);
                if (fabs(laplace_val)>1e-13 && j_col-j0>=1) {
                    cout << "laplace: (" << j0 << ", " << k_row << "), ("
                         << j_col << ", " << k_col << "): " << laplace_val << endl;
                }
                T convection_val = convection_op(XBSpline,j0,k_row,XWavelet,j_col,k_col);
                if (fabs(convection_val)>1e-13 && j_col-j0>=1) {
                    cout << "convection: (" << j0 << ", " << k_row << "), ("
                         << j_col << ", " << k_col << "): " << convection_val << endl;
                }
            }
        }
    }
    for (int j_row=j0; j_row<=J-1; ++j_row) {
        for (int k_row=basis.rangeJ(j_row).firstIndex(); k_row<=basis.rangeJ(j_row).lastIndex(); ++k_row) {

            T val_p0 = sparsemw_p0(j_row,k_row,XWavelet,0);
            T val_p1 = sparsemw_p1(j_row,k_row,XWavelet,0);
            T val_p2 = sparsemw_p2(j_row,k_row,XWavelet,0);
            T val_p3 = sparsemw_p3(j_row,k_row,XWavelet,0);

            if (fabs(val_p0)>1e-16 && k_row>basis.rangeJL(j_row).lastIndex()
                                   && k_row<basis.rangeJR(j_row).firstIndex()) {
                cout << "p0: (" << j_row << ", " << k_row << "): " << val_p0 << endl;
            }
            if (fabs(val_p1)>1e-16 && k_row>basis.rangeJL(j_row).lastIndex()
                                   && k_row<basis.rangeJR(j_row).firstIndex()) {
                cout << "p1: (" << j_row << ", " << k_row << "): " << val_p1 << endl;
            }
            if (fabs(val_p2)>1e-16 && k_row>basis.rangeJL(j_row).lastIndex()
                                   && k_row<basis.rangeJR(j_row).firstIndex()) {
                cout << "p2: (" << j_row << ", " << k_row << "): " << val_p2 << endl;
            }
            if (fabs(val_p3)>1e-16 && k_row>basis.rangeJL(j_row).lastIndex()
                                   && k_row<basis.rangeJR(j_row).firstIndex()) {
                cout << "p3: (" << j_row << ", " << k_row << "): " << val_p3 << endl;
            }

            for (int k_col=basis.mra.rangeI(j0).firstIndex(); k_col<=basis.mra.rangeI(j0).lastIndex(); ++k_col) {
                T identity_val = identity_op(XWavelet,j_row,k_row,XBSpline,j0,k_col);
                if (fabs(identity_val)>1e-13 && j_row-j0>=1) {
                    cout << "identity: (" << j_row << ", " << k_row << "), ("
                         << j0 << ", " << k_col << "): " << identity_val << endl;
                }
            }
            for (int j_col=j0; j_col<=J-1; ++j_col) {
                for (int k_col=basis.rangeJ(j_col).firstIndex(); k_col<=basis.rangeJ(j_col).lastIndex(); ++k_col) {
                    T identity_val = identity_op(XWavelet,j_row,k_row,XWavelet,j_col,k_col);
                    if (fabs(identity_val)>1e-13 && fabs(j_row-j_col)>1) {
                        cout << "identity: (" << j_row << ", " << k_row << "), ("
                             << j_col << ", " << k_col << "): " << identity_val << endl;
                    }
                    T laplace_val = laplace_op(XWavelet,j_row,k_row,XWavelet,j_col,k_col);
                    if (fabs(laplace_val)>1e-13 && fabs(j_row-j_col)>1) {
                        cout << "laplace: (" << j_row << ", " << k_row << "), ("
                             << j_col << ", " << k_col << "): " << laplace_val << endl;
                    }
                    T convection_val = convection_op(XWavelet,j_row,k_row,XWavelet,j_col,k_col);
                    if (fabs(convection_val)>1e-13 && fabs(j_row-j_col)>1) {
                        cout << "convection: (" << j_row << ", " << k_row << "), ("
                             << j_col << ", " << k_col << "): " << convection_val << endl;
                    }
                }
            }
        }
    }

    return 0;
}
