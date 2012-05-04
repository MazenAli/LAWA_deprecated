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
///     Multiwavelet Basis over an interval
typedef Basis<T,Orthogonal,Interval,Multi>                          MWBasis;

///     IdentityOperator in 1D, i.e. for $a(v,u) = \cdot \int(v \cdot u)$
typedef IdentityOperator1D<T, MWBasis>                              IdentityOp;

T
p0(T x) {   return 1.; }

T
p1(T x) {   return x; }

T
p2(T x) {   return x*x; }

T
p3(T x) {   return x*x*x; }

int main()
{
    /// wavelet basis parameters:
    int d = 1;          // d-wavelets
    int j0 = 0;         // minimal level
    int J = 3;         // maximal level
    cout.precision(16);

    /// Basis initialization, using Dirichlet boundary conditions
    MWBasis basis(d,j0);
    if (d>1) basis.enforceBoundaryCondition<DirichletBC>();

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
    ofstream plotfile_scaling("interval_multiscaling.txt");
    for (T x=0.; x<=2.; x+=pow2i<T>(-12)) {
        plotfile_scaling << x;
        for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            plotfile_scaling << " " << basis.generator(XBSpline)(x,j0,k,0);
        }
        plotfile_scaling << endl;
    }
    ofstream plotfile_wavelet("interval_multiwavelet.txt");
    for (T x=0.; x<=1.; x+=pow2i<T>(-12)) {
        plotfile_wavelet << x;
        for (int j=j0; j<=J-1; ++j) {
            for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
                plotfile_wavelet << " " << basis.generator(XWavelet)(x,j,k,0);
            }
        }
        plotfile_wavelet << endl;
    }


    /// Operator and Integral initialization
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


    /// Assembler: Check for orthogonality
    Assembler1D<T, MWBasis> assembler(basis);
    SparseMatrixT   identity_A = assembler.assembleStiffnessMatrix(identity_op, J);
    DenseMatrixT identity_A_dense;
    densify(cxxblas::NoTrans,identity_A,identity_A_dense);
    for (int i=1; i<=identity_A_dense.numRows(); ++i) {
        for (int j=1; j<=identity_A_dense.numCols(); ++j) {
            if (i==j) {
                if (fabs(identity_A_dense(i,j)-1.)>1e-12) {
                    cout << "(" << i << ", " << j << "): " << fabs(identity_A_dense(i,j)) << endl;
                }
            }
            else {
                if (fabs(identity_A_dense(i,j))>1e-12) {
                    cout << "(" << i << ", " << j << "): " << fabs(identity_A_dense(i,j)) << endl;
                }
            }
        }

    }
    cout << identity_A_dense << endl;

    /// Assembler: Check for vanishing moments
    for (int j=j0; j<=J-1; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            T val_p0 = mw_p0(j,k,XWavelet,0);
            T val_p1 = mw_p1(j,k,XWavelet,0);
            T val_p2 = mw_p2(j,k,XWavelet,0);
            T val_p3 = mw_p3(j,k,XWavelet,0);

            if (d>=1 && fabs(val_p0)>1e-16 && k>basis.rangeJL(j).lastIndex()
                                   && k<basis.rangeJR(j).firstIndex()) {
                cout << "p0: (" << j << ", " << k << "): " << val_p0 << endl;
            }
            if (d>=2 && fabs(val_p1)>1e-16 && k>basis.rangeJL(j).lastIndex()
                                   && k<basis.rangeJR(j).firstIndex()) {
                cout << "p1: (" << j << ", " << k << "): " << val_p1 << endl;
            }
            if (d>=3 && fabs(val_p2)>1e-16 && k>basis.rangeJL(j).lastIndex()
                                   && k<basis.rangeJR(j).firstIndex()) {
                cout << "p2: (" << j << ", " << k << "): " << val_p2 << endl;
            }
            if (d>=4 && fabs(val_p3)>1e-16 && k>basis.rangeJL(j).lastIndex()
                                   && k<basis.rangeJR(j).firstIndex()) {
                cout << "p3: (" << j << ", " << k << "): " << val_p3 << endl;
            }
        }
    }

    return 0;
}
