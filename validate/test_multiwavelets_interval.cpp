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

///     HelmholtzOperator in 1D, i.e. for $a(v,u) = \int(v_x \cdot u_x) + c \cdot \int(v \cdot u)$
typedef HelmholtzOperator1D<T, MWBasis>                             HelmholtzOp;


int main()
{
    /// wavelet basis parameters:
    int d = 2;          // d-wavelets
    int j0 = 2;         // minimal level
    int J = 10;         // maximal level
    cout.precision(16);

    /// Basis initialization, using Dirichlet boundary conditions
    MWBasis basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    /// Plot multi-scaling function and multi-wavelets

    ofstream plotfile_scaling("multiscaling.txt");
    for (T x=0.; x<=1.; x+=pow2i<T>(-7)) {
        plotfile_scaling << x;
        for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            plotfile_scaling << " " << basis.generator(XBSpline)(x,j0,k,0);
        }
        plotfile_scaling << endl;
    }
    ofstream plotfile_wavelet("multiwavelet.txt");
    for (T x=0.; x<=1.; x+=pow2i<T>(-7)) {
        plotfile_wavelet << x;
        for (int k=basis.rangeJ(j0).firstIndex(); k<=basis.rangeJ(j0).lastIndex(); ++k) {
            plotfile_wavelet << " " << basis.generator(XWavelet)(x,j0,k,0);
        }
        plotfile_wavelet << endl;
    }


    /// Operator initialization
    IdentityOp   identity_op(basis);
    HelmholtzOp  helmholtz_op(basis, 1.);


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


    return 0;
}
