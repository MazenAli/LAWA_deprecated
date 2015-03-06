#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                          DenseVectorT;

int main(){

    FullColMatrixT A(4, 4);

    A = 1,1,1,1,
        1,1,1,1,
        1,1,1,1,
        1,1,1,1;
    cout << "A = " << A << endl;
    cout << endl;

    DenseVectorT u1(4);
    u1 = 1,1,1,1;
    cout << "u1 = " << u1 << endl;

    DenseVectorT Au1(4);
    Au1 = A * u1;
    cout << "A * u1 = " << Au1 << endl;
    T u1Au1 = u1 * Au1;
    cout << "u1 * A * u1 = " << u1Au1 << endl;
    cout << endl << endl;

    cout << "A = " << A << endl;

    /*
     // Assertion fails: (x.length()==((trans==cxxblas::NoTrans) ? A.numCols() : A.numRows())), function mv, file /Users/ksteih/MyLibraries/LAWA/Flens-lite_LAWA/flens/blas/level2/mv.tcc, line 61.
    DenseVectorT u2(2);
    u2 = 2,2;
    cout << "u2 = " << u2 << endl;

    DenseVectorT Au2 = A * u2;
    cout << "A * u2 = " << Au2 << endl;
    T u2Au2 = u2 * Au2;
    cout << "u2 * A * u2 = " << u2Au2 << endl;
     */

    SparseMatrixT A_sparse(6,4);
    A_sparse(1,1) = 1;
    A_sparse(2,4) = 1;
    A_sparse(3,1) = 1;
    A_sparse(3,4) = 1;
    A_sparse(4,4) = 1;
    A_sparse(5,2) = 1;
    A_sparse(5,1) = 1;
    A_sparse(6,3) = 1;
    A_sparse.finalize();

    DenseVectorT u_sparse1(6);
    DenseVectorT u_sparse2(6);

    flens::mv(cxxblas::NoTrans, 1., A_sparse, u1, 0., u_sparse1);
    cout << "(mv) A_sparse * u1 = " << u_sparse1 << endl;
    u_sparse2 = A_sparse * u1;
    cout << "A_sparse * u1 = " << u_sparse2 << endl;

    return 0;
}
