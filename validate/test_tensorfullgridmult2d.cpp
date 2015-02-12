#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

//  Typedefs for Flens data types:
typedef double T;
typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

//  Typedefs for Basis:
typedef Basis<T, Primal, Interval, Dijkema>                         Basis1D;
typedef TensorBasis2D<Uniform, Basis1D,Basis1D>                     Basis2D;

//  Typedefs for Operators
typedef IdentityOperator1D<T, Basis1D>                              IdentityOp;
typedef PDEOperator1D<T, Basis1D>                                   PDEOp;

typedef UniformTensorMatrix2D<T,Basis2D,PDEOp,IdentityOp,
                              IdentityOp,PDEOp>                     TensorMatrix2D;

typedef IndexSet<Index2D>::const_iterator                           set_it;


int main (int argc, char *argv[]) {
    if (argc!=6) {
        cout << "usage: " << argv[0] << " d d_ j0 Jx Jy" << endl;
        exit(1);
    }

    int d =atoi(argv[1]);
    int d_=atoi(argv[2]);
    int j0=atoi(argv[3]);
    int Jx=atoi(argv[4]);
    int Jy=atoi(argv[5]);

    Basis1D basis(d, d_, j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis,basis);
    IdentityOp  identity_op(basis);
    PDEOp       pde_op(basis,0.5,0.,1.);
    Timer time;

    cout << "Started assembling of uniform tensor matrix." << endl;
    TensorMatrix2D fg_A(basis2d,pde_op,identity_op,identity_op,pde_op,Jx,Jy);
    cout << "Finished assembling of uniform tensor matrix." << endl;


    IndexSet<Index2D> Lambda = fg_A.getIndexSet();
    SparseMatrixT A(Lambda.size(), Lambda.size());
    DenseVectorT  v(Lambda.size());


    cout << "Started assembling of 2d-stiffness matrix." << endl;
    int row_count=1, col_count=1;
    for (set_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        col_count=1;
        v(row_count) = row_count;
        for (set_it col=Lambda.begin(); col!=Lambda.end(); ++col) {
            T val1_x = pde_op((*row).index1,(*col).index1);
            T val1_y = identity_op((*row).index2,(*col).index2);

            T val2_x = identity_op((*row).index1,(*col).index1);
            T val2_y = pde_op((*row).index2,(*col).index2);

            if (fabs(val1_x*val1_y+val2_x*val2_y)>0) {
                A(row_count,col_count) = val1_x*val1_y+val2_x*val2_y;
            }
            ++col_count;
        }
        ++row_count;
    }
    A.finalize();
    cout << "Finished assembling of 2d-stiffness matrix." << endl;

    cout << "Starting sparse matrix vector." << endl;
    time.start();
    DenseVectorT Av1 = A*v;
    time.stop();
    cout << "Elapsed time sparse matrix vector: " << time.elapsed() << endl << endl;

    cout << "Starting uniform tensor matrix vector." << endl;
    time.start();
    DenseVectorT Av2 = fg_A*v;
    time.stop();
    cout << "Elapsed time uniform tensor matrix vector: " << time.elapsed() << endl << endl;

    for (int i=1; i<=(int)Lambda.size(); ++i) {
        if (fabs(Av1(i)-Av2(i))>1e-8) {
            cout << "Position i=" << i << " : " << Av1(i) << " " << Av2(i)
                 << " " << fabs(Av1(i)-Av2(i)) << endl;
        }
    }

    return 0;
}

