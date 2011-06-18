#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

//  Typedefs for Flens data types:
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

//  Typedefs for Basis:
typedef Basis<T, Primal, Interval, Dijkema>                         Basis1D;
typedef TensorBasis2D<Uniform, Basis1D,Basis1D>                     Basis2D;

typedef IdentityOperator1D<T, Basis1D>                              IdentityOp1D;
typedef HelmholtzOperator1D<T, Basis1D>                             HelmholtzOp1D;

typedef TensorSparseGrid2D<T, Basis2D, HelmholtzOp1D, IdentityOp1D,
                            IdentityOp1D, HelmholtzOp1D>             TensorSparseGrid;

typedef IndexSet<Index2D>::const_iterator                           set_it;

const T c=1.;

int main (int argc, char *argv[]) {
    if (argc!=5) {
        cout << "usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(1);
    }
    cout.precision(8);
    int d =atoi(argv[1]);
    int d_=atoi(argv[2]);
    int j0=atoi(argv[3]);
    int J =atoi(argv[4]);

    Basis1D basis(d, d_, j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis,basis);

    IdentityOp1D   identity_op_x(basis);
    HelmholtzOp1D  helmholtz_op_x(basis,0.5*c);
    IdentityOp1D   identity_op_y(basis);
    HelmholtzOp1D  helmholtz_op_y(basis,0.5*c);

    TensorSparseGrid sg_A(basis2d, helmholtz_op_x,identity_op_y,
                               identity_op_x, helmholtz_op_y, J, 0.);
    int sg_dim = sg_A.getDimension();
    cout << "Dimension of sparse grid: " << sg_dim << endl;

    DenseVectorT sg_x(sg_dim);
    for (int i=1; i<=sg_dim; ++i) {
        sg_x(i) = i;
    }

    DenseVectorT tmp = sg_A * sg_x;


    IndexSet<Index2D> sg_Lambda = sg_A.getIndexSet();
    Coefficients<Lexicographical,T,Index2D> sg_x_coeffs;
    sg_A.toCoefficients(sg_x, sg_x_coeffs);

    cout << "Started assembling of 2d-stiffness matrix." << endl;
    SparseMatrixT A(sg_Lambda.size(), sg_Lambda.size());
    DenseVectorT  x(sg_Lambda.size());
    int row_count=1, col_count=1;
    for (set_it row=sg_Lambda.begin(); row!=sg_Lambda.end(); ++row) {
        col_count=1;
        x(row_count) = sg_x_coeffs[*row];
        for (set_it col=sg_Lambda.begin(); col!=sg_Lambda.end(); ++col) {
            T val1_x = helmholtz_op_x((*row).index1,(*col).index1);
            T val1_y = identity_op_y((*row).index2,(*col).index2);

            T val2_x = identity_op_x((*row).index1,(*col).index1);
            T val2_y = helmholtz_op_y((*row).index2,(*col).index2);

            if (fabs(val1_x*val1_y+val2_x*val2_y)>0) {
                A(row_count,col_count) = val1_x*val1_y+val2_x*val2_y;
            }
            ++col_count;
        }
        ++row_count;
    }
    A.finalize();
    cout << "Finished assembling of 2d-stiffness matrix." << endl;
    DenseVectorT Ax(sg_Lambda.size()), Ax2(sg_Lambda.size());
    Ax = A*x;


    Coefficients<Lexicographical,T,Index2D> sg_Ax_coeffs;
    sg_A.toCoefficients(tmp, sg_Ax_coeffs);
    row_count=1;
    for (set_it row=sg_Lambda.begin(); row!=sg_Lambda.end(); ++row) {
        if (fabs(Ax(row_count) - sg_Ax_coeffs[*row])>1e-10) {
            cout << *row << " : " << Ax(row_count) << " " << sg_Ax_coeffs[*row]
                 << " " << fabs(Ax(row_count) - sg_Ax_coeffs[*row]) << endl;
        }
        ++row_count;
    }



    return 0;
}
