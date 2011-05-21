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

//  Typedefs for Operators
typedef IdentityOperator1D<T, Basis1D>                              IdentityOp;
typedef PDEOperator1D<T, Basis1D>                                   PDEOp;

typedef UniformTensorMatrix2D<T,Basis2D,PDEOp,IdentityOp,
                              IdentityOp,PDEOp>                     TensorMatrix2D;

//  Typedefs for Compression
typedef CompressionPDE2D<T, Basis2D>                                Compression;

typedef IndexSet<Index2D>::const_iterator                           set_it;

IndexSet<Index2D>
getFullTensorIndexSet(const Basis1D &basis, int Jx, int Jy);

void
mv_sparse_dense(const SparseMatrixT &M, const DenseMatrixT &V, DenseMatrixT &res);

void
mv_dense_sparse(const DenseMatrixT &V, const SparseMatrixT &M, DenseMatrixT &res);

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

    Compression Compr(basis2d,false,20);

    IdentityOp  identity_op(basis);
    PDEOp       pde_op(basis,2.,3.,1.);

    IndexSet<Index2D> Lambda = getFullTensorIndexSet(basis, Jx, Jy);

    SparseMatrixT A(Lambda.size(), Lambda.size());
    DenseVectorT  v(Lambda.size());

    Timer time;

    Compr.setParameters(Lambda);

    cout << "Started assembling of 2d-stiffness matrix." << endl;
    int row_count=1, col_count=1;
    for (set_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        col_count=1;
        v(row_count) = 1.;
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

    time.start();
    DenseVectorT Av = A*v;
    time.stop();
    cout << "Elapsed time sparse matrix vector: " << time.elapsed() << endl << endl;
    cout << "Av=" << Av << endl;

    cout << "Started assembling of tensor matrix." << endl;
    TensorMatrix2D TensorA(basis2d,pde_op,identity_op,identity_op,pde_op,Jx,Jy);
    time.start();
    DenseVectorT ret = TensorA*v;
    time.stop();
    cout << "Elapsed time uniform tensor matrix vector (vector): " << time.elapsed() << endl << endl;
    cout << "Av=" << ret << endl;

    return 0;
}

IndexSet<Index2D>
getFullTensorIndexSet(const Basis1D &basis, int Jx, int Jy) {

    IndexSet<Index2D> Lambda;

    for (int k_x=basis.mra.rangeI(basis.j0).firstIndex(); k_x<=basis.mra.rangeI(basis.j0).lastIndex(); ++k_x) {
        for (int k_y=basis.mra.rangeI(basis.j0).firstIndex(); k_y<=basis.mra.rangeI(basis.j0).lastIndex(); ++k_y) {
            Lambda.insert(Index2D(Index1D(basis.j0,k_x,XBSpline),Index1D(basis.j0,k_y,XBSpline)));
        }
        for (int j_y=basis.j0; j_y<=Jy-1; ++j_y) {
            for (int k_y=basis.rangeJ(j_y).firstIndex(); k_y<=basis.rangeJ(j_y).lastIndex(); ++k_y) {
                Lambda.insert(Index2D(Index1D(basis.j0,k_x,XBSpline),Index1D(j_y,k_y,XWavelet)));
            }
        }
    }
    for (int j_x=basis.j0; j_x<=Jx-1; ++j_x) {
        for (int k_x=basis.rangeJ(j_x).firstIndex(); k_x<=basis.rangeJ(j_x).lastIndex(); ++k_x) {
            for (int k_y=basis.mra.rangeI(basis.j0).firstIndex(); k_y<=basis.mra.rangeI(basis.j0).lastIndex(); ++k_y) {
                Lambda.insert(Index2D(Index1D(j_x,k_x,XWavelet),Index1D(basis.j0,k_y,XBSpline)));
            }
            for (int j_y=basis.j0; j_y<=Jy-1; ++j_y) {
                for (int k_y=basis.rangeJ(j_y).firstIndex(); k_y<=basis.rangeJ(j_y).lastIndex(); ++k_y) {
                    Lambda.insert(Index2D(Index1D(j_x,k_x,XWavelet),Index1D(j_y,k_y,XWavelet)));
                }
            }
        }
    }
    return Lambda;
}

void
mv_sparse_dense(const SparseMatrixT &B, const DenseMatrixT &V, DenseMatrixT &res) {
    if (B.numCols()!=V.numRows()) {
        std::cerr << "Dimension of matrices do not fit." << std::endl;
        exit(1);
    }
    res.engine().resize(B.numRows(),V.engine().numCols());
    for (int i=1; i<=V.numRows(); ++i) {
        DenseVectorT tmp = B*V(_,i);
        res(_,i) = tmp;
    }

}

void
mv_dense_sparse(const DenseMatrixT &V, const SparseMatrixT &B, DenseMatrixT &res) {
    if (B.numCols()!=V.numCols()) {
        std::cerr << "Dimension of matrices do not fit." << std::endl;
        exit(1);
    }
    res.engine().resize(B.numRows(),V.engine().numCols());
    for (int i=1; i<=V.numRows(); ++i) {
        DenseVectorT v=V(i,_);
        DenseVectorT tmp = B*v;
        /*
        DenseMatrixT B_dense;
        densify(cxxblas::NoTrans,B,B_dense);
        cout << "B = " << B_dense << endl;
        cout << "V(i,_)=" <<  V(i,_) << endl;
        cout << "B*V(_,i)=" << tmp << endl;
        */
        res(i,_) = tmp;
    }
}

/*
cout << "Started assembling of 1d-stiffness matrices." << endl;
Assembler1D<T, Basis1D> assembler(basis);
SparseMatrixT   A1d = assembler.assembleStiffnessMatrix(pde_op, J);
SparseMatrixT   M1d = assembler.assembleStiffnessMatrix(identity_op, J);
DenseMatrixT    V(basis.mra.cardI(J),basis.mra.cardI(J));

for (int i=1; i<=basis.mra.cardI(J); ++i) {
    for (int j=1; j<=basis.mra.cardI(J); ++j) {
        V(i,j) = 1.;
    }
}
cout << "Finished assembling of 1d-stiffness matrices." << endl;

DenseMatrixT res, res2;
time.start();
mv_sparse_dense(M1d, V, res);

mv_dense_sparse(res, A1d, res2);
time.stop();
cout << "Elapsed time matrix vector (tensor): " << time.elapsed() << endl << endl;
*/
//cout << "Av1=" << res2 << endl;


/*
DenseMatrixT A1d_dense;
densify(cxxblas::NoTrans,A1d,A1d_dense);
DenseMatrixT M1d_dense;
densify(cxxblas::NoTrans,M1d,M1d_dense);
DenseMatrixT A_dense;
densify(cxxblas::NoTrans,A,A_dense);
cout << "A1d = " << A1d_dense << endl;
cout << "M1d = " << M1d_dense << endl;
cout << "A = " << A_dense << endl;
cout << "(IxM1d)v=" << res << endl;

UniformIndex2D<Basis2D> uniformindex(basis2d,J,J);
for (set_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
    cout << *it << ": " << uniformindex((*it).index1.xtype,(*it).index1.j,(*it).index1.k,
                                        (*it).index2.xtype,(*it).index2.j,(*it).index2.k) << endl;
}
*/


