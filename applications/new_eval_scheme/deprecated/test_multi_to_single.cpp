#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/source/loc_single_scale_transforms.h>
#include <lawa/methods/adaptive/algorithms/localrefinement.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef Basis<T,Dual,Interval,Dijkema>                              DualBasis;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;



typedef IntegralF<Gauss, PrimalBasis>                               IntegralFPrimal;


void
computeCoefficientVector(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u);

T
u(T x) {
    return sin(2*M_PI*x);
}

DenseVectorT u_singPts;

int main (int argc, char *argv[]) {

    if (argc!=6) {
        cout << "Usage: " << argv[0] << " d d_ j0 J withDirichletBC" << endl;
        return 0;
    }
    int d                = atoi(argv[1]);
    int d_               = atoi(argv[2]);
    int j0               = atoi(argv[3]);
    int J                = atoi(argv[4]);
    bool withDirichletBC = atoi(argv[5]);

    PrimalBasis basis(d,d_,j0);
    DualBasis   dual_basis(d,d_,j0);

    if (withDirichletBC) {
        basis.enforceBoundaryCondition<DirichletBC>();
        dual_basis.enforceBoundaryCondition<DirichletBC>();
    }

    IndexSet<Index1D> LambdaTree;
    constructRandomGradedTree(basis, J, LambdaTree);

    Coefficients<Lexicographical,T,Index1D> u_multi, u_loc_single, u_multi2, diff_u_multi;

    stringstream filename1, filename2;
    filename1 << "u_multi.dat";
    filename2 << "u_multi_coeff";
    computeCoefficientVector(basis, LambdaTree, u_multi);
    plot(basis, u_multi, filename1);
    plotCoeff<T,PrimalBasis>(u_multi, basis, filename2.str().c_str(), false);

    stringstream filename3, filename4;
    filename3 << "u_loc_single.dat";
    filename4 << "u_loc_single_coeff";
    computeMultiToLocallySingleRepr(basis, u_multi, u_loc_single);
    plot(basis, u_loc_single, filename3);
    plotCoeff<T,PrimalBasis>(u_loc_single, basis, filename4.str().c_str(), true);

    stringstream filename5, filename6;
    filename5 << "u_multi2.dat";
    filename6 << "u_multi2_coeff.dat";
    computeLocallySingleToMultiRepr(dual_basis, u_loc_single, LambdaTree, u_multi2);
    plot(basis, u_multi2, filename5);
    plotCoeff<T,PrimalBasis>(u_multi2, basis, filename6.str().c_str(), false);

    diff_u_multi = u_multi-u_multi2;
    cout << "|| u_multi - u_multi2 ||_2 = " << diff_u_multi.norm(2.) << endl;


    return 0;
}

void
computeCoefficientVector(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u_coeff) {

    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;

    typedef IdentityOperator1D<T, PrimalBasis>                          IdentityOp1D;
    typedef NoCompression<T, Index1D, PrimalBasis>                      Compression1D;
    typedef NoPreconditioner<T,Index1D>                                 NoPreconditioner1D;
    typedef MapMatrix<T, Index1D, IdentityOp1D,
                         Compression1D, NoPreconditioner1D>             DataIdentity1D;



    Compression1D compression1d(basis);
    IdentityOp1D identity_op1d(basis);
    NoPreconditioner1D prec1d;
    DataIdentity1D identity_data1d(identity_op1d, prec1d, compression1d);

    SparseMatrixT A_flens(LambdaTree.size(),LambdaTree.size());
    toFlensSparseMatrix(identity_data1d, LambdaTree, LambdaTree, A_flens);

    Function<T> u_Fct(u,u_singPts);
    IntegralFPrimal integral_u_psi(u_Fct, basis);
    DenseVectorT f(LambdaTree.size()), x(LambdaTree.size());
    Coefficients<Lexicographical,T,Index1D> f_coeff;
    int count=1;
    for (const_set1d_it it=LambdaTree.begin(); it!=LambdaTree.end(); ++it, ++count) {
        x(count) = 0.;
        f(count) = integral_u_psi((*it).j,(*it).k,(*it).xtype,0);
        f_coeff[*it] = integral_u_psi((*it).j,(*it).k,(*it).xtype,0);
    }
    int numOfIterations = cg(A_flens, x, f, 1e-16,100);

    count=1;
    for (const_set1d_it it=LambdaTree.begin(); it!=LambdaTree.end(); ++it, ++count) {
        u_coeff[*it] = x(count);
    }

}
