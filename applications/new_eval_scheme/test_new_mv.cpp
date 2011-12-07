#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/loc_single_scale_transforms.h>
#include <applications/new_eval_scheme/new_eval.h>
#include <applications/new_eval_scheme/treecoefficients1d.h>
#include <applications/new_eval_scheme/localoperator.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;


typedef IntegralF<Gauss, PrimalBasis>                               IntegralFPrimal;
typedef Integral<Gauss, PrimalBasis, PrimalBasis>                   Integral_Psi_Psi;

void
computeCoefficientVector(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u);

void
computeMatrixVectorRef(const PrimalBasis &test_basis, const PrimalBasis &trial_basis,
                       const Coefficients<Lexicographical,T,Index1D> &v_coeff,
                       const IndexSet<Index1D> &LambdaRowTree,
                       Coefficients<Lexicographical,T,Index1D> &Mv_coeff);

T
u(T x) {
    return sin(2*M_PI*x);
}

DenseVectorT u_singPts;

int main (int argc, char *argv[]) {

    cout.precision(6);
    if (argc!=9) {
        cout << "Usage: " << argv[0] << " test_d test_d_ test_j0 trial_d trial_d_ trial_j0 J shift" << endl;
        return 0;
    }
    int test_d   = atoi(argv[1]);
    int test_d_  = atoi(argv[2]);
    int test_j0  = atoi(argv[3]);
    int trial_d  = atoi(argv[4]);
    int trial_d_ = atoi(argv[5]);
    int trial_j0 = atoi(argv[6]);
    int J  = atoi(argv[7]);
    int shift = atoi(argv[8]);

    if (test_j0!=trial_j0) {
        cout << "Different minimal levels not implemented yet." << endl;
        return 0;
    }

    Timer time;

    PrimalBasis test_basis(test_d,test_d_,test_j0);
    PrimalBasis trial_basis(trial_d,trial_d_,trial_j0);

    LocalOperator<PrimalBasis,PrimalBasis> localoperator(test_basis, false, trial_basis, false, 2);

    IndexSet<Index1D> test_LambdaTree;
    constructRandomGradedTree(test_basis, J+shift, test_LambdaTree);
    //cout << "test_LambdaTree = " << test_LambdaTree << endl;

    IndexSet<Index1D> trial_LambdaTree;
    constructRandomGradedTree(trial_basis, J, trial_LambdaTree);
    //cout << "trial_LambdaTree = " << trial_LambdaTree << endl;

    Coefficients<Lexicographical,T,Index1D> trial_u_multi, test_Mu;
    computeCoefficientVector(trial_basis, trial_LambdaTree, trial_u_multi);

    time.start();
    computeMatrixVectorRef(test_basis, trial_basis, trial_u_multi, test_LambdaTree, test_Mu);
    time.stop();
    T time_computeMatrixVectorRef = time.elapsed();

    CoefficientsByLevel<T> d;
    TreeCoefficients1D<T>  c;
    CoefficientsByLevel<T> PhiPiCheck_vs_v;
    TreeCoefficients1D<T>  PsiLambdaCheck_vs_v;

    for (const_coeff1d_it it=trial_u_multi.begin(); it!=trial_u_multi.end(); ++it) {
        int j = (*it).first.j;
        if ((*it).first.xtype == XBSpline) {
            d[(*it).first.k] = (*it).second;
        }
    }
    c = trial_u_multi;

    for (const_set1d_it it=test_LambdaTree.begin(); it!=test_LambdaTree.end(); ++it) {
        int j = (*it).j;
        if ((*it).xtype == XBSpline) {
            PhiPiCheck_vs_v[(*it).k] = 0.;
            continue;
        }
        else {
            PsiLambdaCheck_vs_v[j].operator[]((*it).k) = 0.;
        }
    }

    time.start();
    localoperator.evalA(test_basis.j0, d, c, PhiPiCheck_vs_v, PsiLambdaCheck_vs_v);
    time.stop();
    PsiLambdaCheck_vs_v[test_basis.j0-1] = PhiPiCheck_vs_v;
    cout << "Sizes of trees: " << test_LambdaTree.size()+trial_LambdaTree.size()
         << " " << time.elapsed() << endl;

    PsiLambdaCheck_vs_v -= test_Mu;
    cout << "Error norm = " << PsiLambdaCheck_vs_v.norm(2.) << endl;
    //cout << "PsiLambdaCheck_vs_v = " << PsiLambdaCheck_vs_v << endl;

//    cout << "PhiPiCheck_vs_v     = " << PhiPiCheck_vs_v << endl;
//    cout << "PsiLambdaCheck_vs_v = " << PsiLambdaCheck_vs_v << endl;
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

void
computeMatrixVectorRef(const PrimalBasis &test_basis, const PrimalBasis &trial_basis,
                       const Coefficients<Lexicographical,T,Index1D> &v_coeff,
                       const IndexSet<Index1D> &LambdaRowTree,
                       Coefficients<Lexicographical,T,Index1D> &Mv_coeff)
{
    Integral_Psi_Psi integral(test_basis,trial_basis);
    for (const_set1d_it row=LambdaRowTree.begin(); row!=LambdaRowTree.end(); ++row) {
        int  j_row      = (*row).j;
        long k_row      = (*row).k;
        XType xtype_row = (*row).xtype;
        T val = 0.;
        for (const_coeff1d_it col=v_coeff.begin(); col!=v_coeff.end(); ++col) {
            int  j_col      = (*col).first.j;
            long k_col      = (*col).first.k;
            XType xtype_col = (*col).first.xtype;
            val += integral(j_row,k_row,xtype_row,0, j_col,k_col,xtype_col,0) * (*col).second;
        }
        Mv_coeff[(*row)] = val;
    }
}

