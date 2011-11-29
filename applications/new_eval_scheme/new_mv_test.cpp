#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/loc_single_scale_transforms.h>
#include <applications/new_eval_scheme/new_eval.h>

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

    PrimalBasis test_basis(test_d,test_d_,test_j0);
    PrimalBasis trial_basis(trial_d,trial_d_,trial_j0);

    IndexSet<Index1D> test_LambdaTree;
    constructRandomGradedTree(test_basis, J+shift, test_LambdaTree);
    cout << "test_LambdaTree = " << test_LambdaTree << endl;

    IndexSet<Index1D> trial_LambdaTree;
    constructRandomGradedTree(trial_basis, J, trial_LambdaTree);
    cout << "trial_LambdaTree = " << trial_LambdaTree << endl;

    Coefficients<Lexicographical,T,Index1D> trial_u_multi, test_Mu;
    computeCoefficientVector(trial_basis, trial_LambdaTree, trial_u_multi);
    computeMatrixVectorRef(test_basis, trial_basis, trial_u_multi, test_LambdaTree, test_Mu);


    Coefficients<Lexicographical,T,Index1D>                    d_lM1;
    std::map<int, Coefficients<Lexicographical,T,Index1D> >    c;
    IndexSet<Index1D>                                          SquareCap_lM1;
    std::map<int, IndexSet<Index1D> >                          Lhd;

    for (const_coeff1d_it it=trial_u_multi.begin(); it!=trial_u_multi.end(); ++it) {
        int j = (*it).first.j;
        if ((*it).first.xtype == XBSpline) {
            d_lM1[(*it).first] = (*it).second;
            continue;
        }

        if (c.count(j)==0) {
            Coefficients<Lexicographical,T,Index1D> c_multi_j;
            c_multi_j[(*it).first] = (*it).second;
            c[j] = c_multi_j;
        }
        else {
            c[j].operator[]((*it).first) = (*it).second;
        }
    }

    for (const_set1d_it it=test_LambdaTree.begin(); it!=test_LambdaTree.end(); ++it) {
        int j = (*it).j;
        if ((*it).xtype == XBSpline) {
            SquareCap_lM1.insert((*it));
            continue;
        }

        if (Lhd.count(j)==0) {
            IndexSet<Index1D> Lhd_j;
            Lhd_j.insert((*it));
            Lhd[j] = Lhd_j;
        }
        else {
            Lhd[j].insert((*it));
        }
    }


    cout << "d_lM1 = " << d_lM1 << endl;
    for (int j=trial_basis.j0; j<=J; ++j) {
        cout << "c_" << j << " = " << c[j] << endl;
    }
    cout << "SquareCap_lM1 = " << SquareCap_lM1 << endl;
    for (int j=test_basis.j0; j<=J+shift; ++j) {
        cout << "Lhd_" << j << " = " << Lhd[j] << endl;
    }
    Coefficients<Lexicographical,T,Index1D> upsilon_vs_v, theta_vs_v, test_Mu2, diff;
    eval<T, PrimalBasis, PrimalBasis>(test_basis.j0, test_basis, trial_basis, c, d_lM1,
                                      Lhd, SquareCap_lM1, upsilon_vs_v, theta_vs_v);

    test_Mu2 = upsilon_vs_v + theta_vs_v;
    diff = test_Mu - test_Mu2;
    //cout << "Mu  = " << test_Mu << endl;
    //cout << "Mu2 = " << test_Mu2 << endl;
    cout << "||Mu-Mu2||_2 = " << diff.norm(2.) << endl;

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

