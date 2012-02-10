#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/source/loc_single_scale_transforms.h>
#include <applications/new_eval_scheme/source/new_eval.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>
#include <applications/new_eval_scheme/source/localoperator.h>

using namespace std;
using namespace lawa;

typedef double T;


typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;

//typedef HelmholtzOperator1D<T,PrimalBasis>                          BilinearForm;
typedef IdentityOperator1D<T,PrimalBasis>                          BilinearForm;
//typedef WeightedLaplaceOperator1D<T,Primal,Gauss>                   BilinearForm;

typedef DiagonalMatrixPreconditioner1D<T,PrimalBasis,BilinearForm>  Preconditioner;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;

typedef IntegralF<Gauss, PrimalBasis>                               IntegralFPrimal;
typedef Integral<Gauss, PrimalBasis, PrimalBasis>                   Integral_Psi_Psi;

void
computeCoefficientVector(const IndexSet<Index1D> &LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u);

void
computeEvalARef(const BilinearForm &Bil, const Preconditioner &Prec,
                const Coefficients<Lexicographical,T,Index1D> &v_coeff,
                const IndexSet<Index1D> &LambdaRowTree,
                Coefficients<Lexicographical,T,Index1D> &EvalA_coeff);

void
computeEvalURef(const BilinearForm &Bil, const Preconditioner &Prec,
                const Coefficients<Lexicographical,T,Index1D> &v_coeff,
                const IndexSet<Index1D> &LambdaRowTree,
                Coefficients<Lexicographical,T,Index1D> &EvalL_coeff);

T
u(T x) {
    return sin(2*M_PI*x);
}

DenseVectorT u_singPts;

int main (int argc, char *argv[]) {
    Timer reftime;
    T total_time = 0.;
    reftime.start();
    cout.precision(20);
    if (argc!=6) {
        cout << "Usage: " << argv[0] << " d d_ j0 J shift" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);
    int J  = atoi(argv[4]);
    int shift = atoi(argv[5]);

    Timer time;
    T time_EvalA, time_EvalL, time_EvalU;

    PrimalBasis     trial_basis(d,d_,j0);
    //BilinearForm    Bil(trial_basis,0.);
    BilinearForm    Bil(trial_basis);
    Preconditioner  Prec(Bil);
    int offset=5;
    if (d==2 && d_==2) {
        offset=2;
    }
    LocalOperator<PrimalBasis,PrimalBasis, BilinearForm, Preconditioner> localoperator(trial_basis, false, trial_basis, false, offset, Bil, Prec);

    stringstream ct_filename;
    ct_filename << "comptime_fast_eval1d_" << d << "_" << d_ << "_" << j0 << "_" << J << ".dat";
    ofstream ct_file(ct_filename.str().c_str());

    for (int j=J; j<=J; ++j) {

        cout << "**** j = " << j << " ****" << endl;
        IndexSet<Index1D> test_LambdaTree, trial_LambdaTree;
        Coefficients<Lexicographical,T,Index1D> trial_u_multi;
        /*
        constructRandomGradedTree(trial_basis, j+shift, test_LambdaTree);
        constructRandomGradedTree(trial_basis, j, trial_LambdaTree);
        */

        for (int k=trial_basis.mra.rangeI(j0).firstIndex(); k<=trial_basis.mra.rangeI(j0).lastIndex(); ++k) {
            trial_LambdaTree.insert(Index1D(j0,k,XBSpline));
            test_LambdaTree.insert(Index1D(j0,k,XBSpline));
        }
        for (int j1=j0; j1<=j; ++j1) {
            for (int k=trial_basis.rangeJ(j1).firstIndex(); k<=trial_basis.rangeJ(j1).lastIndex(); ++k) {
                trial_LambdaTree.insert(Index1D(j1,k,XWavelet));
                test_LambdaTree.insert(Index1D(j1,k,XWavelet));
            }
        }

        int N = test_LambdaTree.size()+trial_LambdaTree.size();
        computeCoefficientVector(trial_LambdaTree, trial_u_multi);

        stringstream filename1;
        filename1 << "u_multi_" << j;
        plotCoeff<T,PrimalBasis>(trial_u_multi, trial_basis, filename1.str().c_str(), false);


        size_t hashmap_size = COEFFBYLEVELSIZE;
        CoefficientsByLevel<T> PhiPiCheck_vs_v(j0,hashmap_size), U_PhiPiCheck_vs_v(j0,hashmap_size);
        TreeCoefficients1D<T>  c(hashmap_size);
        TreeCoefficients1D<T>  PsiLambdaCheck_vs_v(hashmap_size),
                               L_PsiLambdaCheck_vs_v(hashmap_size),
                               U_PsiLambdaCheck_vs_v(hashmap_size);
        TreeCoefficients1D<T>  PsiLambdaCheck_vs_v_sca(hashmap_size),
                               L_PsiLambdaCheck_vs_v_sca(hashmap_size),
                               U_PsiLambdaCheck_vs_v_sca(hashmap_size);

        c = trial_u_multi;

        for (const_set1d_it it=test_LambdaTree.begin(); it!=test_LambdaTree.end(); ++it) {
            int j = (*it).j;
            if ((*it).xtype == XBSpline) {
                PhiPiCheck_vs_v.map[(*it).k] = 0.;
                U_PhiPiCheck_vs_v.map[(*it).k] = 0.;
                continue;
            }
            else {
                PsiLambdaCheck_vs_v[j].map[((*it).k)] = 0.;
                L_PsiLambdaCheck_vs_v[j].map[((*it).k)] = 0.;
                U_PsiLambdaCheck_vs_v[j].map[((*it).k)] = 0.;
            }
        }

        Coefficients<Lexicographical,T,Index1D> EvalA_coeff;
        //computeEvalARef(Bil, Prec, trial_u_multi, test_LambdaTree, EvalA_coeff);
        TreeCoefficients1D<T>  c1(hashmap_size);
        CoefficientsByLevel<T> d1(j0,hashmap_size);
        localoperator.scale_wrt_trialbasis(c,c1);
        d1 = c1[j0-1];
        reftime.stop();
        total_time += reftime.elapsed();
        time.start();
        localoperator.evalA(j0, d1, c1, PhiPiCheck_vs_v, PsiLambdaCheck_vs_v);
        time.stop();
        time_EvalA = time.elapsed();
        reftime.start();
        PsiLambdaCheck_vs_v[j0-1] = PhiPiCheck_vs_v;
        localoperator.scale_wrt_trialbasis(PsiLambdaCheck_vs_v,PsiLambdaCheck_vs_v_sca);
        PsiLambdaCheck_vs_v_sca -= EvalA_coeff;
        cout << "EvalA: error = " << PsiLambdaCheck_vs_v_sca.norm(2.) << endl;
        cout << "       " << N << " " << time_EvalA << endl;

        Coefficients<Lexicographical,T,Index1D> EvalU_coeff;
        //computeEvalURef(Bil, Prec, trial_u_multi, test_LambdaTree, EvalU_coeff);
        TreeCoefficients1D<T>  c2(hashmap_size);
        CoefficientsByLevel<T> d2(j0,hashmap_size);
        localoperator.scale_wrt_trialbasis(c,c2);
        d2 = c2[j0-1];
        reftime.stop();
        total_time += reftime.elapsed();
        time.start();
        localoperator.evalU(j0, d2, c2, U_PhiPiCheck_vs_v, U_PsiLambdaCheck_vs_v);
        time.stop();
        reftime.start();
        T time_EvalU = time.elapsed();
        U_PsiLambdaCheck_vs_v[j0-1] = U_PhiPiCheck_vs_v;
        localoperator.scale_wrt_trialbasis(U_PsiLambdaCheck_vs_v,U_PsiLambdaCheck_vs_v_sca);
        U_PsiLambdaCheck_vs_v_sca -= EvalU_coeff;
        cout << "EvalU: error = " << U_PsiLambdaCheck_vs_v_sca.norm(2.) << endl;
        cout << "       " << N << " " << time_EvalU << endl;

        Coefficients<Lexicographical,T,Index1D> EvalL_coeff;
        EvalL_coeff = EvalA_coeff;
        EvalL_coeff-= EvalU_coeff;
        TreeCoefficients1D<T>  c3(hashmap_size);
        CoefficientsByLevel<T> d3(j0,hashmap_size);
        localoperator.scale_wrt_trialbasis(c,c3);
        d3 = c3[j0-1];
        reftime.stop();
        total_time += reftime.elapsed();
        time.start();
        localoperator.evalL(j0, d3, c3, L_PsiLambdaCheck_vs_v);
        time.stop();
        T time_EvalL = time.elapsed();
        reftime.start();
        localoperator.scale_wrt_trialbasis(L_PsiLambdaCheck_vs_v,L_PsiLambdaCheck_vs_v_sca);
        L_PsiLambdaCheck_vs_v_sca -= EvalL_coeff;
        cout << "EvalL: error = " << L_PsiLambdaCheck_vs_v_sca.norm(2.) << endl;
        cout << "       " << N << " " << time_EvalL << endl;

        ct_file << N << " " << time_EvalA << " " << time_EvalL << " "  << time_EvalU << endl;

        cout << endl;
        reftime.stop();
        total_time += reftime.elapsed();
        cout << "  Reference time: " << total_time << endl;
    }
}

void
computeCoefficientVector(const IndexSet<Index1D> &LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u_coeff) {
    for (const_set1d_it it=LambdaTree.begin(); it!=LambdaTree.end(); ++it) {
        u_coeff[*it] = (T)rand() / RAND_MAX;
    }
}


void
computeEvalARef(const BilinearForm &Bil, const Preconditioner &Prec,
                const Coefficients<Lexicographical,T,Index1D> &v_coeff,
                const IndexSet<Index1D> &LambdaRowTree,
                Coefficients<Lexicographical,T,Index1D> &EvalA_coeff)
{
    for (const_set1d_it row=LambdaRowTree.begin(); row!=LambdaRowTree.end(); ++row) {
        long double val = 0.L;
        long double prec_row_val = Prec(*row);
        for (const_coeff1d_it col=v_coeff.begin(); col!=v_coeff.end(); ++col) {
            val +=  (prec_row_val * Bil(*row,(*col).first) * Prec((*col).first) ) * (*col).second;
        }
        EvalA_coeff[(*row)] = (T)val;
    }
}


void
computeEvalURef(const BilinearForm &Bil, const Preconditioner &Prec,
                const Coefficients<Lexicographical,T,Index1D> &v_coeff,
                const IndexSet<Index1D> &LambdaRowTree,
                Coefficients<Lexicographical,T,Index1D> &EvalL_coeff)
{
    for (const_set1d_it row=LambdaRowTree.begin(); row!=LambdaRowTree.end(); ++row) {
        long double val = 0.L;
        long double prec_row_val = Prec(*row);
        for (const_coeff1d_it col=v_coeff.begin(); col!=v_coeff.end(); ++col) {
            if ( ((*row).xtype==XBSpline) ||
                 (((*row).xtype==XWavelet && (*col).first.xtype==XWavelet  && (*row).j<=(*col).first.j)) ) {
                val +=  (prec_row_val * Bil(*row,(*col).first) * Prec((*col).first) ) * (*col).second;
            }
        }
        EvalL_coeff[(*row)] = (T)val;
    }
}

/*
 *     typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;

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
    Integral_Psi_Psi integral(basis,basis);
 */
