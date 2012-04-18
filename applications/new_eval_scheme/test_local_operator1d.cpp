/* TEST LOCAL OPERATOR
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/source/localoperator1d.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

///  Typedefs for problem components:

///  Wavelet basis over an interval
//typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;
typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

///  Underlying bilinear form
typedef HelmholtzOperator1D<T,PrimalBasis>                          BilinearForm;
typedef HelmholtzOperator1D<T,RefinementBasis>                      RefinementBilinearForm;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm>                     LocOp1D;

typedef CoefficientsByLevel<T>::const_it                            const_coeffbylevel_it;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;

void
constructRandomTree(const PrimalBasis &basis, int J, bool withRandomValues,
                    TreeCoefficients1D<T> &LambdaTree, int seed);

void
computeEvalLRef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Lv_tree);

void
computeEvalURef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Uv_tree);

void
computeEvalARef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Av_tree);

int main(int argc, char*argv[])
{
    cout.precision(20);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d j0 J seed" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d    = atoi(argv[1]);
    int j0   = atoi(argv[2]);
    int J    = atoi(argv[3]);
    int seed = atoi(argv[4]);


    Timer time;

    /// Basis initialization, using Dirichlet boundary conditions
    //PrimalBasis basis(d, j0);           // For L2_orthonormal and special MW bases
    PrimalBasis basis(d, d, j0);      // For biorthogonal wavelet bases
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis &refinementbasis = basis.refinementbasis;

    BilinearForm              Bil(basis,1.);
    RefinementBilinearForm    RefinementBil(refinementbasis,1.);

    LocOp1D localOperator1D(basis,basis,RefinementBil);

    TreeCoefficients1D<T> v_tree(COEFFBYLEVELSIZE), Av_tree(COEFFBYLEVELSIZE),
                          Uv_tree(COEFFBYLEVELSIZE), Lv_tree(COEFFBYLEVELSIZE);
    TreeCoefficients1D<T> Av_ref_tree(COEFFBYLEVELSIZE), Uv_ref_tree(COEFFBYLEVELSIZE),
                          Lv_ref_tree(COEFFBYLEVELSIZE);
    constructRandomTree(basis, J, true, v_tree, seed);
    constructRandomTree(basis, J+1, false, Av_tree, seed+37);
    Av_ref_tree = Av_tree;
    Uv_tree     = Av_tree;
    Uv_ref_tree = Av_tree;
    Lv_tree     = Av_tree;
    Lv_ref_tree = Av_tree;



    cout << " ********* Testing evalA **********" << endl;
    cout << "N = " << Av_tree.size() + v_tree.size() << endl;
    time.start();
    computeEvalARef(Bil, basis, v_tree, Av_ref_tree);
    time.stop();
    cout << "Elapsed time for EvalARef: " << time.elapsed() << endl;
    time.start();
    localOperator1D.eval(v_tree, Av_tree,"A");
    time.stop();
    cout << "Elapsed time for evalA: " << time.elapsed() << endl;
    //cout << "v_tree = " << v_tree << endl;
    //cout << "Av_tree_ref = " << Av_ref_tree << endl;
    //cout << "Av_tree = " << Av_tree << endl;
    Av_ref_tree -= Av_tree;
    cout << "Error norm: " << Av_ref_tree.norm(2.) << endl;
    cout << " **********************************" << endl << endl;


    cout << " ********* Testing evalU **********" << endl;
    cout << "N = " << Uv_tree.size() + v_tree.size() << endl;
    time.start();
    computeEvalURef(Bil, basis, v_tree, Uv_ref_tree);
    time.stop();
    time.start();
    localOperator1D.eval(v_tree, Uv_tree,"U");
    time.stop();
    cout << "Elapsed time for evalU: " << time.elapsed() << endl;
    //cout << "v_tree = " << v_tree << endl;
    //cout << "Uv_tree_ref = " << Uv_ref_tree << endl;
    //cout << "Uv_tree = " << Uv_tree << endl;
    Uv_ref_tree -= Uv_tree;
    cout << "Error norm: " << Uv_ref_tree.norm(2.) << endl;
    cout << " **********************************" << endl << endl;

    cout << " ********* Testing evalL **********" << endl;
    cout << "N = " << Lv_tree.size() + v_tree.size() << endl;
    time.start();
    computeEvalLRef(Bil, basis, v_tree, Lv_ref_tree);
    time.stop();
    time.start();
    localOperator1D.eval(v_tree, Lv_tree,"L");
    time.stop();
    cout << "Elapsed time for evalL: " << time.elapsed() << endl;
    //cout << "v_tree = " << v_tree << endl;
    //cout << "Uv_tree_ref = " << Uv_ref_tree << endl;
    //cout << "Uv_tree = " << Uv_tree << endl;
    Lv_ref_tree -= Lv_tree;
    cout << "Error norm: " << Lv_ref_tree.norm(2.) << endl;
    cout << " **********************************" << endl << endl;

    return 0;
}

void
constructRandomTree(const PrimalBasis &basis, int J, bool withRandomValues,
                    TreeCoefficients1D<T> &LambdaTree, int seed)
{
    srand ( seed );
    T val = 0.;
    /*
    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        val = withRandomValues ? (rand() / (T)RAND_MAX) : 0.;
        LambdaTree[0].map.operator[](k) = val;
    }
    for (int j=basis.j0; j<=J; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            val = withRandomValues ? (rand() / (T)RAND_MAX) : 0.;
            LambdaTree.bylevel[j-basis.j0+1].map.operator[](k) = val;
        }
    }
    */
    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
        LambdaTree[0].map.operator[](k) = val;
    }
    for (int j=basis.j0; j<=J; ++j) {
        val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
        int random_k1 = rand() % basis.cardJ(j) + 1;
        LambdaTree[j-basis.j0+1].map.operator[](random_k1) = val;
        val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
        int random_k2 = rand() % basis.cardJ(j) + 1;
        LambdaTree[j-basis.j0+1].map.operator[](random_k2) = val;
    }
    for (int i=J-basis.j0+1; i>=2; --i) {
        CoefficientsByLevel<T> *currentlevel;
        currentlevel = &(LambdaTree.bylevel[i]);
        int j=basis.j0+i-1;
        for (const_coeffbylevel_it it=(*currentlevel).map.begin(); it!=(*currentlevel).map.end(); ++it) {
            long k = (*it).first;
            long k_first = std::max((int)k / 2 - 30, basis.rangeJ(j-1).firstIndex());
            long k_last  = std::min((int)k / 2 + 30, basis.rangeJ(j-1).lastIndex());
            for (int k1=k_first; k1<=k_last; ++k1) {
                if (overlap(basis.psi.support(j,k),basis.psi.support(j-1,k1))>0) {
                    val = withRandomValues ? ((T)rand() / RAND_MAX) : 0.;
                    LambdaTree[i-1].map.operator[](k1) = val;
                }
            }
        }
    }
}

void
computeEvalARef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Av_tree)
{
    Coefficients<Lexicographical,T,Index1D> v, Av;
    fromTreeCofficientsToCofficients(basis, v_tree, v);
    fromTreeCofficientsToCofficients(basis, Av_tree, Av);

    for (coeff1d_it row=Av.begin(); row!=Av.end(); ++row) {
        long double val = 0.L;
        for (const_coeff1d_it col=v.begin(); col!=v.end(); ++col) {
            val +=  Bil((*row).first,(*col).first) * (*col).second;
        }
        (*row).second = val;
    }
    fromCofficientsToTreeCofficients(basis, Av, Av_tree);
}

void
computeEvalURef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Uv_tree)
{
    Coefficients<Lexicographical,T,Index1D> v, Uv;
    fromTreeCofficientsToCofficients(basis, v_tree, v);
    fromTreeCofficientsToCofficients(basis, Uv_tree, Uv);

    for (coeff1d_it row=Uv.begin(); row!=Uv.end(); ++row) {
        long double val = 0.L;
        for (const_coeff1d_it col=v.begin(); col!=v.end(); ++col) {
            if ( ((*row).first.xtype==XBSpline) ||
                 (((*row).first.xtype==XWavelet && (*col).first.xtype==XWavelet
                                          && (*row).first.j<=(*col).first.j)) ) {
                val +=  Bil((*row).first,(*col).first) * (*col).second;
            }
        }
        (*row).second = val;
    }
    fromCofficientsToTreeCofficients(basis, Uv, Uv_tree);
}

void
computeEvalLRef(const BilinearForm &Bil, const PrimalBasis &basis,
                const TreeCoefficients1D<T> &v_tree, TreeCoefficients1D<T> &Lv_tree)
{
    Coefficients<Lexicographical,T,Index1D> v, Lv;
    fromTreeCofficientsToCofficients(basis, v_tree, v);
    fromTreeCofficientsToCofficients(basis, Lv_tree, Lv);

    for (coeff1d_it row=Lv.begin(); row!=Lv.end(); ++row) {
        long double val = 0.L;
        for (const_coeff1d_it col=v.begin(); col!=v.end(); ++col) {
            if ( !( (   (*row).first.xtype==XBSpline) ||
                    ( ( (*row).first.xtype==XWavelet && (*col).first.xtype==XWavelet
                                          && (*row).first.j<=(*col).first.j) ) ) ) {
                val +=  Bil((*row).first,(*col).first) * (*col).second;
            }
        }
        (*row).second = val;
    }
    fromCofficientsToTreeCofficients(basis, Lv, Lv_tree);
}
