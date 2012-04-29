/* TEST LOCAL OPERATOR
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/source/localoperator1d.h>
#include <applications/new_eval_scheme/source/localoperator2d_new.h>
#include <applications/new_eval_scheme/source/multitreeoperations.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

///  Typedefs for problem components:

///  Wavelet basis over an interval
typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;
//typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
bool isL2Orthonormal_x = true;
bool isL2Orthonormal_y = true;

typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

///  Underlying bilinear form
typedef IdentityOperator1D<T,PrimalBasis>                           BilinearForm_x;
typedef RefinementBasis::IdentityOperator1D                         RefinementBilinearForm_x;
typedef IdentityOperator1D<T,PrimalBasis>                           BilinearForm_y;
typedef RefinementBasis::IdentityOperator1D                         RefinementBilinearForm_y;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm_x>                   LocalOp1D_x;
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm_y>                   LocalOp1D_y;
typedef LocalOperator2DNew<LocalOp1D_x, LocalOp1D_y>                LocalOp2D;

///  Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

typedef AlignedCoefficients<T,Index2D,Index1D,Index1D>              alignedCoefficients;

void
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma=0.);

void
getRandomCoefficientVector(const IndexSet<Index2D> &Lambda,
                           Coefficients<Lexicographical,T,Index2D> &coeff);

void
refComputationIAv(const BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &IAv);

void
refComputationLIIAv(const BilinearForm_x &Bil_y, const Coefficients<Lexicographical,T,Index2D> &IAv,
                    Coefficients<Lexicographical,T,Index2D> &LIIAv);

void
refComputationUIv(const BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &UIv);
void
refComputationIAUIv(const BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &UIv,
                    Coefficients<Lexicographical,T,Index2D> &IAUIv);

void
refComputationAAv(const BilinearForm_x &Bil_x, const BilinearForm_y &Bil_y,
                  const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AAv);

int main (int argc, char *argv[]) {

#ifdef TRONE
    cout << "using tr1." << endl;
#else
    cout << "using gnu_cxx." << endl;
#endif
    cout.precision(6);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);

    int numOfIter=J;
    bool withDirichletBC=true;
    bool useSparseGrid=true;
    bool calcRefSol=true;

    /// Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, j0);           // For L2_orthonormal and special MW bases
    //PrimalBasis basis(d, d, j0);      // For biorthogonal wavelet bases
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    BilinearForm_x    Bil_x(basis);
    BilinearForm_y    Bil_y(basis);

    LocalOp1D_x localOperator_x(basis,basis,refinementbasis.IdentityOp1D);
    LocalOp1D_y localOperator_y(basis,basis,refinementbasis.IdentityOp1D);

    LocalOp2D   localop2d(localOperator_x,localOperator_y);

    Timer time;

    ofstream file2("comptimes_mv2d.dat");

    T old_time = 1.;
    T old_N = 1.;
    T time_intermediate1=0., time_intermediate2=0.,
                  time_IAv1=0., time_IAv2=0., time_LIv=0., time_UIv=0.;
    T time_intermediate1_old=0., time_intermediate2_old=0.,
      time_IAv1_old=0., time_IAv2_old=0., time_LIv_old=0., time_UIv_old=0.;
    int N = 0, N_old = 0;

    for (int j=numOfIter; j<=numOfIter; ++j) {

        //localop2d.setJ(J+3);
        IndexSet<Index2D> checkLambda, checkLambda2, Lambda, Lambda2;

        getSparseGridIndexSet(basis,checkLambda,j,0.2);
        getSparseGridIndexSet(basis,checkLambda2,j,0.2);
        getSparseGridIndexSet(basis,Lambda,j,0.2);
        getSparseGridIndexSet(basis,Lambda2,j,0.2);

        cout << "#checkLambda  = " << checkLambda.size() << endl;
        Index1D index1_x(j0+j+2,25,XWavelet);
        Index1D index1_y(j0+j+2,12,XWavelet);
        Index2D new_index1(index1_x,index1_y);
        extendMultiTree( basis2d,new_index1,checkLambda);
        extendMultiTree2(basis2d,new_index1,20,checkLambda2);
        cout<< "#checkLambda1 = " << checkLambda.size() << endl;
        cout<< "#checkLambda2 = " << checkLambda2.size() << endl;


        cout << "#Lambda  = " << Lambda.size() << endl;
        Index1D index2_x(j0+j+2,25,XWavelet);
        Index1D index2_y(j0+j+2,12,XWavelet);
        Index2D new_index2(index2_x,index2_y);
        extendMultiTree( basis2d,new_index2,Lambda);
        extendMultiTree2(basis2d,new_index2,20,Lambda2);
        cout<< "#Lambda1 = " << Lambda.size() << endl;
        cout<< "#Lambda2 = " << Lambda2.size() << endl;

        /*
        IndexSet<Index2D> C_Lambda =  C(Lambda, 1., basis2d);
        for (const_set2d_it it=C_Lambda.begin(); it!=C_Lambda.end(); ++it) {
            std::cerr << "CHECKING " << (*it) << std::endl;
            extendMultiTree2(basis2d,(*it),offset,Lambda);
        }
        std::cerr << "Size of Lambda: " << Lambda.size() << ", size of C(Lambda): " << C_Lambda.size() << std::endl;
        */

        Coefficients<Lexicographical,T,Index2D> v(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index2D> intermediate(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index2D> LIIAv(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index2D> IAUIv(SIZEHASHINDEX2D);

        getRandomCoefficientVector(Lambda,v);

        T time_evalAA1 = 0.;
        Coefficients<Lexicographical,T,Index2D> IAv_ref, LIIAv_ref, UIv_ref, IAUIv_ref, AAv_ref;
        IndexSet<Index1D> checkLambda_x;
        for (const_set2d_it it=checkLambda.begin(); it!=checkLambda.end(); ++it) {
            checkLambda_x.insert((*it).index1);
            LIIAv[*it] = 0.;
            LIIAv_ref[*it] = 0.;
            IAUIv[*it] = 0.;
            IAUIv_ref[*it] = 0.;
            AAv_ref[*it] = 0.;
        }
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            for (const_set2d_it row=checkLambda.begin(); row!=checkLambda.end(); ++row) {
                Index1D row_x = (*row).index1;
                Index1D row_y = (*row).index2;
                if (     (row_x.xtype==XWavelet && col_x.xtype==XBSpline)
                      || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j > col_x.j)) {
                    Support<T> col_supp_x = basis.generator(col_x.xtype).support(col_x.j,col_x.k);
                    Support<T> row_supp_x = basis.generator(row_x.xtype).support(row_x.j,row_x.k);
                    if (overlap(col_supp_x,row_supp_x)>0) {
                        Index2D index(col_x,row_y);
                        IAv_ref[index] = 0.;
                    }
                }
            }
        }
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            for (const_set1d_it row=checkLambda_x.begin(); row!=checkLambda_x.end(); ++row) {
                Index1D row_x = (*row);
                if (     (row_x.xtype==XBSpline)
                      || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j <= col_x.j)) {
                    Support<T> col_supp_x = basis.generator(col_x.xtype).support(col_x.j,col_x.k);
                    Support<T> row_supp_x = basis.generator(row_x.xtype).support(row_x.j,row_x.k);
                    if (overlap(col_supp_x,row_supp_x)>0) {
                        Index2D index(row_x,col_y);
                        UIv_ref[index] = 0.;
                    }
                }
            }
        }
        cout << "Size of checkLambda: " << checkLambda.size() << endl;
        cout << "Size of Lambda:      " << Lambda.size() << endl;
        cout << "Size of IAv:         " << IAv_ref.size() << endl;
        cout << "Size of UIv:         " << UIv_ref.size() << endl;
        cout << "Size of v:           " << v.size() << endl;

        cout << "Reference calculation started..." << endl;
        refComputationIAv(Bil_y, v, IAv_ref);
        cout << "IAv_ref finished." << endl;
        refComputationLIIAv(Bil_x, IAv_ref, LIIAv_ref);
        cout << "LIIAv_ref finished." << endl;
        refComputationUIv(Bil_x, v, UIv_ref);
        cout << "UIv_ref finished." << endl;
        refComputationIAUIv(Bil_y, UIv_ref, IAUIv_ref);
        cout << "IAUIv_ref finished." << endl;
        refComputationAAv(Bil_x,Bil_y, v, AAv_ref);
        cout << "AAv_ref finished." << endl;
        cout << "Reference calculation finished." << endl;
        cout << "New scheme started..." << endl;
        time.start();
        localop2d.debug_evalAA(v, intermediate, LIIAv, IAUIv, IAv_ref, LIIAv_ref, UIv_ref, IAUIv_ref, AAv_ref);
        time.stop();
        time_evalAA1 = time.elapsed();
        cout << "New scheme finished." << endl;

    }


    return 0;
}


void
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma)
{
    int j0 = basis.j0;
    for (long k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        for (long k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D row(j0,k1,XBSpline);
            Index1D col(j0,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0+i2-1;
            for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
                Lambda.insert(Index2D(col,row));
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0+i1-1;
        for (int i2=1; i1+i2<=j; ++i2) {
            if (T(i1+i2)-gamma*max(i1,i2)>(1-gamma)*j) continue;
            int j2=j0+i2-1;
            for (long k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    Lambda.insert(Index2D(row,col));
                }
            }
        }
    }
    return;
}

void
getRandomCoefficientVector(const IndexSet<Index2D> &Lambda,
                           Coefficients<Lexicographical,T,Index2D> &coeff)
{
    for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        coeff[*it] = T(rand()) / T(RAND_MAX);
    }
    return;
}

void
refComputationIAv(const BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &IAv)
{
    for (coeff2d_it row=IAv.begin(); row!=IAv.end(); ++row) {
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                if (isL2Orthonormal_y) {
                    if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                        val +=  (*col).second;
                    }
                }
                else {
                    val +=  Bil_y(row_y,col_y) * (*col).second;
                }
            }
        }
        (*row).second = val;
    }
    return;
}


void
refComputationLIIAv(const BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &IAv,
                    Coefficients<Lexicographical,T,Index2D> &LIIAv)
{
    if (isL2Orthonormal_y) return;
    for (coeff2d_it row=LIIAv.begin(); row!=LIIAv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        for (const_coeff2d_it col=IAv.begin(); col!=IAv.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                if (     (row_x.xtype==XWavelet && col_x.xtype==XBSpline)
                      || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j > col_x.j)) {
                    val +=   Bil_x(row_x,col_x) * (*col).second;
                }
            }
        }
        (*row).second = val;
    }
}


void
refComputationUIv(const BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &UIv)
{
    for (coeff2d_it row=UIv.begin(); row!=UIv.end(); ++row) {
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                if (isL2Orthonormal_x) {
                    if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                        val += (*col).second;
                    }
                }
                else {
                    if (    (row_x.xtype==XBSpline) || ((row_x.xtype==XWavelet && col_x.xtype==XWavelet
                                     && row_x.j<=col_x.j)) ) {
                        val += Bil_x(row_x,col_x) * (*col).second;
                    }
                }
            }
        }
        (*row).second = val;
    }
    return;
}


void
refComputationIAUIv(const BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &UIv,
                    Coefficients<Lexicographical,T,Index2D> &IAUIv)
{
    for (coeff2d_it row=IAUIv.begin(); row!=IAUIv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        for (const_coeff2d_it col=UIv.begin(); col!=UIv.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                if (isL2Orthonormal_y) {
                    if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                        val +=   (*col).second;
                    }
                }
                else {
                    val +=   Bil_y(row_y,col_y) * (*col).second;
                }
            }
        }
        (*row).second = val;
    }
}


void
refComputationAAv(const BilinearForm_x &Bil_x, const BilinearForm_y &Bil_y,
                  const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AAv)
{
    for (coeff2d_it row=AAv.begin(); row!=AAv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (isL2Orthonormal_x && isL2Orthonormal_y) {
                if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k &&
                    row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                    val +=   (*col).second;
                }
            }
            else {
                val +=   Bil_x(row_x,col_x) * Bil_y(row_y,col_y) * (*col).second;
            }
        }
        (*row).second = val;
    }
}
