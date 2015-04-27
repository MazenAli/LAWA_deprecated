/* TEST UNIDIRECTIONAL OPERATOR
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

typedef double T;

///  Wavelet basis over an interval
typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

///  Underlying bilinear form
typedef LaplaceOperator1D<T,PrimalBasis>                           BilinearForm;
typedef RefinementBasis::LaplaceOperator1D                         RefinementBilinearForm1D;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementBilinearForm1D>                  LocalOp1D;

typedef UniDirectionalLocalOperator<Index2D,XOne,LocalOp1D,
                                           NotXOne,Index1D>        UniDirectionalLocalOpXOne2D;
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalOp1D,
                                           NotXTwo,Index1D>        UniDirectionalLocalOpXTwo2D;

typedef Coefficients<Lexicographical,T,Index2D>::iterator          coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator    const_coeff2d_it;

void
refComputationAIv(BilinearForm &Bil, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AIv);

void
refComputationIAv(BilinearForm &Bil, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &IAv);

int main (int argc, char *argv[]) {

    cout.precision(16);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);

    /// Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, j0);           // For L2_orthonormal and special MW bases
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    BilinearForm                   Bil(basis);
    LocalOp1D                      localOp1D(basis,basis,refinementbasis.LaplaceOp1D);
    UniDirectionalLocalOpXOne2D    uniDirectionalOpXOne2D(localOp1D);
    UniDirectionalLocalOpXTwo2D    uniDirectionalOpXTwo2D(localOp1D);

    Coefficients<Lexicographical,T,Index2D> v(SIZEHASHINDEX2D), IAv(SIZEHASHINDEX2D),
                                            AIv(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> IAv_ref(SIZEHASHINDEX2D),
                                            AIv_ref(SIZEHASHINDEX2D);


    getSparseGridVector(basis2d, v,  J, (T)0.);
    getSparseGridVector(basis2d, IAv, J, (T)0.);
    getSparseGridVector(basis2d, IAv_ref, J, (T)0.);
    getSparseGridVector(basis2d, AIv, J, (T)0.);
    getSparseGridVector(basis2d, AIv_ref, J, (T)0.);

    Index1D index1_x(j0+J+2,5,XWavelet);
    Index1D index1_y(j0+J+3,2,XWavelet);
    Index2D new_index1(index1_x,index1_y);
    completeMultiTree(basis2d, new_index1, v);
    Index1D index2_x(j0+J+3,3,XWavelet);
    Index1D index2_y(j0+J+2,4,XWavelet);
    Index2D new_index2(index2_x,index2_y);
    completeMultiTree(basis2d, new_index2, IAv);
    completeMultiTree(basis2d, new_index2, IAv_ref);
    completeMultiTree(basis2d, new_index2, AIv);
    completeMultiTree(basis2d, new_index2, AIv_ref);

    for (coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second = rand() / (T)RAND_MAX;
    }

    cout << "Reference calculation AIv started with #supp v = " << v.size() << endl;
    refComputationAIv(Bil, v, AIv_ref);
    cout << "... finished." << endl;

    cout << "Uni-directional evaluation AIv started with #supp v = " << v.size() << endl;
    uniDirectionalOpXOne2D.eval(v, AIv);
    cout << "... finished." << endl;
    AIv_ref -= AIv;
    cout << "Error: " << AIv_ref.norm((T)2.) << endl;

    cout << "Reference calculation IAv started with #supp v = " << v.size() << endl;
    refComputationIAv(Bil, v, IAv_ref);
    cout << "... finished." << endl;

    cout << "Uni-directional evaluation IAv started with #supp v = " << v.size() << endl;
    uniDirectionalOpXTwo2D.eval(v, IAv);
    cout << "... finished." << endl;
    IAv_ref -= IAv;
    cout << "Error: " << IAv_ref.norm((T)2.) << endl;

    return 0;

}

void
refComputationAIv(BilinearForm &Bil, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AIv)
{
    for (coeff2d_it row=AIv.begin(); row!=AIv.end(); ++row) {
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                val +=  Bil(row_x,col_x) * (*col).second;
            }
        }
        (*row).second = val;
    }
    return;
}

void
refComputationIAv(BilinearForm &Bil, const Coefficients<Lexicographical,T,Index2D> &v,
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
                val +=  Bil(row_y,col_y) * (*col).second;
            }
        }
        (*row).second = val;
    }
    return;
}
