#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Periodic, CDF>           PrimalBasis;

typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;
typedef GeMatrix<FullStorage<T,ColMajor> >		FullColMatrixT;

int main (int argc, char *argv[]) {

    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d d_ j0" << endl;
        return 0;
    }
    int d  = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    int j0 = atoi(argv[3]);
    cout.precision(18);

    PrimalBasis basis(d,d_,j0);

    cout << "************** MRA ***************" << endl;
    cout << basis.mra.M0.band << endl;
    cout << "Number of left MRA parts:  " << basis.mra.cardIL(j0) << endl;
    cout << "Number of right MRA parts: " << basis.mra.cardIR(j0) << endl;

    cout << endl;
    FullColMatrixT M0(basis.mra.cardI(j0+1), basis.mra.cardI(j0));
    for(int i = 0; i < basis.mra.cardI(j0); ++i){
        int pos = basis.mra.rangeI(j0).firstIndex()+i;
        DenseVectorT x(basis.mra.cardI(j0)), y, ytrans;
        x(i+1) = 1;
        mv(NoTrans, 1, basis.mra.M0, x, 0., y);
        M0(flens::_, i+1) = y;
    }
    cout << "M0: " << M0 << endl;

    FullColMatrixT M0_Trans(basis.mra.cardI(j0), basis.mra.cardI(j0+1));
    for(int i = 0; i < basis.mra.cardI(j0+1); ++i){
        int pos = basis.mra.rangeI(j0+1).firstIndex()+i;
        DenseVectorT x(basis.mra.cardI(j0+1)), ytrans;
        x(i+1) = 1;
        mv(Trans, 1, basis.mra.M0, x, 0., ytrans);
        M0_Trans(flens::_, i+1) = ytrans;
    }
    cout << "M0_T: " << M0_Trans << endl;

    cout << "**********************************" << endl << endl;


    cout << "*********** Wavelets *************" << endl;
    cout << basis.M1.band << endl;
    cout << "Number of left Basis parts:  " << basis.cardJL(j0) << endl;
    cout << "Number of right Basis parts: " << basis.cardJR(j0) << endl;

    cout << "Wavelet range: " << basis.rangeJ(j0) << endl;
    cout << "Scaling range: " << basis.mra.rangeI(j0+1) << endl;

    cout << endl;
    FullColMatrixT M1(basis.cardJ(j0+1), basis.cardJ(j0));
    for(int i = 0; i < basis.cardJ(j0); ++i){
        int pos = basis.rangeJ(j0).firstIndex()+i;
        DenseVectorT x(basis.cardJ(j0)), y;
        x(i+1) = 1;
        mv(NoTrans, 1, basis.M1, x, 0., y);
        M1(flens::_, i+1) = y;
    }
    cout << "M1: " << M1 << endl;

    FullColMatrixT M1_Trans(basis.cardJ(j0), basis.cardJ(j0+1));
    for(int i = 0; i < basis.cardJ(j0+1); ++i){
        int pos = basis.rangeJ(j0+1).firstIndex()+i;
        DenseVectorT x(basis.cardJ(j0+1)), y;
        x(i+1) = 1;
        mv(Trans, 1, basis.M1, x, 0., y);
        M1_Trans(flens::_, i+1) = y;
    }
    cout << "M1_T: " << M1_Trans << endl;

    /*GeMatrix<FullStorage<T,ColMajor> > Mj1, Mj1_;
    initial_stable_completion(basis.mra,basis.mra_,Mj1,Mj1_);
    cout << "Wavelet range: " << basis.rangeJ(j0) << endl;
    cout << "Scaling range: " << basis.mra.rangeI(j0+1) << endl;
    cout << "Refinement matrix M1 = " << Mj1 << endl;
    GeMatrix<FullStorage<T,ColMajor> > t_Mj1;
    //copy(cxxblas::Trans, Mj1, t_Mj1);
    //cout << "Transposed Refinement matrix M1^T = " << t_Mj1 << endl;

*/
    return 0;
}
