#include <iostream>
#include <lawa/lawa.h>
#include <lawa/methods/adaptive/algorithms/localrefinement.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef Basis<T,Dual,Interval,Dijkema>                              DualBasis;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

int main (int argc, char *argv[]) {

    if (argc!=6) {
        cout << "Usage: " << argv[0] << " d d_ j0 J bc" << endl;
        return 0;
    }
    int d  = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    int j0 = atoi(argv[3]);
    int J  = atoi(argv[4]);
    bool withDirichletBC = atoi(argv[5]);
    cout.precision(18);
    int j = J;

    PrimalBasis basis(d,d_,j0);
    if (withDirichletBC) basis.enforceBoundaryCondition<DirichletBC>();

    LocalRefinement<PrimalBasis> LocalRefine(basis,withDirichletBC);
    /*
    for (int k=basis.mra.rangeI(j).firstIndex(); k<=basis.mra.rangeI(j).lastIndex(); ++k) {
        LocalRefine.test_reconstruct(j,k,XBSpline);
        cout << "Hit enter" << endl;
        getchar();
    }
    for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
        LocalRefine.test_reconstruct(j,k,XWavelet);
        cout << "Hit enter" << endl;
        getchar();
    }
    */
    cout << "Number of left Basis parts:  " << basis.cardJL(j0) << endl;
    cout << "Number of right Basis parts: " << basis.cardJR(j0) << endl;
    GeMatrix<FullStorage<T,ColMajor> > Mj1, Mj1_;
    initial_stable_completion(basis.mra,basis.mra_,Mj1,Mj1_);
    cout << "Refinement matrix without boundary conditions:" << endl;
    cout << "Wavelet range: " << basis.rangeJ(j0) << endl;
    cout << "Scaling range: " << basis.mra.rangeI(j0+1) << endl;
    cout << "Refinement matrix M1 = " << Mj1 << endl;
    GeMatrix<FullStorage<T,ColMajor> > t_Mj1;
    //copy(cxxblas::Trans, Mj1, t_Mj1);
    //cout << "Transposed Refinement matrix M1^T = " << t_Mj1 << endl;

    cout << basis.mra.M0.rightband << endl;
    cout << "Number of left MRA parts:  " << basis.mra.cardIL(j0) << endl;
    cout << "Number of right MRA parts: " << basis.mra.cardIR(j0) << endl;
    for (int i=0; i<basis.mra.M0.left.length(); ++i) {
        int pos = basis.mra.M0.left.firstIndex()+i;
        cout << "Left part nr. " << i << " : " <<  basis.mra.M0.left(pos) << endl;
    }
    return 0;
}


