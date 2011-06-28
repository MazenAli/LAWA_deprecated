#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

// Basis definitions
typedef Basis<T,Primal,R,CDF>                   Basis1D;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

T
weight(T x)
{
    return x*x;
}

int main(int argc, char *argv[]) {

    cout.precision(10);
    if (argc != 5) {
        cout << "usage " << argv[0] << " j1 k1 j2 k2" << endl; exit(1);
    }
    int d =2;
    int d_=2;
    int j0=0;
    int j1=atoi(argv[1]);
    int k1 =atoi(argv[2]);
    int j2=atoi(argv[3]);
    int k2 =atoi(argv[4]);

    cout.precision(16);
    Basis1D basis(d, d_, j0);


    DenseVectorT singpts;
    Function<T> weightFct(weight, singpts);

    // 1D: works fine

    IntegralF<Gauss, Basis1D> integral1(weightFct, basis);
    cout << "<f,psi>: " << integral1(j1, k1, XWavelet, 0) << endl;
    cout << "<f,d_psi>: " << integral1(j1, k1, XWavelet, 1) << endl;

    IntegralF<Gauss, Basis1D> integral2(weightFct, basis, basis);
    cout << "<f,psi1*psi2>: " << integral2(j1, k1, XWavelet, 0, j2, k2, XWavelet, 0) << endl;
    cout << "<f,d_psi1*d_psi2>: " << integral2(j1, k1, XWavelet, 1, j2, k2, XWavelet, 1) << endl;


    return 0;
}
