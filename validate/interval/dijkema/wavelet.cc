#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

int
main(int argc, char *argv[])
{
    if (argc!=6) {
        cerr << "usage: " << argv[0] << " d d_ j k deriv" << endl;
        exit(-1);
    }
    int d     = atoi(argv[1]);
    int d_    = atoi(argv[2]);
    int j     = atoi(argv[3]);
    int k     = atoi(argv[4]);
    int deriv = atoi(argv[5]);

    Basis<T,Primal,Interval,Dijkema> basis(d,d_);
    Wavelet<T,Primal,Interval,Dijkema> psi(basis,deriv);
    cout.precision(18);
    for (T x = 0.; x<=1.; x+=1./10240) {
        cout << x << " " << psi(x,j,k) << endl;
    }
    return 0;
}
