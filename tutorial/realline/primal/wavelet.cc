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

    Wavelet<T,Primal,R,CDF> psi(d,d_,deriv);

    Support<T> supp = psi.support(j,k);
    std::cerr << supp << std::endl;
    for (T x = supp.l1; x<=supp.l2; x+=1./1024) {
        cout.precision(18);
        cout << x << " " << psi(x,j,k) << endl;
    }

    return 0;
}
