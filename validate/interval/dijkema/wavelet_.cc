#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

int
main(int argc, char *argv[])
{
    if (argc!=5) {
        cerr << "usage: " << argv[0] << " d d_ j k" << endl;
        exit(-1);
    }
    int d     = atoi(argv[1]);
    int d_    = atoi(argv[2]);
    int j     = atoi(argv[3]);
    int k     = atoi(argv[4]);

    Basis<T,Dual,Interval,Dijkema> basis_(d,d_);
    Wavelet<T,Dual,Interval,Dijkema> psi_(basis_);
    for (int i=basis_.rangeJ_(j).firstIndex(); i<=basis_.rangeJ_(j).lastIndex(); ++i) {
        std::cerr << basis_.M1_(j,_,i) << std::endl;
    }
    cout.precision(18);
    for (T x = 0.; x<=1.; x+=1./10240) {
        cout << x << " " << psi_(x,j,k) << endl;
    }
    return 0;
}
