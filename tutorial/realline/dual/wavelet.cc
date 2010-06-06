#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

int
main(int argc, char *argv[])
{
    if (argc!=3) {
        cerr << "usage: " << argv[0] << " d d_" << endl;
        exit(-1);
    }
    int d     = atoi(argv[1]);
    int d_    = atoi(argv[2]);

    Wavelet<double,Dual,R,CDF> psi_(d,d_);

    cout << psi_.support(0,0) << endl;

    return 0;
}
