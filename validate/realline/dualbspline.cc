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

    cout.precision(18);
    cerr.precision(18);

    BSpline<T,Dual,R,CDF> phi_(d,d_);

    for (T x = phi_.support(j,k).l1; x<=phi_.support(j,k).l2; x+=1./1024) {
        cout << x << " " << phi_(x,j,k) << endl;
    }
    return 0;
}
