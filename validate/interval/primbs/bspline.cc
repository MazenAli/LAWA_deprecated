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
        cerr << "usage: " << argv[0] << " d d_ j deriv" << endl;
        exit(-1);
    }
    int d      = atoi(argv[1]);
    int d_     = atoi(argv[2]);
    int j      = atoi(argv[3]);
    int deriv  = atoi(argv[4]);

    MRA<T,Primal,Interval,Primbs> mra(d,j);
    BSpline<T,Primal,Interval,Primbs> phi(mra, deriv);

    for (int k=mra.rangeI(j).firstIndex(); k<=mra.rangeI(j).lastIndex(); ++k) {
        T x = 0.0, 
          h = 1./1024;
        for (; x<=1; x += h) {
            cout << x << " " << phi(x,j,k) << endl;
        }
        cout << endl << endl;
    }
    
    return 0;
}
