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
        cerr << "usage: " << argv[0] << " d j k deriv" << endl;
        exit(-1);
    }
    int d     = atoi(argv[1]);
    int j     = atoi(argv[2]);
    int k     = atoi(argv[3]);
    int deriv = atoi(argv[4]);

    BSpline<T,Primal,R> phi(d,deriv);

    Support<T> supp = phi.support(j,k);
    cerr << supp << endl;
    for (T x = supp.l1; x<=supp.l2; x+=1./1024) {
        cout << x << " " << phi(x,j,k) << endl;
    }

    return 0;
}
