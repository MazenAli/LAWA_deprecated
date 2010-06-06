#include <cstdlib>
#include <iostream>
#include <lawa/flensforlawa.h>
#include <lawa/lawa.h>

using namespace flens;
using namespace lawa;
using namespace std;

typedef double  T;

int
main(int argc, char *argv[])
{
    if (argc!=4) {
        cerr << "usage: " << argv[0] << " d d_ J" << endl;
        exit(-1);
    }
    int d  = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    int J  = atoi(argv[3]);

    BSpline<T,Primal,R,CDF> phi(d);
    BSpline<T,Dual,R,CDF> phi_(d,d_);
    DenseVector<Array<T> > spline;

    cout.precision(18);
    subdivide(phi, J, spline);
    cout << phi.a << endl;
    cout << spline << endl;


    subdivide(phi_, J, spline);
    cout << phi_.a_ << endl;
    cout << spline << endl;

    return 0;
}
