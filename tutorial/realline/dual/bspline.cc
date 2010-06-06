#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

int
main(int argc, char *argv[])
{
    if (argc!=5) {
        cerr << "usage: " << argv[0] << " d d_ j k" << endl;
        exit(-1);
    }
    int d  = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    int j  = atoi(argv[3]);
    int k  = atoi(argv[4]);

    BSpline<double,Dual,R> phi_(d,d_);
    
    Support<double> supp = phi_.support(j,k);
    for (double x = supp.l1; x<=supp.l2; x+=supp.length()/1024.) {
        cout << x << " " << phi_(x,j,k) << endl;
    }

    return 0;
}