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

    MRA<T,Dual,Interval,Dijkema> mra_(d,d_,j);  
    BSpline<T,Dual,Interval,Dijkema> phi_(mra_);
    std::cerr << mra_.R_Left.rows() << "x" << mra_.R_Left.cols() << std::endl;
    std::cerr << mra_.R_Right.rows() << "x" << mra_.R_Right.cols() << std::endl;
    std::cerr << mra_.R_Left << mra_.R_Right << std::endl;

    for (T x = 0.; x<=1.; x+=1./1024) {
        cout << x << " " << phi_(x,j,k) << endl;
    }
    return 0;
}
