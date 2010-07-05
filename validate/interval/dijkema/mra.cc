#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

int
main(int argc, char *argv[])
{
    if (argc!=4) {
        cerr << "usage: " << argv[0] << " d d_ j" << endl;
        exit(-1);
    }
    int d     = atoi(argv[1]);
    int d_    = atoi(argv[2]);
    int j     = atoi(argv[3]);

    MRA<T,Dual,Interval,Dijkema> mra_(d,d_,j);
    mra_.enforceBoundaryCondition<DirichletBC>();
//    GeMatrix<FullStorage<T,ColMajor> > DM0_;
//    densify(cxxblas::NoTrans, mra_.M0_,DM0_);
//    std::cerr << DM0_ << std::endl;
/*
    BSpline<T,Dual,Interval,Dijkema> phi_(mra_);
    
    cout.precision(18);

    cout << "plot '-' index 0:" << mra.cardI(j)-1 << " using 1:2 w l t''" << endl;
    for (int k=mra.rangeI(j).firstIndex(); k<=mra.rangeI(j).lastIndex(); ++k) {
        for (T x = 0.; x<=1.; x+=1./1024) {
            cout << x << " " << phi(x,j,k) << endl;
        }
        cout << endl << endl;
    }
*/
    return 0;
}
