#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

int
main(int argc, char *argv[])
{
    if (argc!=3) {
        cerr << "usage: " << argv[0] << " d d_" << endl;
        exit(-1);
    }
    int d  = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    
    GeMatrix<FullStorage<T,ColMajor> > DM0, DM0_;
    
    MRA<T,Primal,Interval,DKU> mra(d,d_);
    densify(cxxblas::NoTrans, mra.M0, DM0, mra.M0.firstRow(), mra.M0.firstCol());
    cout << "M0 = [" << DM0 << "];" << endl;

    MRA<T,Dual,Interval,DKU> mra_(d,d_);
    densify(cxxblas::NoTrans, mra_.M0_, DM0_, mra_.M0_.firstRow(), mra_.M0_.firstCol());
    cout << "M0_ = [" << DM0_ << "];" << endl;
    return 0;
}