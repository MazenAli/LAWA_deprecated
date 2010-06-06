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
    
    GeMatrix<FullStorage<T,ColMajor> > DM0;
    
    MRA<T,Primal,Interval,Primbs> mra(d,d_);
    densify(mra.M0, DM0, mra.M0.firstRow(), mra.M0.firstCol());
    cout << DM0.rows() << "x" << DM0.cols() << DM0 << endl;
    mra.enforceBoundaryCondition<DirichletBC>();
    densify(mra.M0, DM0, mra.M0.firstRow(), mra.M0.firstCol());
    cout << DM0.rows() << "x" << DM0.cols() << DM0 << endl;
    return 0;
}