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

    cout.precision(18);
    cerr.precision(18);

    MRA<T,Dual,Interval,Dijkema> mra_(d,d_,j);  
    BSpline<T,Dual,Interval,Dijkema> phi_(mra_);
    cout << "rangeI_ = " << mra_.rangeI_(j) << ";" << std::endl;
    cout << "R_Left = [" << mra_.R_Left << "];" << std::endl;
    cout << "R_Right = [" << mra_.R_Right << "];" << std::endl;
    
    cout << "x" << d << d_ << " = [";
    for (T x = 0.; x<=1.; x+=1./1024) {
        cout << x << " ";
    }
    cout << "];" << std::endl;

    for (int k=mra_.rangeI_(j).firstIndex(); k<=mra_.rangeI_(j).lastIndex(); ++k) {
        cout << "phi_" << d << d_ << j << k << " = [";
        for (T x = 0.; x<=1.; x+=1./1024) {
            cout << phi_(x,j,k) << " ";
        }
        cout << "];" << endl;
    }
    return 0;
}
