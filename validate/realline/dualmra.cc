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

    BSpline<T,Dual,R,CDF> phi_(d,d_);

    cout << "x = [";
    for (T x = 0; x<=1; x+=1/10240.) {
        cout << x << " ";
    }
    cout << "];" << std::endl;
    for (int k=-15;k<=15;++k) {//-phi_.a.lastIndex()-1; k<=-phi_.a.firstIndex()+1; ++k) {
        cout << "phir_" << d << d_ << j << ((k<0) ? "m" : "") << abs(k) << " = [";
        for (T x = 0; x<=1; x+=1/10240.) {
            cout << phi_(x,j,k) << " ";
        }
        cout << "];" << endl;
    }
    
    return 0;
}
