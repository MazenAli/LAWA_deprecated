/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Typedef for double precision
typedef double T;

///  Typedefs for basis construction:
///     CDF Basis over for the realline
typedef Basis<T,Dual,R,CDF> Basis1D;


int main (int argc, char *argv[]) {

    if(argc != 3){
        cerr << "Usage: " << argv[0] << " d d_" << endl;
        exit(-1);
    }
    cout.precision(20);
    /// wavelet basis parameters:
    int d = atoi(argv[1]);
    int d_ =atoi(argv[2]);

    Basis1D basis(d,d_,0);


    ofstream plotfile_scaling("primal_scaling_generator.txt");
    T a=basis.generator(XBSpline).support(0,0).l1;
    T b=basis.generator(XBSpline).support(0,0).l2;
    for (T x=a; x<=b; x+=pow2i<T>(-12)) {
        plotfile_scaling << x << " " << basis.generator(XBSpline)(x,0,0,0) << endl;
    }

    ofstream plotfile_wavelet("primal_wavelet_generator.txt");
    a=basis.generator(XWavelet).support(0,0).l1;
    b=basis.generator(XWavelet).support(0,0).l2;
    for (T x=a; x<=b; x+=pow2i<T>(-12)) {
        plotfile_wavelet << x << " " << basis.generator(XWavelet)(x,0,0,0) << endl;
    }

    return 0;
}

