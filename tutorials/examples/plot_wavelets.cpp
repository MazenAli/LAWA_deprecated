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
///     Dijkema Basis over an interval
typedef Basis<T,Primal,Interval,Dijkema> Basis1D;
///     L2 orthonormal Basis over an interval
//typedef Basis<T,Orthogonal,Interval,Multi> Basis1D;

///  Typedefs for integral:
typedef Integral<Gauss,Basis1D,Basis1D> IntegralBasis1DvsBasis1D;

template<FunctionSide Side, Construction Cons>
Range<int>
getRange(const Basis<T,Side,Interval,Cons>& basis, int j, XType type);


int main (int argc, char *argv[]) {

    if(argc != 5){
        cerr << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(-1);
    }
    /// wavelet basis parameters:
    int d = atoi(argv[1]);
    int d_ =atoi(argv[2]);
    int j0 = atoi(argv[3]);
    int J = atoi(argv[4]);

    Basis1D basis(d,d_,j0);
    //Basis1D basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    IntegralBasis1DvsBasis1D integral(basis,basis);


    Range<int> scalingrange = getRange(basis, j0, XBSpline);
    /// Plot scaling function and wavelets
    ofstream plotfile_scaling("scaling.txt");
    for (T x=0.; x<=1.; x+=pow2i<T>(-12)) {
        plotfile_scaling << x;
        for (int k=scalingrange.firstIndex(); k<=scalingrange.lastIndex(); ++k) {
            plotfile_scaling << " " << basis.generator(XBSpline)(x,j0,k,0);
            if (x==0.) {
                cout << "||phi_{" << j0 << "," << k << "}||_L2 = "
                     << integral(j0,k,XBSpline,0, j0,k,XBSpline,0) << endl;
            }
        }
        plotfile_scaling << endl;
    }
    ofstream plotfile_wavelet("wavelet.txt");
    for (T x=0.; x<=1.; x+=pow2i<T>(-7)) {
        plotfile_wavelet << x;
        for (int j=j0; j<=J-1; ++j) {
            Range<int> waveletrange = getRange(basis, j, XWavelet);
            for (int k=waveletrange.firstIndex(); k<=waveletrange.lastIndex(); ++k) {
                plotfile_wavelet << " " << basis.generator(XWavelet)(x,j,k,0);
            }
        }
        plotfile_wavelet << endl;
    }


    return 0;
}


template<FunctionSide Side, Construction Cons>
Range<int>
getRange(const Basis<T,Side,Interval,Cons>& basis, int j, XType type)
{
    if (type==XBSpline) {
        return basis.mra.rangeI(j);
    }
    else {
        return basis.rangeJ(j);
    }
}

