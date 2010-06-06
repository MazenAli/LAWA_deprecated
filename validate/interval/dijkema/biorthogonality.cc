#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

const int maxD  = 4;
const int maxD_ = 9; 

int
main(int argc, char *argv[])
{

    Param<BSpline<T,Dual,R,CDF> >::resolution = 12;
{ // primal B-splines vs. dual B-splines
    for (int d=2; d<=maxD; ++d) {
        MRA<T,Primal,Interval,Primbs> mra(d);
        for (int d_=d; d_<=maxD_; d_+=2) {
            MRA<T,Dual,Interval,Dijkema> mra_(d,d_);
            int J = mra_.min_j0;
            BSpline<T,Primal,Interval,Primbs> phi(mra);
            BSpline<T,Dual,Interval,Dijkema>  phi_(mra_);

            Integral<T,CompositeTrapezoidal,
                     BSpline<T,Primal,Interval,Primbs>,
                     BSpline<T,Dual,Interval,Dijkema>
                    > integral(phi,phi_);

            for (int k1=mra.rangeI(J).firstIndex(); k1<=mra.rangeI(J).lastIndex(); ++k1) {
                for (int k2=mra_.rangeI_(J).firstIndex(); k2<=mra_.rangeI_(J).lastIndex(); ++k2) {
                    T value = integral(J,k1,J,k2);
                    std::cerr << "(" << k1 << "," << k2  << "): " << ((value < 1e-15) ? 0 : value) << std::endl;
                }
                std::cerr << std::endl;
            }
        }
    }
}

{ // primal wavelets vs. dual wavelets
    for (int d=4; d<=maxD; ++d) {
        for (int d_=d; d_<=maxD_; d_+=2) {
            std::cerr << "---- " << d << "," << d_ << " -------------" << std::endl;
            Basis<T,Primal,Interval,Dijkema> basis(d,d_);
            Basis<T,Dual,Interval,Dijkema> basis_(d,d_);
            int J = basis_.min_j0;
            Wavelet<T,Primal,Interval,Dijkema> psi(basis);
            Wavelet<T,Dual,Interval,Dijkema>  psi_(basis_);

            Integral<T,CompositeTrapezoidal,
                     Wavelet<T,Primal,Interval,Dijkema>,
                     Wavelet<T,Dual,Interval,Dijkema>
                    > integral(psi,psi_);

            for (int k1=basis.rangeJ(J).firstIndex(); k1<=basis.rangeJ(J).lastIndex(); ++k1) {
                for (int k2=basis_.rangeJ_(J).firstIndex(); k2<=basis_.rangeJ_(J).lastIndex(); ++k2) {
                    T value = integral(J,k1,J,k2);
                    std::cerr << "(" << k1 << "," << k2  << "): " << ((value < 1e-15) ? 0 : value) << std::endl;
                }
                std::cerr << std::endl;
            }
        }
    }
}

    return 0;
}
