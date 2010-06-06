#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

const int J = 3;
const int maxD  = 4;
const int maxD_ = 10; 

int
main(int argc, char *argv[])
{
/*
{ // primal B-splines vs. dual B-splines
    Param<BSpline<T,Dual,R,CDF> >::resolution = 10;
    for (int d=2; d<=maxD; ++d) {
        for (int d_=d; d_<=maxD_; d_+=2) {
            BSpline<T,Primal,R,CDF> phi(d);
            BSpline<T,Dual,R,CDF>   phi_(d,d_);

            Integral<T,CompositeTrapezoidal,
                     BSpline<T,Primal,R,CDF>,
                     BSpline<T,Dual,R,CDF> 
                    > integral(phi,phi_);

            int k2 = 0;
            T error = 0;
            int first = phi_.a.firstIndex()-phi.a.lastIndex()+1;
            int last  = phi_.a.lastIndex()-phi.a.firstIndex()-1;
            for (int k1=first; k1<=last; ++k1) {
                error += fabs(integral(J,k1,J,k2) - (k1==k2));
            }
            std::cerr << "(" << d << "," << d_  << "): " << error << std::endl;
        }
    }
}
*/

{ // primal wavelets vs. dual wavelets
    for (int d=2; d<=maxD; ++d) {
        for (int d_=d; d_<=maxD_; d_+=2) {
            BSpline<T,Primal,R,CDF> phi(d);
            Wavelet<T,Dual,R,CDF> psi_(d,d_);

            Integral<T,CompositeTrapezoidal,
                     BSpline<T,Primal,R,CDF>,
                     Wavelet<T,Dual,R,CDF> 
                    > integral(phi,psi_);

            int k2 = 0;
            T error = 0;
            int first = psi_.b_.firstIndex()-phi.a.lastIndex()+1;
            int last  = psi_.b_.lastIndex()-phi.a.firstIndex()-1;
            for (int k1=first; k1<=last; ++k1) {
                error += fabs(integral(J,k1,J,k2));
            }
            std::cerr << "(" << d << "," << d_  << "): " << error << std::endl;
        }
    }

}

/*
{ // primal wavelets vs. dual wavelets
    Param<BSpline<T,Dual,R,CDF> >::resolution = 10;
    for (int d=2; d<=maxD; ++d) {
        for (int d_=d; d_<=maxD_; d_+=2) {
            Wavelet<T,Primal,R,CDF> psi(d,d_);
            Wavelet<T,Dual,R,CDF> psi_(d,d_);

            Integral<T,CompositeTrapezoidal,
                     Wavelet<T,Primal,R,CDF>,
                     Wavelet<T,Dual,R,CDF> 
                    > integral(psi,psi_);

            int k2 = 0;
            T error = 0;
            int first = psi_.b_.firstIndex()-psi.b.lastIndex()+1;
            int last  = psi_.b_.lastIndex()-psi.b.firstIndex()-1;
            for (int k1=first; k1<=last; ++k1) {
                error += fabs(integral(J,k1,J,k2) - (k1==k2));
            }
            std::cerr << "(" << d << "," << d_  << "): " << error << std::endl;
        }
    }

}
*/
    return 0;
}
