#include <cmath>
#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

double
f(double x)
{
    return std::sin(x);
}



int
main(int argc, char *argv[])
{
    typedef BSpline<double,Primal,R,CDF> PrimalSpline;
    typedef Wavelet<double,Primal,R,CDF> PrimalWavelet;
    PrimalSpline phi1(2), phi2(3);
    PrimalWavelet psi1(3,5), psi2(3,5);

    DenseVector<Array<double> > singularPoints;
    Function<double> Sin(f,singularPoints);
    
    cout.precision(18);

    Integral<double,Gauss,PrimalSpline,PrimalSpline> integral1(phi1, phi2);
    cout << integral1(0,1,0,0) << endl;

    Integral<double,Gauss,PrimalWavelet,PrimalSpline> integral2(psi1, phi2);
    cout << integral2(1,1,0,2) << endl;
    
    Integral<double,CompositeTrapezoidal,PrimalSpline,Function<double> > integralf(phi1,Sin);
    cout << integralf(0,0) << endl;

    return 0;
}