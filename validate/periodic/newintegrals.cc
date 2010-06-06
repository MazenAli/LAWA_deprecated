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

typedef double T;

typedef BSpline<T,Primal,Periodic,CDF> PrimalSpline;
typedef BSpline<T,Dual,Periodic,CDF>   DualSpline;
typedef Wavelet<T,Primal,Periodic,CDF> PrimalWavelet;
typedef Wavelet<T,Dual,Periodic,CDF>   DualWavelet;


int
main(int argc, char *argv[])
{
    PrimalSpline phi1(3);
    DualSpline phi2(3,5);
    PrimalWavelet psi1(3,5);
    DualWavelet psi2(3,5);

    DenseVector<Array<double> > singularPoints;
    Function<double> Sin(f,singularPoints);
    
    cout.precision(18);
{
    Integral<double,Gauss,PrimalSpline,PrimalSpline> integral(phi1, phi1);
    cout << integral(0,1,0,0) << endl;
}
{
    Integral<double,CompositeTrapezoidal,DualSpline,DualSpline > integral(phi2, phi2);
    cout << integral(0,1,0,0) << endl;
}
{
    Integral<double,Gauss,PrimalWavelet,PrimalWavelet> integral(psi1, psi1);
    cout << integral(0,1,0,0) << endl;
}
{
    Integral<double,CompositeTrapezoidal,DualWavelet,DualWavelet> integral(psi2, psi2);
    cout << integral(0,1,0,0) << endl;
}
{
    Integral<double,CompositeTrapezoidal,PrimalSpline,DualSpline> integral(phi1, phi2);
    cout << integral(0,1,0,0) << endl;
}
{
    Integral<double,CompositeTrapezoidal,PrimalWavelet,DualWavelet> integral(psi1, psi2);
    cout << integral(0,1,0,0) << endl;
}

{
    Integral<double,CompositeTrapezoidal,PrimalSpline,Function<double> > integral(phi1, Sin);
    cout << integral(0,1) << endl;
}

{
    Integral<double,CompositeTrapezoidal,PrimalWavelet,Function<double> > integral(psi1, Sin);
    cout << integral(0,0) << endl;
}
{
    Integral<double,CompositeTrapezoidal,DualSpline,Function<double> > integral(phi2, Sin);
    cout << integral(0,1) << endl;
}
{
    Integral<double,CompositeTrapezoidal,DualWavelet,Function<double> > integral(psi2, Sin);
    cout << integral(1,0) << endl;
}
    return 0;
}