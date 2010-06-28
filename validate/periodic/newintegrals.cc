#include <cmath>
#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

double
f(double x)
{
    return std::sin(2*M_PI*x);
}

double
f2(double x)
{
    return x;
}

typedef double T;

typedef BSpline<T,Primal,Periodic,CDF> PrimalSpline;
typedef BSpline<T,Primal,R,CDF> PrimalRSpline;


typedef BSpline<T,Dual,Periodic,CDF>   DualSpline;
typedef Wavelet<T,Primal,Periodic,CDF> PrimalWavelet;
typedef Wavelet<T,Primal,R,CDF> PrimalRWavelet;
typedef Wavelet<T,Dual,Periodic,CDF>   DualWavelet;

void printSpline( PrimalSpline phi,int j,int k, const char* filename){
    ofstream plotFile(filename);
    for(double x = phi.support(j,k).l1; x <= phi.support(j,k).l2; x += 0.001){
        plotFile << x << " " << phi(x, j, k) << endl;
    }
    plotFile.close();
}

template <typename FctType>
void printFunction(FctType F, const char* filename, T a, T b, T dx){
    ofstream plotFile(filename);
    for(T x = a; x <= b; x += dx ){
        plotFile << x << " " << F(x) << endl;
    }
    plotFile.close();
}

template <typename IntegralType>
void printIntegrand(IntegralType I, const char* filename, T a, T b, T dx){
    
    ofstream plotFile(filename);
    for(T x = a; x <= b; x += dx ){
        plotFile << x << " " << I.integrand(x) << endl;
    }
    plotFile.close();
}

int
main(int argc, char *argv[])
{
    
    PrimalSpline phi1(2);
    DualSpline phi2(3,5);
    PrimalWavelet psi1(3,5);
    DualWavelet psi2(3,5);

    DenseVector<Array<double> > singularPoints(2);
    singularPoints = 0., 1.;
    Function<double> Sin(f,singularPoints);
    cout << "SingPoints Sin: " << Sin.singularPoints << endl;
    
    Function<double> X(f2, singularPoints);
    
    cout.precision(18);
    
{
    Integral<double, CompositeTrapezoidal, PrimalSpline, Function<double> > integral(phi1, X);
    cout << integral(0,0) << endl << endl;
    Integral<double, Gauss, PrimalSpline, Function<double> > integral2(phi1, X);
    cout << integral2(0,0) << endl << endl;
}    
    

{
    Integral<double,Gauss,PrimalSpline,PrimalSpline> integral(phi1, phi1);
    printSpline(phi1, 1, 1, "Phi1_1_1.txt");
    printSpline(phi1, 1, 0, "Phi1_1_0.txt");
    cout << integral(1,1,1,0) << endl;
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
    Integral<double, CompositeTrapezoidal, PrimalRSpline, Function<double> > integralR(phi1.phiR, Sin);
    printSpline(phi1, 0, 1,  "Phi1_0_1.txt");
    printFunction<Function<T> >(Sin, "Sin.txt", 0, 1, 0.001);
    printIntegrand<Integral<double,CompositeTrapezoidal,PrimalSpline,Function<double> > >(integral, "Integrand.txt", 0, 1, 0.001);

    cout << "Spline * Sin     : " << integral(1,1) << endl;
    cout << "Spline * Sin on R: " << integralR(1,1) << endl;

}

{
    Integral<double,CompositeTrapezoidal,PrimalWavelet,Function<double> > integral(psi1, Sin);
    Integral<double,CompositeTrapezoidal,PrimalRWavelet,Function<double> > integralR(psi1.psiR, Sin);
    cout << "Wavelet * Sin     : " << integral(0,1) << endl;
    cout << "Wavelet * Sin on R: " << integralR(0,1) << endl;
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