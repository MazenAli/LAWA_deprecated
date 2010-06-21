#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;
typedef BSpline<T, Primal, Periodic, CDF> PrimalSpline;

int main(int argc, char* argv[]){
    
    int d = 2;
    PrimalSpline phi1(d), phi2(d), d_phi1(d,1), d_phi2(d,1);
    
    Integral<T, Gauss, PrimalSpline, PrimalSpline> integral_sfsf(phi1, phi2),
                                                dd_integral_sfsf(d_phi1, d_phi2);
    
    /* VERSION A */
    T a_dd = dd_integral_sfsf(3, 0, 3, 1);
    T a = integral_sfsf(3,0,3,1);
    
    
    /* VERSION B */
    // T a = integral_sfsf(3,0,3,1);
    // T a_dd = dd_integral_sfsf(3, 0, 3, 1);    
    
    cout << "Integral Abl vs Abl : " << a_dd << endl;
    cout << "Integral Spline vs Spline : " << a << endl;

	return 0;
}