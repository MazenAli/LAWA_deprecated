#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

int
main()
{
    Basis<double,Primal,Periodic,CDF> basis(3,5);
    Basis<double,Dual,Periodic,CDF> basis_(3,5);
    DenseVectorT x(2,1),y;
   // x = 1,2,3,4,5,2,1,1;
    x = 1,0;
    cout << x << endl;
    fwt(x,basis_,0,y);
    cout << y << endl;
    ifwt(y,basis,0,x);
    cout << x << endl;
    
    return 0;
}