#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef BSpline<double,Primal,R,CDF> Spline;

int
main()
{
    int n = 10;
    Spline *AllSplines = (Spline *)malloc(n*sizeof(Spline));
    for (int i=0; i<n; ++i) {
        Spline *spline = new (AllSplines+i) Spline(i+1);
    }
    cout << AllSplines[8].d << endl;


    // oder aufgehübscht ...
    const Spline &currentSpline = AllSplines[3];
    cout << currentSpline.d << endl;


    // natürlich obliegt dir dann die Speicherverwaltung ...
    for (int i=1; i<n; ++i) {
        ((Spline *)AllSplines+i)->~Spline();
    }
    free(AllSplines);
    
    return 0;
}
