#include <iostream>
#include <lawa/lawa.h>
#include <ctime>

using namespace std;
using namespace lawa;

typedef double T;


int
main(int argc, char* argv[]){
    
    int d = 2;
    BSpline<T, Primal, Periodic, CDF> phi(d);
    
    int j = 3;
    int k = 1;
    double x = 0.5;
    int N = 10000000;
    
    clock_t start;
    clock_t stop;
    
    cout << "j: " << j << ", k: " << k << endl << endl;
    
    cout << "=== Original ====== " << endl;
    cout << "Support: " << phi.phiR.support(j,k) << endl;
    cout << "Sing Supp: " << phi.phiR.singularSupport(j,k) << endl;
    cout << "Phi(" << x << ") = " << phi.phiR(x, j, k) << endl;
    
    cout << endl;
    start = clock();
    for(int i = 1; i <= N; i++){
        phi.phiR.singularSupport(j,k);
    }
    stop = clock();
    cout << N << " SingularSupports: " << double(stop - start) / CLOCKS_PER_SEC << " seconds" << endl;
    
    cout << endl << endl;
    
    cout << "=== Periodisch ==== " << endl;
    cout << "Support: " << phi.support(j,k) << endl;
    cout << "Sing Supp: " << phi.singularSupport(j,k) << endl;
    cout << "Phi(" << x << ") = " << phi(x, j, k) << endl;
    
    cout << endl;
    start = clock();
    for(int i = 1; i <= N; i++){
        phi.singularSupport(j,k);
    }
    stop = clock();
    cout << N << " SingularSupports: " << double(stop - start) / CLOCKS_PER_SEC << " seconds" << endl;
    
    return 0;
}