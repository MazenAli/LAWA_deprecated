#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>
#include <memory>
#include <extensions/flens/jacobi_transform.h>

using namespace lawa;
using namespace std;

typedef double T;

const int M = 5;
int
main()
{
    T **A;
    A = (T **)malloc(M*sizeof(T *));
    for (int i=0; i<M; ++i) {
        A[i] = (T *)malloc(M*sizeof(T));
    }
    for (int i=0; i<M; ++i) {
        std::uninitialized_fill_n(A[i],M,T(0));
    }
    A[1][1] = 1;
    A[1][2] = 2;
    A[1][3] = 3;
    A[1][4] = 4;
    A[2][1] = 2;
    A[2][2] = -1;
    A[2][3] = 9;
    A[2][4] = 22;
    A[3][1] = 3;
    A[3][2] = 9;
    A[3][3] = 0;
    A[3][4] = 8;
    A[4][1] = 4;
    A[4][2] = 22;
    A[4][3] = 8;
    A[4][4] = 4;

    T d[M];
    
    T **v;
    v = (T **)malloc(M*sizeof(T *));
    for (int i=0; i<M; ++i) {
        v[i] = (T *)malloc(M*sizeof(T));
    }
    for (int i=0; i<M; ++i) {
        std::uninitialized_fill_n(v[i],M,T(0));
    }
    for (int i=0; i<M; ++i) {
        for (int j=0; j<M; ++j) {
            v[i][j] = 0.0;
        }
    }
    int n;
    jacobi(A,M-1,d,v,&n);

    cout.precision(256);
    for (int i=1; i<M; ++i) {
        cout << d[i] << "\t";
    }
    cout << endl << endl;
    for (int i=1; i<M; ++i) {
        for (int j=1; j<M; ++j) {
            cout << v[i][j] << "\t";
        }
        cout << endl;
    }
    
    return 0;
}
