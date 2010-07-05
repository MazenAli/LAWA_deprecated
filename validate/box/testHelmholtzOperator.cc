#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef BSpline<T,Primal,Periodic,CDF> PrimalSpline;
typedef Wavelet<T,Primal,Periodic,CDF> PrimalWavelet;
typedef Basis<T, Primal, Periodic, CDF> PrimalBasis;
typedef TensorBasis<PrimalBasis, PrimalBasis> PeriodicTensorBasis;

int main (int argc, char*argv[])
{
    if(argc != 7){
        cerr << "Usage: " << argv[0] << " j0_x J_x j0_y J_y d d_" << endl;
        return 0;
    }
    
    /* PARAMETERS: minimal level etc */
    
    int j0_x = atoi(argv[1]);
    int J_x = atoi(argv[2]);
    int j0_y = atoi(argv[3]);
    int J_y = atoi(argv[4]);
    int d = atoi(argv[5]);
    int d_ = atoi(argv[6]);
    
    PrimalBasis b1(d, d_, j0_x);
    PrimalBasis b2(d, d_, j0_y);
    PeriodicTensorBasis basis(b1, b2);
    
    T c = 2.0;
    
    HelmholtzOperator<T, PeriodicTensorBasis> a(basis, c);
    
    cout << "a( Phi_[0,0] x Phi_[0,0], Phi_[0,0] x Phi_[0,0]) = " 
        << a(true, true, true, true, j0_x, 1, j0_x, 1, j0_y, 1, j0_y, 1) << endl;
        
    cout << "a( Phi_[0,0] x Psi_[0,0], Phi_[0,0] x Psi_[0,0]) = " 
        << a(true, false, true, false, j0_x, 1, j0_x, 1, j0_y, 1, j0_y, 1) << endl;
    
    return 0;
}