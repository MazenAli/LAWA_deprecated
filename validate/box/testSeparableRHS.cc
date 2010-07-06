#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef flens::DenseVector<flens::Array<T> > DenseVectorT;
typedef Basis<T, Primal, Periodic, CDF> PrimalBasis;
typedef TensorBasis<PrimalBasis, PrimalBasis> PeriodicTensorBasis;


#define c 2.0

T f_x(T x)
{
    return (4.*M_PI*M_PI + c) * std::cos(2.*M_PI*x);
}

T f_y(T y)
{
    return 2.0;
}

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
	
    DenseVectorT singularSupport;
    Function<T> F_x(f_x, singularSupport);
    Function<T> F_y(f_y, singularSupport);
    
    
    cout << "========= Test SeparableFunction ===========" << endl;
    SeparableFunction<T> SepFct1(F_x, F_y);
    SeparableFunction<T> SepFct2(f_x, singularSupport, f_y, singularSupport);
    
    T x = 0;
    T y = 0;
    cout << "F1("<< x <<"," << y << ") = " << SepFct1.F_x(x)<< " * " 
         << SepFct1.F_y(y) << " = " << SepFct1(x,y) << endl;

    cout << "F2("<< x << "," << y << ") = " << SepFct2.F_x(x)<< " * " 
         << SepFct2.F_y(y) << " = " << SepFct2(x,y) << endl;
	
    cout << endl;
    
    
    x = 1;
    y = 1;
    cout << "F1("<< x << "," << y << ") = " << SepFct1.F_x(x)<< " * " 
         << SepFct1.F_y(y) << " = " << SepFct1(x,y) << endl;

    cout << "F2("<< x << "," << y << ") = " << SepFct2.F_x(x)<< " * " 
         << SepFct2.F_y(y) << " = " << SepFct2(x,y) << endl;
         
         
    cout << "========= Test SeparableRHS ===========" << endl;
    
    SeparableRHS<T, PeriodicTensorBasis> SepRHS(basis, SepFct2);
    
    cout << "RHS(Phi[0,0] x Phi[0,0]) = " << SepRHS(true, j0_x, 1, true, j0_y, 1) << endl;
    cout << "RHS(Psi[0,0] x Phi[0,0]) = " << SepRHS(false, j0_x , 1, true, j0_y, 1) << endl;
	
	return 0;
}