
#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

// Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>                                IntervalBasis;
//typedef Basis<T,Primal,Periodic,CDF>                                    PeriodicBasis;
typedef TensorBasis2D<Adaptive,IntervalBasis,IntervalBasis>             Basis2D;

// Operator Definitions
typedef NoPreconditioner<T, Index2D>                                    NoPrec2D;

// Data type definitions
typedef Coefficients<Lexicographical,T,Index2D>                         Coeffs;

T
sol_dummy(T t, T x){
    return 0;
}

int main(int argc, char* argv[]) {
    
    /* PARAMETERS: */

    if (argc  != 9) {
        cerr << "Usage " << argv[0] << " dt d_t j0_t dx d_x j0_x coeff_file outfile" << endl; 
        exit(1);
    }

    int dt       = atoi(argv[1]);
    int d_t      = atoi(argv[2]);
    int j0_t    = atoi(argv[3]);
    int dx       = atoi(argv[4]);
    int d_x      = atoi(argv[5]);
    int j0_x    = atoi(argv[6]);
    
    /* Basis Initialization */
    IntervalBasis   basis_t(dt, d_t, j0_t);
    IntervalBasis   basis_x(dx, d_x, j0_x);
    basis_x.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2d(basis_t, basis_x);
    
    NoPrec2D               noprec;

    
    Coeffs repr;
    readCoeffVector2D(repr, argv[7], false);
    plot2D(basis2d, repr, noprec, sol_dummy , 0., 1., 0., 1., 0.01, argv[8]);



    return 0;
}
