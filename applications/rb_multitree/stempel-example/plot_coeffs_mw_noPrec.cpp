
#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

// Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>                            Basis_X;
typedef Basis<T,Orthogonal,Interval,Multi>							Basis_Y;
typedef TensorBasis2D<Adaptive,Basis_X,Basis_Y>                     Basis2D;

// Operator Definitions
typedef H1NormPreconditioner2D<T, Basis2D>                          Prec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

// Data type definitions
typedef Coefficients<Lexicographical,T,Index2D>                     Coeffs;

T
sol_dummy(T t, T x){
    return 0;
}

int main(int argc, char* argv[]) {
    
    /* PARAMETERS: */

    if (argc  != 7) {
        cerr << "Usage " << argv[0] << " dx d_x j0_x dy coeff_file outfile" << endl; 
        exit(1);
    }

    int dx       = atoi(argv[1]);
    int d_x      = atoi(argv[2]);
    int j0_x    = atoi(argv[3]);
    int dy       = atoi(argv[4]);
    
    /* Basis Initialization */
    Basis_X   basis_x(dx, d_x, j0_x);
    Basis_Y   basis_y(dy, 0);
    basis_y.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2d(basis_x, basis_y);
    
    NoPrec2D            noprec;
    //Prec2D            prec(basis2d);

    
    Coeffs repr;
    readCoeffVector2D(repr, argv[5], false);
    //plot2D(basis2d, repr, prec, sol_dummy , 0., 1., 0., 1., 0.01, argv[6]);
    plot2D(basis2d, repr, noprec, sol_dummy , 0., 1., 0., 1., 0.01, argv[6]);



    return 0;
}
