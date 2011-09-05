/* HELMHOLTZ PROBLEM 2D
 *
 * This example calculates a Helmholtz problem on the two-dimensional domain [0,1]
 * with periodic boundary conditions in the first dimension and Dirichlet conditions
 * in the second dimension, i.e.
 *     - u_xx + c * u = f on (0,1)^2 ,
 *         u(1,y) = u(0,y)     for y in [0,1],
 *         u(x,1) = u(x,0) = 0 for x in [0,1]
 * The solution is obtained using a uniform tensor Wavelet-Galerkin method with a
 * diagonal scaling preconditioner.
 */

/// _The complete example code in snippets:_
#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Typedefs for Flens data types:
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >   FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >     SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                     DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;

/// Typedefs for problem components: +
/// Interval Basis using Dijkema construction +
/// TensorBasis: uniform (= full) basis
typedef Basis<T, Primal, Interval, Dijkema>                          Basis1D;
typedef TensorBasis2D<Uniform, Basis1D, Basis1D>                     Basis2D;

///      HelmholtzOperator in 2D, i.e. for $a(v,u) = \int(v_x \cdot u_x) + c \cdot \int(v \cdot u)$
typedef IdentityOperator1D<T, Basis1D>                               IdentityOp1D;
typedef HelmholtzOperator1D<T, Basis1D>                              HelmholtzOp1D;

typedef TensorSparseGrid2D<T, Basis2D, HelmholtzOp1D, IdentityOp1D,
                            IdentityOp1D, HelmholtzOp1D>             TensorSparseGrid;

typedef RHSWithPeaks1D<T, Basis1D>                                   RhsIntegral1D;

/// The constant 'c' in the Helmholtz equation.
const double c =  1.;


T u1(T x)
{
    return x*x*(1-x)*(1-x);
}

T u2(T y)
{
    return y*y*(1-y)*(1-y);
}

T du1(T x) {
    return (2*x-6*x*x+4*x*x*x);
}

T du2(T y) {
    return (2*y-6*y*y+4*y*y*y);
}

T ddu1(T x)
{
    return (2-12*x+12*x*x);
}

T ddu2(T y)
{
    return (2-12*y+12*y*y);
}

T f1(T x) {
    return -ddu1(x) + 0.5*c* u1(x);
}

T f2(T y) {
    return -ddu2(y) + 0.5*c* u2(y);
}

T sol(T x, T y)
{
    return u1(x) * u2(y);
}

T dx_sol(T x, T y) {
    return du1(x)*u2(y);
}

T dy_sol(T x, T y) {
    return u1(x)*du2(y);
}

/// again we provide a function that writes the results into a file for later visualization.

DenseVectorT
printUandCalcError(const DenseVectorT u, const TensorSparseGrid &sg,
                   const char* filename, const double deltaX=0.01, const double deltaY=0.01)
{
    T L2error = 0.;
    T H1error = 0.;
    ofstream file(filename);
    for(double x = 0.; x <= 1.; x += deltaX){
        for(double y = 0; y <= 1.; y += deltaY){
            T u_approx =    sg.evaluate(u,x,y,0,0);
            T dx_u_approx = sg.evaluate(u,x,y,1,0);
            T dy_u_approx = sg.evaluate(u,x,y,0,1);
            file << x << " " << y << " " << u_approx << " " << sol(x,y) << " "
                 << dx_u_approx << " " << dy_u_approx << endl;
            T factor = deltaX * deltaY;
            if((x == 0) || (x == 1.)){
                factor *= 0.5;
            }
            if((y == 0) || (y == 1.)){
                factor *= 0.5;
            }
            L2error += factor * (u_approx -  sol(x, y)) * (u_approx - sol(x,y));
            H1error += factor * (dx_u_approx - dx_sol(x,y))*(dx_u_approx - dx_sol(x,y))
                    +  factor * (dy_u_approx - dy_sol(x,y))*(dy_u_approx - dy_sol(x,y));
        }
        file << endl;
    }
    file.close();

    H1error += L2error;
    L2error = sqrt(L2error);
    H1error = sqrt(H1error);
    DenseVectorT error(2);
    error = L2error, H1error;
    return error;
}

/// First we check our arguments ...
int main (int argc, char*argv[])
{
    /* PARAMETERS: */
    if(argc != 5){
        cerr << "Usage: " << argv[0] << " d d_ j0 I" << endl;
        exit(-1);
    }
    /// ... order of wavelets,
    int d = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    ///  ... minimal level
    int j0 = atoi(argv[3]);
    /// ...  maximum increment level.
    int I = atoi(argv[4]);

    /// Basis initialization, setting Dirichlet boundary conditions in the second dimension
    Basis1D basis(d, d_, j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis, basis);

    Timer time, total_time;

    total_time.start();

    /// Operator initialization
    IdentityOp1D   identity_op_x(basis);
    HelmholtzOp1D  helmholtz_op_x(basis,0.5*c);
    IdentityOp1D   identity_op_y(basis);
    HelmholtzOp1D  helmholtz_op_y(basis,0.5*c);

    /// Assemble A of linear system Au=f
    cout << "Assembly of tensor sparse grid matrix started." << endl;
    time.start();
    TensorSparseGrid sg_A(basis2d, helmholtz_op_x,identity_op_y,
                                   identity_op_x, helmholtz_op_y, I, 0.25);
    time.stop();
    cout << "Assembly of tensor sparse grid matrix finished after " << time.elapsed() << endl;

    /// Assemble Preconditioner for linear system Au=f
    cout << "Assembly of tensor sparse grid preconditioner started." << endl;
    time.start();
    DiagonalMatrixT sg_P = sg_A.assembleDiagonalMatrixPreconditioner();
    time.stop();
    cout << "Assembly of tensor sparse grid preconditioner finished after " << time.elapsed() << endl;


    /// Right Hand Side:
    ///     Singular Supports in both dimensions
    DenseVectorT sing_support_x;
    DenseVectorT sing_support_y;
    Function<T> f1Fct(f1, sing_support_x);
    Function<T> f2Fct(f2, sing_support_y);
    Function<T> u1Fct(u2, sing_support_x);
    Function<T> u2Fct(u2, sing_support_y);
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    RhsIntegral1D rhs_f1(basis, f1Fct, nodeltas, 4);
    RhsIntegral1D rhs_f2(basis, f2Fct, nodeltas, 4);
    RhsIntegral1D rhs_u1(basis, u1Fct, nodeltas, 4);
    RhsIntegral1D rhs_u2(basis, u1Fct, nodeltas, 4);

    /// Assemble f of linear system Au=f
    cout << "Assembly of tensor sparse grid rhs started." << endl;
    time.start();
    DenseVectorT sg_f = sg_A.assembleRHS(rhs_f1, rhs_u2);
    sg_f = sg_f + sg_A.assembleRHS(rhs_u1, rhs_f2);
    time.stop();
    cout << "Assembly of tensor sparse grid rhs finished after " << time.elapsed() << endl;



    /// Solve system
    DenseVectorT sg_u(sg_A.getDimension());
    cout << "Preconditioned cg solver started." << endl;
    int iterations = pcg(sg_P, sg_A, sg_u, sg_f);

    total_time.stop();

    cout << "Preconditioned cg solver finished after " << iterations << " pcg iterations." << endl;

    /// Generate output for gnuplot visualization
    cout << "Post processing started." << endl;
    DenseVectorT error = printUandCalcError(sg_u, sg_A, "tensorsparsegridhelmholtz2d.txt");
    cout << "Post processing finished (N,L2-error,H1-error,time): " << I
         << " " << sg_A.getDimension() << " " << error(1) << " " << error(2)
         << " "<< total_time.elapsed() << endl;

    return 0;
}
