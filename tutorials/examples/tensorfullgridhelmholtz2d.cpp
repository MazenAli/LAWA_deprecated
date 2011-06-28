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
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

/// Typedefs for problem components: +
/// Interval Basis using Dijkema construction +
/// TensorBasis: uniform (= full) basis
typedef Basis<T, Primal, Interval, Dijkema>                         Basis1D;
typedef TensorBasis2D<Uniform, Basis1D, Basis1D>                    Basis2D;

///      HelmholtzOperator in 2D, i.e. for $a(v,u) = \int(v_x \cdot u_x) + c \cdot \int(v \cdot u)$
typedef IdentityOperator1D<T, Basis1D>                              IdentityOp1D;
typedef HelmholtzOperator1D<T, Basis1D>                             HelmholtzOp1D;
typedef UniformTensorMatrix2D<T,Basis2D,HelmholtzOp1D,IdentityOp1D,
                              IdentityOp1D,HelmholtzOp1D>           TensorMatrix2D;
typedef HelmholtzOperator2D<T, Basis2D>                             HelmholtzOp2D;

///      Preconditioner: diagonal scaling with norm of operator
typedef DiagonalMatrixPreconditioner2D<T, Basis2D, HelmholtzOp2D>   DiagonalPrec;
///      Right Hand Side (RHS): basic 2D class for rhs integrals of the form $\int f \cdot v$,
///      with $f(x,y) = f_1(x) \cdot f_2(y)$, i.e. the RHS is separable.
//typedef SeparableRHS2D<T, Basis2D>                                  RHS2D;
//typedef SumOfTwoRHSIntegrals<T, Index2D, RHS2D, RHS2D>              SumOf2RHS;

typedef RHSWithPeaks1D<T, Basis1D>                                   RHS1D;

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
printUandCalcError(const DenseVectorT u, const Basis2D& basis, const int J_x, const int J_y,
                   const char* filename, const double deltaX=0.01, const double deltaY=0.01)
{
    T L2error = 0.;
    T H1error = 0.;
    ofstream file(filename);
    for(double x = 0.; x <= 1.; x += deltaX){
        for(double y = 0; y <= 1.; y += deltaY){
            T u_approx =    evaluate(basis, J_x, J_y, u, x, y, 0, 0);
            T dx_u_approx = evaluate(basis, J_x, J_y, u, x, y, 1, 0);
            T dy_u_approx = evaluate(basis, J_x, J_y, u, x, y, 0, 1);
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
    if(argc != 7){
        cerr << "Usage: " << argv[0] << " d d_ j0_x J_x j0_y J_y" << endl;
        exit(-1);
    }
    /// ... order of wavelets,
    int d = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    ///  ... minimal levels
    int j0_x = atoi(argv[3]);
    int j0_y = atoi(argv[5]);
    /// ... maximal levels.
    int J_x = atoi(argv[4]);
    int J_y = atoi(argv[6]);

    /// Basis initialization, setting Dirichlet boundary conditions in the second dimension
    Basis1D basis_x(d, d_, j0_x);
    Basis1D basis_y(d, d_, j0_y);
    basis_x.enforceBoundaryCondition<DirichletBC>();
    basis_y.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis_x, basis_y);

    cout << "Cardinality of basis: " << basis2d.dim(J_x, J_y) << endl;

    Timer time;
    time.start();

    /// Operator initialization
    cout << "Assembly of tensor matrix started." << endl;
    IdentityOp1D   identity_op_x(basis_x);
    HelmholtzOp1D  helmholtz_op_x(basis_x,0.5*c);
    IdentityOp1D   identity_op_y(basis_y);
    HelmholtzOp1D  helmholtz_op_y(basis_y,0.5*c);
    TensorMatrix2D TensorA(basis2d,helmholtz_op_x,identity_op_y,
                                   identity_op_x,helmholtz_op_y,J_x,J_y);

    /// Assemble A of linear system Au=f.
    cout << "Assembly of tensor matrix finished." << endl;

    /// Integral initialization for Preconditioner
    HelmholtzOp2D   a(basis2d, c);
    DiagonalPrec    p(a);
    Integral<Gauss,Basis1D,Basis1D> integral_x(basis_x,basis_x);
    Integral<Gauss,Basis1D,Basis1D> integral_y(basis_y,basis_y);

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
    RHS1D rhs_f1(basis_x, f1Fct, nodeltas, 10);
    RHS1D rhs_f2(basis_y, f2Fct, nodeltas, 10);
    RHS1D rhs_u1(basis_x, u1Fct, nodeltas, 10);
    RHS1D rhs_u2(basis_x, u1Fct, nodeltas, 10);

    DenseVectorT f(basis2d.dim(J_x,J_y));
    DenseVectorT P_vec(basis2d.dim(J_x,J_y));


    /// Assemble f of linear system Au=f and preconditioner P.

    cout << "Assembly of f and P started." << endl;
    int count=1;
    for (int k_x=basis_x.mra.rangeI(basis_x.j0).firstIndex(); k_x<=basis_x.mra.rangeI(basis_x.j0).lastIndex(); ++k_x) {
        T val_rhs_f1 = rhs_f1(XBSpline,basis_x.j0,k_x);
        T val_rhs_u1 = rhs_u1(XBSpline,basis_x.j0,k_x);

        T dd_x = integral_x(basis_x.j0,k_x,XBSpline,1,basis_x.j0,k_x,XBSpline,1);
        T id_x = integral_x(basis_x.j0,k_x,XBSpline,0,basis_x.j0,k_x,XBSpline,0);
        for (int k_y=basis_y.mra.rangeI(basis_y.j0).firstIndex(); k_y<=basis_y.mra.rangeI(basis_y.j0).lastIndex(); ++k_y) {
            T val_rhs_f2 = rhs_f2(XBSpline,basis_y.j0,k_y);
            T val_rhs_u2 = rhs_u2(XBSpline,basis_y.j0,k_y);

            T dd_y = integral_y(basis_y.j0,k_y,XBSpline,1,basis_y.j0,k_y,XBSpline,1);
            T id_y = integral_y(basis_y.j0,k_y,XBSpline,0,basis_y.j0,k_y,XBSpline,0);

            f(count) = val_rhs_f1*val_rhs_u2 + val_rhs_u1*val_rhs_f2;
            P_vec(count) = 1./(std::sqrt(dd_x*id_y + id_x*dd_y + id_x*id_y));
            ++count;
        }
        for (int j_y=basis_y.j0; j_y<J_y; ++j_y) {
            for (int k_y=basis_y.rangeJ(j_y).firstIndex(); k_y<=basis_y.rangeJ(j_y).lastIndex(); ++k_y) {
                T val_rhs_f2 = rhs_f2(XWavelet,j_y,k_y);
                T val_rhs_u2 = rhs_u2(XWavelet,j_y,k_y);

                T dd_y = integral_y(j_y,k_y,XWavelet,1,j_y,k_y,XWavelet,1);
                T id_y = integral_y(j_y,k_y,XWavelet,0,j_y,k_y,XWavelet,0);

                f(count) = val_rhs_f1*val_rhs_u2 + val_rhs_u1*val_rhs_f2;
                P_vec(count) = 1./(std::sqrt(dd_x*id_y + id_x*dd_y + id_x*id_y));
                ++count;
            }
        }
    }
    for (int j_x=basis_x.j0; j_x<J_x; ++j_x) {
        for (int k_x=basis_x.rangeJ(j_x).firstIndex(); k_x<=basis_x.rangeJ(j_x).lastIndex(); ++k_x) {
            T val_rhs_f1 = rhs_f1(XWavelet,j_x,k_x);
            T val_rhs_u1 = rhs_u1(XWavelet,j_x,k_x);

            T dd_x = integral_x(j_x,k_x,XWavelet,1,j_x,k_x,XWavelet,1);
            T id_x = integral_x(j_x,k_x,XWavelet,0,j_x,k_x,XWavelet,0);
            for (int k_y=basis_y.mra.rangeI(basis_y.j0).firstIndex(); k_y<=basis_y.mra.rangeI(basis_y.j0).lastIndex(); ++k_y) {
                T val_rhs_f2 = rhs_f2(XBSpline,basis_y.j0,k_y);
                T val_rhs_u2 = rhs_u2(XBSpline,basis_y.j0,k_y);

                T dd_y = integral_y(basis_y.j0,k_y,XBSpline,1,basis_y.j0,k_y,XBSpline,1);
                T id_y = integral_y(basis_y.j0,k_y,XBSpline,0,basis_y.j0,k_y,XBSpline,0);

                f(count) = val_rhs_f1*val_rhs_u2 + val_rhs_u1*val_rhs_f2;
                P_vec(count) = 1./(std::sqrt(dd_x*id_y + id_x*dd_y + id_x*id_y));
                ++count;
            }
            for (int j_y=basis_y.j0; j_y<J_y; ++j_y) {
                for (int k_y=basis_y.rangeJ(j_y).firstIndex(); k_y<=basis_y.rangeJ(j_y).lastIndex(); ++k_y) {
                    T val_rhs_f2 = rhs_f2(XWavelet,j_y,k_y);
                    T val_rhs_u2 = rhs_u2(XWavelet,j_y,k_y);

                    T dd_y = integral_y(j_y,k_y,XWavelet,1,j_y,k_y,XWavelet,1);
                    T id_y = integral_y(j_y,k_y,XWavelet,0,j_y,k_y,XWavelet,0);

                    f(count) = val_rhs_f1*val_rhs_u2 + val_rhs_u1*val_rhs_f2;
                    P_vec(count) = 1./(std::sqrt(dd_x*id_y + id_x*dd_y + id_x*id_y));
                    ++count;
                }
            }
        }
    }
    DiagonalMatrixT P(P_vec);
    cout << "Assembly of f and P finished." << endl;


    /// Solve system
    DenseVectorT u(basis2d.dim(J_x, J_y));
    cout << "Preconditioned cg solver started." << endl;
    int iterations = pcg(P, TensorA, u, f);

    time.stop();

    cout << "Preconditioned cg solver finished after " << iterations << " pcg iterations." << endl;

    /// Generate output for gnuplot visualization
    cout << "Post processing started." << endl;
    DenseVectorT error(2);
    //= printUandCalcError(u, basis2d, J_x, J_y, "u_helmholtz_pi.txt");
    cout << "Post processing finished (N,L2-error,H1-error,time): " << std::min(J_x-j0_x,J_y-j0_y)
         << " " << basis2d.dim(J_x,J_y) << " " << error(1) << " " << error(2)
         << " "<< time.elapsed() << endl;

    return 0;
}
