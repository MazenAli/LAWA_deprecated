#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef BSpline<T,Primal,Periodic,CDF> PrimalSpline;
typedef Wavelet<T,Primal,Periodic,CDF> PrimalWavelet;
typedef Basis<T, Primal, Periodic, CDF> PrimalBasis;
typedef TensorBasis<PrimalBasis, PrimalBasis> PeriodicTensorBasis;
typedef HelmholtzOperator<T, PeriodicTensorBasis> HelmholtzOp;
typedef SeparableRHS<T, PeriodicTensorBasis> PeriodicRHS;
typedef ScalFactorPreconditioner<T> Prec;
typedef BoxProblem<T, PeriodicTensorBasis, HelmholtzOp, PeriodicRHS, Prec> PeriodicBoxProblem; 

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

#define c  2.0

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
        
    HelmholtzOperator<T, PeriodicTensorBasis> a(basis, c);
    
    DenseVectorT singularSupport;
    SeparableFunction<T> rhs_fct(f_x, singularSupport, f_y, singularSupport);
    SeparableRHS<T, PeriodicTensorBasis> rhs(basis, rhs_fct);
    
    Prec P;
    
    PeriodicBoxProblem problem(basis, a, rhs, P);
    
    FullColMatrixT A = problem.getStiffnessMatrix(J_x, J_y);
    ofstream Afile("A.txt");
    Afile << A << endl;
    Afile.close();
    
    DenseVectorT f = problem.getRHS(J_x, J_y);
    ofstream ffile("f.txt");
    ffile << f << endl;
    ffile.close();
    
    FullColMatrixT D = problem.getPreconditioner(J_x, J_y);
    ofstream Dfile("D.txt");
    Dfile << D << endl;
    Dfile.close();
    
    
    
    return 0;
}