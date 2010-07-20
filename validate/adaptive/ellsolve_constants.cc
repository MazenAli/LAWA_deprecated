#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
const Construction Cons=Dijkema;
typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > FullColMatrix;

int main(int argc, char *argv[])
{
    cout << endl << argv[0] << " start" << endl << endl;

    if (argc < 5) {
        cout << "usage: " << argv[0] << " d d_ J example [beta]" << endl;
        exit(1);
    }

    int d=atoi(argv[1]), d_=atoi(argv[2]), J_max=atoi(argv[3]);
    Basis<T,Primal,Interval,Cons> basis(d,d_);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis<T,Dual,Interval,Cons> basis_(d,d_);
    basis_.enforceBoundaryCondition<DirichletBC>();
    const int &j0=basis.j0;


    HelmholtzExamples<T>::setExample(atoi(argv[4]),argc>5 ? atof(argv[5]):0, 
                                                   argc>6 ? atof(argv[6]):0);
    Helmholtz<T,Cons> problem(basis,basis_);
/*
    HypersingularExamples<T>::setExample(atoi(argv[4]));
    Hypersingular<T,Cons> problem(basis);
*/
    cout << problem.name() << ", d = " << d << ", d_ = " << d_
        << ", J = " << J_max  << endl << endl;
    cout << argv[0] << " ready..." << endl; getchar();

    IndexSet<T,Cons> Lambda(basis), LambdaCheck(basis);
    for (int J=j0; J<=J_max; ++J) {
        // scaling functions first ...
        if (J==j0) {
            for (WaveletIndex<T,Cons> lambda(basis); lambda.xtype==XBSpline; ++lambda) {
                Lambda.insert(lambda);
                LambdaCheck.insert(lambda);
            }
        }
        // ... then wavelets
        if (J>j0) {
            WaveletIndex<T,Cons> lambda(basis,J-1,basis.rangeJ(J-1).firstIndex(),XWavelet);
            for (; lambda.j==J-1; ++lambda) {
                Lambda.insert(lambda);
                LambdaCheck.insert(lambda);
            }
        }
        cout << "Matrix setup for d = " << d << ", d_ = " << d_
            << ", j0 = " << j0 << ", J = " << setw(2) << J << flush;
        problem.LambdaCheck=LambdaCheck;
        problem.Lambda=Lambda;
        FullColMatrix A;
        
        int m = problem.A.numRows();
        int n = problem.A.numCols();
        A.engine().resize(m,n,1,1);
        DenseVector<Array<T> > e(n);
        for (int i=1; i<=n; ++i) {
            e(i) = 1.;
            A(_,i) = problem.A*e;
            e(i) = 0.;
        }
        std::cerr << "A = [" << A << "];" << std::endl;
/*        
//        densify(cxxblas::NoTrans,problem.A,A);
        cout << " done, N = " << setw(6) << A.numCols() << flush;

        DenseVector<Array<T> > s;
        GeMatrix<FullStorage<T,ColMajor> > U, VT;
        svd(A, s, U, VT);
        T cA = std::abs(s(s.firstIndex())), CA = std::abs(s(s.firstIndex()));
        for (int i = s.firstIndex(); i <= s.lastIndex(); ++i) {
            CA = std::max(std::abs(s(i)), CA);
            cA = std::min(std::abs(s(i)), cA);
        }
        
        cout << ", CA    = " << setw(10) << CA
             << ", cA    = " << setw(10) << cA
             << ", kappa = " << setw(10) << CA/cA << endl << flush;
*/
    }
}

