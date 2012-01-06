#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> > DenseVectorT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >        SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >      DenseMatrixT;

//Iterator definitions
typedef IndexSet<Index1D>::const_iterator                               const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator         const_coeff1d_it;

typedef Basis<T,Primal,RPlus,SparseMulti> Basis1D;
typedef IntegralF<Gauss,Basis1D>                                      IntegralFRPlus;
typedef Integral<Gauss,Basis1D,Basis1D>                       IntegralRPlus;

T
U(T x)
{
    if (0<=x && x<=1) {
        return -x*(x-2);
    }
    else if (x>1) {
        return exp(-(x-1)*(x-1));
    }
    else {
        cerr << "   error in u" << endl;
        exit(1);
    }
}

int main (int argc, char *argv[]) {

    if (argc!=3) {
        cout << "usage " << argv[0] << " j0 J" << endl; exit(1);
    }
    cout.precision(10);
    int d=4;
    int j0=atoi(argv[1]);
    int J =atoi(argv[2]);
    int deriv=0;

    Basis1D basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    ofstream file1("sparsemulti_bspline_rplus.dat");
    T l1=basis.mra.phi.support(j0,0).l1;
    T l2=basis.mra.phi.support(j0,7).l2;
    for (T x=l1; x<=l2; x+=pow2i<T>(-9)) {
        file1 << x << " " << basis.mra.phi(x,j0,0,deriv) << " " << basis.mra.phi(x,j0,1,deriv)
                  << " " << basis.mra.phi(x,j0,2,deriv) << " " << basis.mra.phi(x,j0,3,deriv)
                  << " " << basis.mra.phi(x,j0,4,deriv) << " " << basis.mra.phi(x,j0,5,deriv)
                  << " " << basis.mra.phi(x,j0,6,deriv) << " " << basis.mra.phi(x,j0,7,deriv) << endl;
    }
    file1.close();

    for (int j=j0; j<=2; ++j) {
        for (int k=1; k<=16; ++k) {
            ofstream file2("sparsemulti_wavelet_rplus.dat");
            l1=basis.psi.support(j,k).l1;
            l2=basis.psi.support(j,k).l2;
            cout << "Support psi_(" << j << ", " << k << ") = " << basis.psi.singularSupport(j,k) << endl;
            for (T x=l1; x<=l2; x+=pow2i<T>(-9)) {
                file2 << x << " " << basis.psi(x,j,k,deriv) << endl;
            }
            file2.close();
            getchar();
        }
    }

    /*
    IndexSet<Index1D> Lambda;
    for (int k=0; k<=80; ++k) {
        Lambda.insert(Index1D(j0,k,XBSpline));
    }
    for (int j=j0; j<=J; ++j) {
        for (int k=1; k<=80; ++k) {
            Lambda.insert(Index1D(j,k,XWavelet));
        }
    }

    int N = Lambda.size();
    cout << "N = " << N << endl;
    SparseMatrixT A(N,N);

    DenseVectorT rhs(N), u(N);
    DenseVectorT singpts(1);
    singpts= 1.;
    Function<T> f_fct(U,singpts);
    IntegralFRPlus integralf(f_fct,basis);
    integralf.quadrature.setOrder(20);
    IntegralRPlus integral(basis,basis);

    int row_count=1, col_count=1;
    for (const_set1d_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
        col_count=1;
        rhs(row_count) = integralf((*row).j,(*row).k,(*row).xtype,0);
        for (const_set1d_it col=Lambda.begin(); col!=Lambda.end(); ++col,++col_count) {
            T tmp = integral((*row).j,(*row).k,(*row).xtype,0, (*col).j,(*col).k,(*col).xtype,0);
            if (fabs(tmp)>0) {
                A(row_count,col_count) = tmp;
            }
        }
    }
    A.finalize();
    int number_of_iterations = lawa::cg(A, u, rhs, 1e-12, 100);

    ofstream plotfile("u.txt");
    for (T x=0.; x<=10; x+=0.0001) {
        row_count = 1;
        T ret = 0.;
        for (const_set1d_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
            if ((*row).xtype==XBSpline) {
                ret += u(row_count) * basis.mra.phi(x,(*row).j,(*row).k,0);
            }
            else {
                ret += u(row_count) * basis.psi(x,(*row).j,(*row).k,0);
            }
        }
        plotfile << x << " " << U(x) << " " << ret << endl;
    }
    plotfile.close();
    */
    return 0;
}


