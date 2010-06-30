#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;
typedef flens::DenseVector<flens::Array<T> > DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > FullColMatrix;


int
main()
{
    {
        int nj = 8, njp1 = 16;
        MRA<double,Primal,Periodic,CDF> mra(4,6);
        FullColMatrix DM0(njp1,nj,0,0);
        DenseVectorT e(nj,0);
        for (int c=0; c<nj; ++c) {
            e(c) = 1;
            DM0(_,c) = mra.M0*e;
            e(c) = 0;
        }
        cout << DM0 << endl;
    }

    {
        int nj = 8, njp1 = 16;
        MRA<double,Dual,Periodic,CDF> mra_(4,6);
        FullColMatrix DM0_(njp1,nj,0,0);
        DenseVectorT e(nj,0);
        for (int c=0; c<nj; ++c) {
            e(c) = 1;
            DM0_(_,c) = mra_.M0_*e;
            e(c) = 0;
        }
        cout << DM0_ << endl;
    }

    {
        int nj = 32, njp1 = 64;
        Basis<double,Primal,Periodic,CDF> basis(4,6);
        FullColMatrix DM1(nj,njp1,0,0);
        DenseVectorT e(njp1,0);
        for (int c=0; c<njp1; ++c) {
            e(c) = 1;
            DM1(_,c) = transpose(basis.M1)*e;
            e(c) = 0;
        }
        cout << DM1 << endl;
    }
    
    {
        int nj = 32, njp1 = 64;
        Basis<double,Dual,Periodic,CDF> basis_(4,6);
        FullColMatrix DM1_(nj,njp1,0,0);
        DenseVectorT e(njp1,0);
        for (int c=0; c<njp1; ++c) {
            e(c) = 1;
            DM1_(_,c) = transpose(basis_.M1_)*e;
            e(c) = 0;
        }
        cout << DM1_ << endl;
    }
    return 0;
}