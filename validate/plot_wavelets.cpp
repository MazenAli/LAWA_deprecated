#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

typedef Basis<T,Primal,Interval,Dijkema> Basis1D;

Range<int>
getRange(const Basis1D& basis, int j, XType type);


int main (int argc, char *argv[]) {

    cout.precision(16);
    int d =2;
    int d_=2;
    int j0=2;

    Basis1D basis(d,d_,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    int j = j0+1;

    Range<int> scalingrange = getRange(basis, j, XBSpline);
    for (int k=scalingrange.firstIndex(); k<=scalingrange.lastIndex(); ++k) {
        stringstream filename;
        filename << "dijkema_bspline_" << k << ".dat";
        ofstream file(filename.str().c_str());
        T l1=basis.mra.phi.support(j,k).l1;
        T l2=basis.mra.phi.support(j,k).l2;
        for (T x=l1; x<=l2; x+=pow2i<T>(-8)) {
            file << x << " " << basis.mra.phi(x,j,k,0) << " " << basis.mra.phi(x,j,k,1) <<  endl;
        }
        file.close();
        getchar();
    }

    Range<int> waveletrange = getRange(basis, j, XWavelet);
    for(int k=waveletrange.firstIndex(); k<=waveletrange.lastIndex(); ++k) {
        stringstream filename;
        filename << "dijkema_wavelet_" << k << ".dat";
        ofstream file(filename.str().c_str());
        T l1=basis.psi.support(j,k).l1;
        T l2=basis.psi.support(j,k).l2;
        for (T x=l1; x<=l2; x+=pow2i<T>(-9)) {
            file << x << " " << basis.psi(x,j,k,0) << " " << basis.psi(x,j,k,1)  << endl;
        }
        file.close();
        getchar();
    }


    return 0;
}


Range<int>
getRange(const Basis1D& basis, int j, XType type)
{
    if (type==XBSpline) {
        return basis.mra.rangeI(j);
    }
    else {
        return basis.rangeJ(j);
    }
}

