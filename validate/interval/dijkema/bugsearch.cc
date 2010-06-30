#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;
const Construction Cons = Dijkema;

typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

int
main(int argc, char *argv[])
{
    if (argc!=6) {
        cerr << "usage: " << argv[0] << " d d_ j0 J k" << endl;
        exit(-1);
    }
    int d  = atoi(argv[1]), d_ = atoi(argv[2]),
        j0 = atoi(argv[3]), J  = atoi(argv[4]),
        k  = atoi(argv[5]);

    Basis<T,Primal,Interval,Dijkema> basis(d,d_);
    Basis<T,Dual,Interval,Dijkema> basis_(d,d_);
    basis.enforceBoundaryCondition<DirichletBC>();
    basis_.enforceBoundaryCondition<DirichletBC>();

    BSpline<T,Primal,R,CDF> phi(d);
    BSpline<T,Dual,R,CDF> phi_(d,d_);
    Wavelet<T,Primal,R,CDF> psi(d,d_);
    Wavelet<T,Dual,R,CDF> psi_(d,d_);

    std::cerr << basis.mra.M0.leftband << phi.a << std::endl;
    std::cerr << basis_.mra_.M0_.leftband << phi_.a_ << std::endl;
    std::cerr << basis.M1.leftband << psi.b << std::endl;
    std::cerr << basis_.M1_.leftband << psi_.b_ << std::endl;

    GeMatrix<FullStorage<T,ColMajor> > DM1, DM1_, I;
    densify(cxxblas::NoTrans, basis.M1, DM1);
    densify(cxxblas::NoTrans, basis_.M1_, DM1_);
    blas::mm(Trans,NoTrans,1.,DM1,DM1_,0.,I);
    std::cerr << I << std::endl;

    DenseVectorT coeffs(basis.mra.rangeI(J)), dec;
    coeffs(k) = 1.;
    T x = 0.0, h = 1./1024.;
    for (T x=0.0; x<=1.; x+=h) {
        cout << evaluate(basis.mra,J,coeffs,x,0) << endl;
    }
    cout << endl << endl;
    std::cerr << coeffs << std::endl;
        
    fwt(coeffs,basis_,J-1,dec);
    std::cerr << dec << std::endl;
    dec(dec.firstIndex()+15) = 0.;
    std::cerr << dec << std::endl;
    x = 0.0, h = 1./1024.;
    for (T x=0.0; x<=1.; x+=h) {
        cout << evaluate(basis,J,dec,x,0) << endl;
    }
    cout << endl << endl;
    ifwt(dec,basis,J-1,coeffs);
    std::cerr << coeffs << std::endl;

    x = 0.0, h = 1./1024.;
    for (T x=0.0; x<=1.; x+=h) {
        cout << evaluate(basis.mra,J,coeffs,x,0) << endl;
    }
    cout << endl << endl;

    return 0;
}
