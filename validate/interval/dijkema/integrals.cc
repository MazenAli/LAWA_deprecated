#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef GeMatrix<FullStorage<T,cxxblas::ColMajor> > FullColMatrix;

const Construction Cons = Dijkema;


int
main(int argc, char *argv[])
{
    if (argc!=5) {
        cerr << "usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(-1);
    }
    int d  = atoi(argv[1]), d_ = atoi(argv[2]),
        j0 = atoi(argv[3]), J  = atoi(argv[4]);

    Basis<T,Primal,Interval,Cons> basis(d,d_,j0);
    const MRA<T,Primal,Interval,Primbs> &mra = basis.mra;

    FullColMatrix A(mra.rangeI(J),mra.rangeI(J));
    FullColMatrix Ad(mra.rangeI(J),mra.rangeI(J));

    cout << "SF*SF begin ------------------" << endl;
    BSpline<T,Primal,Interval,Primbs> phi1(mra,0), d_phi1(mra,1),
                                      phi2(mra,0), d_phi2(mra,1);
    Integral<T,Gauss,
             BSpline<T,Primal,Interval,Primbs>,
             BSpline<T,Primal,Interval,Primbs> > integralsfsf(phi1, phi2),
                                              dd_integralsfsf(d_phi1, d_phi2);
                                              
    for (int k1=mra.rangeI(j0).firstIndex(); k1<=mra.rangeI(j0).lastIndex(); ++k1) {
        for (int k2=mra.rangeI(j0).firstIndex(); k2<=mra.rangeI(j0).lastIndex(); ++k2) {
            Ad(k2,k1) = dd_integralsfsf(j0,k1,j0,k2);
            A(k2,k1)  = integralsfsf(j0,k1,j0,k2);
	}
    }
    cout << "SF*SF end --------------------" << endl;
    
    cout << "W*SF begin -------------------" << endl;
    Wavelet<T,Primal,Interval,Cons> psi1(basis,0), d_psi1(basis,1);
    Integral<T,Gauss,
             Wavelet<T,Primal,Interval,Cons>,
             BSpline<T,Primal,Interval,Primbs> > integralwsf(psi1,   phi1),
                                            dd_integralwsf(d_psi1, d_phi1);
    for (int k1=mra.rangeI(j0).firstIndex(); k1<=mra.rangeI(j0).lastIndex(); ++k1) {
        for (int j2=j0; j2 <= J-1; ++j2) {
            for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                Ad(mra.rangeI(j2).lastIndex() + k2, k1) = dd_integralwsf(j2,k2,j0,k1);
                A(mra.rangeI(j2).lastIndex() + k2, k1) =    integralwsf(j2,k2,j0,k1);

            }
        }
    }
    cout << "W*SF end ---------------------" << endl;

    cout << "SF*W begin -------------------" << endl;
    for (int k2=mra.rangeI(j0).firstIndex(); k2<=mra.rangeI(j0).lastIndex(); ++k2) {
        for (int j1=j0; j1<=J-1; ++j1) {
            for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                Ad(k2, mra.rangeI(j1).lastIndex() + k1) = dd_integralwsf(j1,k1,j0,k2);
                A(k2, mra.rangeI(j1).lastIndex() + k1) =    integralwsf(j1,k1,j0,k2);
            }
        }
    }
    cout << "SF*W end ---------------------" << endl;

    cout << "W*W begin --------------------" << endl;
    Wavelet<T,Primal,Interval,Cons> psi2(basis,0), d_psi2(basis,1);
    Integral<T,Gauss,
             Wavelet<T,Primal,Interval,Cons>,
             Wavelet<T,Primal,Interval,Cons> > integralww(psi1,   psi2),
                                            dd_integralww(d_psi1, d_psi2);
    for (int j1=j0; j1<=J-1; ++j1) {
        cout << "j1 = " << j1 << ", j2 = " << flush;
        for (int j2=j0; j2<=J-1; ++j2) {
            cout << j2 << ", " << flush;
            for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                    Ad(mra.rangeI(j2).lastIndex() + k2,
                       mra.rangeI(j1).lastIndex() + k1) = dd_integralww(j1,k1,j2,k2);
                    A(mra.rangeI(j2).lastIndex() + k2,
                      mra.rangeI(j1).lastIndex() + k1) =    integralww(j1,k1,j2,k2);
                }
            }
        }
        cout << endl;
    }
    cout << "W*W end ----------------------" << endl;

    cout << "Ad = [" << Ad << "];" << endl;
    cout << "A = [" << A << "];" << endl;

    return 0;
}
