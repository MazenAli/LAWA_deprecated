#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

int
main(int argc, char *argv[])
{
    if (argc!=5) {
        cerr << "usage: " << argv[0] << " d d_ j deriv" << endl;
        exit(-1);
    }
    int d     = atoi(argv[1]);
    int d_    = atoi(argv[2]);
    int j     = atoi(argv[3]);
    int deriv = atoi(argv[4]);

    Basis<T,Primal,Interval,Dijkema> basis(d,d_,j);
    BSpline<T,Primal,Interval,Primbs> phi(basis.mra,deriv);
    Wavelet<T,Primal,Interval,Dijkema> psi(basis,deriv);

    std::cerr << "--- checking supports of B-splines ... ";
    // --- primal B-splines
    const MRA<T,Primal,Interval,Primbs> &mra = basis.mra;
    // ... supports are large enough
    for (int k=mra.rangeI(j).firstIndex(); k<=mra.rangeI(j).lastIndex(); ++k) {
        for (T x=0.; x<=1.; x+=1./1024) {
            if ((x<phi.support(j,k).l1 || x>phi.support(j,k).l2) 
             && (phi(x,j,k)!=0.)) {
                std::cerr << "non-zero outside support: "<< k << ": " 
                          << x << " " << phi.support(j,k) << std::endl;                
                assert(0);
            }
        }
    }
    // ... supports are not too large
    for (int k=mra.rangeI(j).firstIndex(); k<=mra.rangeI(j).lastIndex(); ++k) {
        T x;
        for (x=0.; phi(x,j,k)==0; x+=1./4096) ;
        if (fabs(phi.support(j,k).l1-x)>pow2i<T>(-j-2)) {
            std::cerr << "support too large (left) : " << k << ": " << x << " " 
                      << phi.support(j,k) << std::endl;
            assert(fabs(phi.support(j,k).l1-x)<=pow2i<T>(-j-2));
        }
        for (x = 1.; phi(x,j,k)==0; x-=1./4096) ;
        if (fabs(phi.support(j,k).l2-x)>pow2i<T>(-j-2)) {
            std::cerr << "support too large (right): " << k << ": " << x << " " 
                      << phi.support(j,k) << std::endl;
            assert(fabs(phi.support(j,k).l2-x)<=pow2i<T>(-j-2));
        }
    }
    std::cerr << "done. ---" << std::endl;
    std::cerr << "--- checking supports of wavelets ... ";
    //--- primal wavelets
    // ... supports are large enough
    for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
        for (T x = 0.; x<=1.; x+=1./1024) {
            if ((x<psi.support(j,k).l1 || x>psi.support(j,k).l2) 
             && (psi(x,j,k)!=0.)) {
                std::cerr << "non-zero outside support: "<< k << ": " 
                          << x << " " << psi.support(j,k) << std::endl;                
                assert(0);
            }
        }
    }
    // ... supports are not too large
    for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
        T x;
        for (x = 0.; psi(x,j,k)==0; x+=1./4096) ;
        if (fabs(psi.support(j,k).l1-x)>pow2i<T>(-j-2)) {
            std::cerr << "support too large (left) : " << k << ": " << x << " " 
                      << psi.support(j,k) << std::endl;
            assert(fabs(psi.support(j,k).l1-x)<=pow2i<T>(-j-2));
        }
        for (x = 1.; psi(x,j,k)==0; x-=1./4096) ;
        if (fabs(psi.support(j,k).l2-x)>pow2i<T>(-j-2)) {
            std::cerr << "support too large (right): " << k << ": " << x << " " 
                      << psi.support(j,k) << std::endl;
            assert(fabs(psi.support(j,k).l2-x)<=pow2i<T>(-j-2));
        }
    }
    std::cerr << "done. ---" << std::endl;

    Basis<T,Dual,Interval,Dijkema> basis_(d,d_,j);
    BSpline<T,Dual,Interval,Dijkema> phi_(basis_.mra_);
//    Wavelet<T,Dual,Interval,Dijkema> psi_(basis_,deriv);

    std::cerr << "--- checking supports of B-splines ... ";
    // --- dual B-splines
    const MRA<T,Dual,Interval,Dijkema> &mra_ = basis_.mra_;
    // ... supports are large enough
    for (int k=mra_.rangeI_(j).firstIndex(); k<=mra_.rangeI_(j).lastIndex(); ++k) {
        for (T x=0.; x<=1.; x+=1./1024) {
            if ((x<phi_.support(j,k).l1 || x>phi_.support(j,k).l2) 
             && (phi_(x,j,k)!=0.)) {
                std::cerr << "non-zero outside support: "<< k << ": " 
                          << x << " " << phi_.support(j,k) << std::endl;                
                assert(0);
            }
        }
    }
    // ... supports are not too large
    for (int k=mra_.rangeI_(j).firstIndex(); k<=mra_.rangeI_(j).lastIndex(); ++k) {
        T x;
        for (x=0.; phi_(x,j,k)==0; x+=1./4096) ;
        if (fabs(phi_.support(j,k).l1-x)>pow2i<T>(-j-2)) {
            std::cerr << "support too large (left) : " << k << ": " << x << " " 
                      << phi_.support(j,k) << std::endl;
            assert(fabs(phi_.support(j,k).l1-x)<=pow2i<T>(-j-2));
        }
        for (x = 1.; phi_(x,j,k)==0; x-=1./4096) ;
        if (fabs(phi_.support(j,k).l2-x)>pow2i<T>(-j-2)) {
            std::cerr << "support too large (right): " << k << ": " << x << " " 
                      << phi_.support(j,k) << std::endl;
            assert(fabs(phi_.support(j,k).l2-x)<=pow2i<T>(-j-2));
        }
    }
    std::cerr << "done. ---" << std::endl;
/*    std::cerr << "--- checking supports of wavelets ... ";
    //--- dual wavelets
    // ... supports are large enough
    for (int k=basis_.rangeJ(j).firstIndex(); k<=basis_.rangeJ(j).lastIndex(); ++k) {
        for (T x = 0.; x<=1.; x+=1./1024) {
            if ((x<psi_.support(j,k).l1 || x>psi_.support(j,k).l2) 
             && (psi_(x,j,k)!=0.)) {
                std::cerr << "non-zero outside support: "<< k << ": " 
                          << x << " " << psi_.support(j,k) << std::endl;                
                assert(0);
            }
        }
    }
    // ... supports are not too large
    for (int k=basis_.rangeJ(j).firstIndex(); k<=basis_.rangeJ(j).lastIndex(); ++k) {
        T x;
        for (x = 0.; psi_(x,j,k)==0; x+=1./4096) ;
        if (fabs(psi_.support(j,k).l1-x)>pow2i<T>(-j-2)) {
            std::cerr << "support too large (left) : " << k << ": " << x << " " 
                      << psi_.support(j,k) << std::endl;
            assert(fabs(psi_.support(j,k).l1-x)<=pow2i<T>(-j-2));
        }
        for (x = 1.; psi_(x,j,k)==0; x-=1./4096) ;
        if (fabs(psi_.support(j,k).l2-x)>pow2i<T>(-j-2)) {
            std::cerr << "support too large (right): " << k << ": " << x << " " 
                      << psi_.support(j,k) << std::endl;
            assert(fabs(psi_.support(j,k).l2-x)<=pow2i<T>(-j-2));
        }
    }
    std::cerr << "done. ---" << std::endl;
*/
    return 0;
}
