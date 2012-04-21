#ifndef LAWA_CONSTRUCTIONS_PERIODIC_DUAL_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_DUAL_WAVELET_TCC 1

#include <cassert>

namespace lawa {

using namespace flens;

template <typename T>
Wavelet<T,Dual,Periodic,CDF>::Wavelet(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1),
      psiR_(_d, _d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);
}

template <typename T>
Wavelet<T,Dual,Periodic,CDF>::Wavelet(const BSpline<T,Primal,Periodic,CDF> &_phi,
                                      const BSpline<T,Dual,Periodic,CDF> &_phi_)
    : d(_phi.d), d_(_phi_.d_), mu(d&1), psiR_(_phi, _phi_)
{
}

template <typename T>
Wavelet<T,Dual,Periodic,CDF>::Wavelet(const Basis<T,Dual,Periodic,CDF> &_basis)
    : d(_basis.d), d_(_basis.d_), mu(d&1), psiR_(d,d_)
{
}


template <typename T>
T
Wavelet<T,Dual,Periodic,CDF>::operator()(T x, int j, Integer k, unsigned short deriv) const
{
    assert(deriv==0);
    
    // maximal support: [0,1]
    if((x < 0.) || (x > 1.)){
        return 0.;
    }

    // sum contributions of original spline on R
    // = 'wrapping' around [0,1]
    T val = 0;
    for(int l = ifloor(psiR_.support(j,k).l1); l < iceil(psiR_.support(j,k).l2); ++l){
        val += psiR_(l+x, j, k);
    }
    return val;
}

template <typename T>
PeriodicSupport<T>
Wavelet<T,Dual,Periodic,CDF>::support(int j, Integer k) const
{
    Support<T> suppR = psiR_.support(j,k);
    if(suppR.length() >= 1){
        return PeriodicSupport<T>(0,1);
    }
    if(suppR.l1 < 0){
        return PeriodicSupport<T>(0,1,suppR.l2, suppR.l1 + 1);
    }
    if(suppR.l2 > 1){
        return PeriodicSupport<T>(0,1,suppR.l2 - 1, suppR.l1);
    }
    return PeriodicSupport<T>(suppR.l1, suppR.l2);
}

template <typename T>
const DenseVector<Array<T> > &
Wavelet<T,Dual,Periodic,CDF>::mask() const
{
    return psiR_.b_;
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Dual,Periodic,CDF>::mask(int d, int d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);

    int mu = d & 1;
    BSpline<T,Primal,R,CDF>  phi(d);
    DenseVector<Array<T> > b_(_(1-(d+mu)/2, 1+(d-mu)/2));
    for (int k=b_.firstIndex(); k<=b_.lastIndex(); ++k) {
        int sign = (k&1) ? -1 : 1;
        b_(k) = sign * phi.a(1-k);
    }
    return b_;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_DUAL_WAVELET_TCC
