#ifndef LAWA_CONSTRUCTIONS_REALLINE_DUAL_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_DUAL_WAVELET_TCC 1

#include <cassert>

namespace lawa {

using namespace flens;

template <typename T>
Wavelet<T,Dual,R,CDF>::Wavelet(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1), l1_((2-(d+d_))/2), l2_((d+d_)/2),
      b_(mask(d,d_)), phi(d), phi_(d,d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);
}

template <typename T>
Wavelet<T,Dual,R,CDF>::Wavelet(const BSpline<T,Primal,R,CDF> &_phi,
                               const BSpline<T,Dual,R,CDF> &_phi_)
    : d(_phi.d), d_(_phi_.d_), mu(d&1), l1_((2-(d+d_))/2), l2_((d+d_)/2),
      b_(mask(d,d_)), phi(_phi), phi_(_phi_)
{
}

template <typename T>
T
Wavelet<T,Dual,R,CDF>::operator()(T x, int j, Integer k, unsigned short deriv) const
{
    assert(deriv==0);
    
    T ret = T(0);
    x = pow2i<T>(j)*x-k;
    for (int i=b_.firstIndex(); i<=b_.lastIndex(); ++i) {
        ret += b_(i)*phi_(2*x-i, 0, 0);
    }
    return pow2ih<T>(j) * ret;
}

template <typename T>
Support<T>
Wavelet<T,Dual,R,CDF>::support(int j, Integer k) const
{
    return pow2i<T>(-j) * Support<T>(l1_+k, l2_+k);
}

template <typename T>
const DenseVector<Array<T> > &
Wavelet<T,Dual,R,CDF>::mask() const
{
    return b_;
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Dual,R,CDF>::mask(int d, int d_)
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

#endif // LAWA_CONSTRUCTIONS_REALLINE_DUAL_WAVELET_TCC
