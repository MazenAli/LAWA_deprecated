#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_WAVELET_TCC 1

#include <cassert>

namespace lawa {

using namespace flens;

template <typename T>
Wavelet<T,Primal,R,CDF>::Wavelet(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1), l1((2-d-d_)/2), l2((d+d_)/2),
      vanishingMoments(_d_), b(mask(d,d_)),
      phi(d), phi_(d,d_)
{
    if (d==2 && d_==2) {
        singularPts.engine().resize(5l);
        singularPts = -1., 0., 0.5, 1., 2.;
    }
    else if (d==2 && d_==4) {
        singularPts.engine().resize(7l);
        singularPts = -2., -1., 0., 0.5, 1., 2., 3.;
    }
    else if (d==3 && d_==3) {
        singularPts.engine().resize(7l);
        singularPts = -2., -1., 0., 0.5, 1., 2., 3.;
    }
    else if (d==3 && d_==5) {
        singularPts.engine().resize(9l);
        singularPts = -3., -2., -1., 0., 0.5, 1., 2., 3., 4.;
    }
    else if (d==3 && d_==7) {
        singularPts.engine().resize(11l);
        singularPts = -4.,-3., -2., -1., 0., 0.5, 1., 2., 3., 4., 5.;
    }
    else {
        std::cout << "Optimized singular points not implemented!" << std::endl;
        singularPts = linspace<T>(l1, l2, 2*(d+d_)-1);
    }
}

template <typename T>
Wavelet<T,Primal,R,CDF>::Wavelet(const BSpline<T,Primal,R,CDF> &_phi,
                                 const BSpline<T,Dual,R,CDF> &_phi_)
    : d(_phi.d), d_(_phi_.d_), mu(d&1), l1((2-d-d_)/2), l2((d+d_)/2),
      vanishingMoments(d_), b(mask(d,d_)),
      phi(_phi), phi_(_phi_)
      
{
    if (d==2 && d_==2) {
        singularPts.engine().resize(5);
        singularPts = -1., 0., 0.5, 1., 2.;
    }
    else if (d==2 && d_==4) {
        singularPts.engine().resize(7);
        singularPts = -2., -1., 0., 0.5, 1., 2., 3.;
    }
    else if (d==3 && d_==3) {
        singularPts.engine().resize(7);
        singularPts = -2., -1., 0., 0.5, 1., 2., 3.;
    }
    else if (d==3 && d_==5) {
        singularPts.engine().resize(9);
        singularPts = -3., -2., -1., 0., 0.5, 1., 2., 3., 4.;
    }
    else if (d==3 && d_==7) {
        singularPts.engine().resize(11);
        singularPts = -4.,-3., -2., -1., 0., 0.5, 1., 2., 3., 4., 5.;
    }
    else {
        std::cout << "Optimized singular points not implemented!" << std::endl;
        singularPts = linspace<T>(l1, l2, 2*(d+d_)-1);
        assert(0);
    }
}

template <typename T>
Wavelet<T,Primal,R,CDF>::Wavelet(const Basis<T,Primal,R,CDF> &basis)
    : d(basis.d), d_(basis.d_), mu(d&1), l1((2-d-d_)/2), l2((d+d_)/2),
      vanishingMoments(basis.d_), b(mask(d,d_)),
      phi(d), phi_(d,d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);
    if (d==2 && d_==2) {
        singularPts.engine().resize(5);
        singularPts = -1., 0., 0.5, 1., 2.;
    }
    else if (d==2 && d_==4) {
        singularPts.engine().resize(7);
        singularPts = -2., -1., 0., 0.5, 1., 2., 3.;
    }
    else if (d==3 && d_==3) {
        singularPts.engine().resize(7);
        singularPts = -2., -1., 0., 0.5, 1., 2., 3.;
    }
    else if (d==3 && d_==5) {
        singularPts.engine().resize(9);
        singularPts = -3., -2., -1., 0., 0.5, 1., 2., 3., 4.;
    }
    else if (d==3 && d_==7) {
        singularPts.engine().resize(11);
        singularPts = -4.,-3., -2., -1., 0., 0.5, 1., 2., 3., 4., 5.;
    }
    else {
        std::cout << "Optimized singular points not implemented!" << std::endl;
        singularPts = linspace<T>(l1, l2, 2*(d+d_)-1);
    }
}

template <typename T>
T
Wavelet<T,Primal,R,CDF>::operator()(T x, int j, Integer k, unsigned short deriv) const
{
    T ret = T(0);
    x = pow2i<T>(j)*x-k;
    for (Integer i=b.firstIndex(); i<=b.lastIndex(); ++i) {
        ret += b(i)*phi(2*x-i, 0, 0, deriv);
    }
    return pow2i<T>(deriv*(j+1)) * pow2ih<T>(j) * ret;
}

template <typename T>
Support<T>
Wavelet<T,Primal,R,CDF>::support(int j, Integer k) const
{
    //required for large negative levels!!
    T scale = pow2i<T>(-j);
    return  Support<T>(scale*l1+scale*k, scale*l2+scale*k);
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Primal,R,CDF>::singularSupport(int j, Integer k) const
{
    return linspace(support(j,k).l1, support(j,k).l2, 2*(d+d_)-1);
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Primal,R,CDF>::optim_singularSupport(int j, Integer k) const
{
    //return linspace(support(j,k).l1, support(j,k).l2, 2*(d+d_)-1);
    flens::DenseVector<flens::Array<T> > x(singularPts.length());
    for (int i=1; i<=x.length(); ++i) {
        x(i) = pow2i<T>(-j)*(singularPts(i)+T(k));
    }
    return x;
}

template <typename T>
T
Wavelet<T,Primal,R,CDF>::tic(int j) const
{
    return pow2i<T>(-(j+1));
}

template <typename T>
const DenseVector<Array<T> > &
Wavelet<T,Primal,R,CDF>::mask() const
{
    return b;
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Primal,R,CDF>::mask(int d, int d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);

    int mu = d & 1;
    BSpline<T,Dual,R,CDF> phi_(d,d_);
    DenseVector<Array<T> > b(_(2-(d+mu)/2-d_, (d-mu)/2+d_));
    for (Integer k=b.firstIndex(); k<=b.lastIndex(); ++k) {
        int sign = (k&1) ? -1 : 1;
        b(k) = sign * phi_.a_(1-k);
    }
    return b;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_WAVELET_TCC
