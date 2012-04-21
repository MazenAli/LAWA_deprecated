#ifndef LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BSPLINE_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BSPLINE_TCC 1

#include <cassert>
#include <lawa/flensforlawa.h>
#include <lawa/math/ceil.h>
#include <lawa/math/floor.h>
#include <lawa/math/polynomial.h>
#include <lawa/settings/param.h>
#include <lawa/constructions/realline/subdivision.h>

namespace lawa {

using namespace flens;

template <typename T>
BSpline<T,Dual,Periodic,CDF>::BSpline(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1), phiR_(_d, _d_)
{
    assert(d>0);
    assert(d_>=d);
    assert(((d+d_)&1)==0);
}

template <typename T>
BSpline<T,Dual,Periodic,CDF>::BSpline(const MRA<T,Dual,Periodic,CDF> &mra)
    : d(mra.d), d_(mra.d_), mu(d&1), phiR_(d,d_)
{
    assert(d>0);
    assert(d_>=d);
    assert(((d+d_)&1)==0);
}

template <typename T>
BSpline<T,Dual,Periodic,CDF>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Dual,Periodic,CDF>::operator()(T x, int j, Integer k, unsigned short deriv) const
{
    assert(deriv==0);

    if((x < 0.) || (x > 1.)){
        return 0.;
    }

    T val = 0;
    for(int l = ifloor<int>(phiR_.support(j,k).l1); l < iceil<int>(phiR_.support(j,k).l2); ++l){
        val += phiR_(l+x, j, k);
    }
    return val;
}


template <typename T>
PeriodicSupport<T>
BSpline<T,Dual,Periodic,CDF>::support(int j, Integer k) const
{    
    Support<T> suppR = phiR_.support(j,k);
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
BSpline<T,Dual,Periodic,CDF>::mask() const
{
    return phiR_.a_;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BSPLINE_TCC
