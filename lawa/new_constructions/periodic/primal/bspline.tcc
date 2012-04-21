#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_TCC 1

#include <cassert>
#include <cmath>
#include <lawa/flensforlawa.h>

#include <lawa/math/math.h>
#include <lawa/constructions/realline/primal/bspline.h>
#include <extensions/extensions.h>

namespace lawa {

using namespace flens;

template <typename T>
BSpline<T,Primal,Periodic,CDF>::BSpline(int _d)
    : d(_d), mu(d&1), phiR(_d)
{
    assert(_d>0);
}

template <typename T>
BSpline<T,Primal,Periodic,CDF>::BSpline(const MRA<T,Primal,Periodic,CDF> &mra)
    : d(mra.d), mu(d&1), phiR(d)
{
    assert(d>0);
}

template <typename T>
BSpline<T,Primal,Periodic,CDF>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Primal,Periodic,CDF>::operator()(T x, int j, Integer k, unsigned short deriv) const
{
    // maximal support: [0,1]
    if((x < 0.) || (x > 1.)){
        return 0.;
    }
    
    // add contributions of original spline on R
    // = 'wrapping' around [0,1]
    T val = 0;
    for(int l = ifloor<int>(phiR.support(j,k).l1); l < iceil<int>(phiR.support(j,k).l2); ++l){
        val += phiR(l+x, j, k, deriv);
    }
    return val;
    
}

template <typename T>
PeriodicSupport<T>
BSpline<T,Primal,Periodic,CDF>::support(int j, Integer k) const
{
    Support<T> suppR = phiR.support(j,k);
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
DenseVector<Array<T> >
BSpline<T,Primal,Periodic,CDF>::singularSupport(int j, Integer k) const
{   
    if((phiR.support(j,k).l1 >= 0) && (phiR.support(j,k).l2 <= 1)){
         return linspace(support(j,k).l1, support(j,k).l2, d+1);
    }
    
    std::list<T> temp;
    DenseVector<Array<T> > singSuppR = linspace(phiR.support(j,k).l1, phiR.support(j,k).l2, d+1);
    temp.push_back(0.);
    temp.push_back(1.);
    for( Integer  i = singSuppR.firstIndex(); i <= singSuppR.lastIndex(); ++i){
        temp.push_back(singSuppR(i) - ifloor<int>(singSuppR(i)));
    }
    temp.sort();
    temp.unique();
    
    DenseVector<Array<T> > singSupp(temp.size());
    int i = 1;
    for (typename std::list<T>::const_iterator it = temp.begin(); it != temp.end(); ++it, ++i) {
        singSupp(i) = *it;
    }
    
    return singSupp;
}

template <typename T>
T
BSpline<T,Primal,Periodic,CDF>::tic(int j) const
{
    return pow2i<T>(-j);
}

template <typename T>
const DenseVector<Array<T> > &
BSpline<T,Primal,Periodic,CDF>::mask() const
{
    return phiR.a;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BSPLINE_TCC
