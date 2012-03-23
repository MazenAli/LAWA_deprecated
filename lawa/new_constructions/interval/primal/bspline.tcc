#ifndef LAWA_CONSTRUCTIONS_INTERVAL_PRIMAL_BSPLINE_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_PRIMAL_BSPLINE_TCC 1

#include <cassert>
#include <algorithm>
#include <lawa/math/math.h>

namespace lawa {

template <typename T>
    T
    _evaluateUnitBSpline(int d, T x, int j,  Integer k, unsigned short deriv);

//------------------------------------------------------------------------------

template <typename T, Construction Cons>
BSpline<T,Primal,Interval,Cons>::BSpline(const MRA<T,Primal,Interval,Cons> &_mra)
    : mra(_mra)
{
}

template <typename T, Construction Cons>
T
BSpline<T,Primal,Interval,Cons>::operator()(T x, int j,  Integer k, unsigned short deriv) const
{
    assert(j>=mra.j0);
    assert(k>=mra.rangeI(j).firstIndex());
    assert(k<=mra.rangeI(j).lastIndex());
    return _evaluateUnitBSpline(mra.d, x, j, k, deriv);
}

template <typename T, Construction Cons>
Support<T>
BSpline<T,Primal,Interval,Cons>::support(int j,  Integer k) const
{
    assert(j>=mra.j0);
    assert(k>=mra.rangeI(j).firstIndex());
    assert(k<=mra.rangeI(j).lastIndex());
    return pow2i<T>(-j) * Support<T>(std::max(0L,k-mra.d),
                                     std::min(k,pow2i< Integer>(j)));
}

template <typename T, Construction Cons>
DenseVector<Array<T> >
BSpline<T,Primal,Interval,Cons>::singularSupport(int j,  Integer k) const
{
    const  Integer tics = (k<mra.d) ? k+1 : (k>pow2i<Integer>(j)) ? pow2i<Integer>(j)+mra.d-1-k+2 : mra.d+1;
    return linspace(pow2i<T>(-j) * std::max(0L,k-mra.d),
                    pow2i<T>(-j) * std::min(k,pow2i<Integer>(j)),
                    tics);
}

template <typename T, Construction Cons>
T
BSpline<T,Primal,Interval,Cons>::tic(int j) const
{
    return pow2i<T>(-j);
}

//--- evaluate B-spline --------------------------------------------------------

template <typename T>
T
_evaluateUnitBSpline(int d, T x, int j,  Integer k, unsigned short deriv)
{
    assert(x>=0.0);
    assert(x<=1.0);
    
    if (deriv>=d) {
        return 0;
    }
    
    if (deriv==0) {
        // "if" needed for calculation of derivatives
        //  (otherwise assertion would be correct).
        if ((k>pow2i<Integer>(j)+d-1) || (k<1)) {
            return 0;
        }

         Integer twoj = pow2i< Integer>(j); 
        x *= twoj;
         Integer pos = ifloor<Integer>(x) - (x==twoj);
        // we are not inside the support.
        if ((pos<k-d) || (pos>k-1)) {
            return 0;
        }

        Array<T> values(d,0);
        // initialize correct 'slot'.
        if ((pos<=k) && (pos>=k-d)) {
            values(pos-(k-d)) = 1.;
        }

        // utilizing left multiplicities.
        if (k<d) {
            for (int m=2; m<=d; ++m) {
                for (int i=0; i<=d-m; ++i) {
                    if (m+i+k-1>d) {
                        values(i) = (x-std::max(0L,i+k-d))*values(i) / (m+i+k-1-d-std::max(k+i-d,0L));                    
                    }
                    if (m+i+k>d) {
                        values(i) += (m+i+k-d-x)*values(i+1) / (m+i+k-d-std::max(k+i+1-d,0L));                    
                    }
                }
            }
            return pow2ih<T>(j)*values(0);
        }

        // utilizing right multiplicities.
         Integer t = twoj;
        if (k>t) {
            for (int m=2; m<=d; ++m) {
                for (int i=0; i<=d-m; ++i) {
                    if (k+i-d<t) {
                        values(i) = (x-(i+k-d))*values(i) / (std::min(m+i+k-1-d,t)-(k+i-d));                    
                    }
                    if (k+i+1-d<t) {
                        values(i) += (std::min(t,m+i+k-d)-x)*values(i+1) / (std::min(t,m+i+k-d)-(k+i+1-d));                    
                    }
                }
            }
            return pow2ih<T>(j)*values(0);
        }

        // 'inner' B-Splines.
        t = k-d;
        for (int m=2; m<=d; ++m) {
            for (int i=0; i<=d-m; ++i) {
                values(i) =  ((x-(t+i))*values(i) + ((t+m+i)-x)*values(i+1))/(m-1);
            }
        }
        return pow2ih<T>(j)*values(0);
    } else {
        assert(k>=1);
        assert(k<=pow2i<T>(j)+d-1);

        T value = 0.;

        T twomj = pow2i<T>(-j);
        T a=twomj*(k-1), b=twomj*(k), c=twomj*(k-d), e=twomj*(k-d+1);
        if (k<=d) {
            c = 0;
        }
        if (k+d>=pow2i<T>(j)+d) {
            b = 1.;
        }
        if (k+d-1<=d) {
            a = 0.;
        }
        if (k+d-1>=pow2i<T>(j)+d) {
            a = 1.;
        }
        if (k+1<=d) {
            e = 0.;
        }

        // remark: the index k is shifted by -1 since for d-1 the implicit knot
        //         vector is implicitely shifted to the left by one.
        if (a!=c) {
            value  = _evaluateUnitBSpline(d-1,x,j,k-1,deriv-1)   / (a-c);                    
        }
        if (b!=e) {
            value -= _evaluateUnitBSpline(d-1,x,j,k,deriv-1) / (b-e);                    
        }

        return (d-1)*value;
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_PRIMAL_BSPLINE_H
