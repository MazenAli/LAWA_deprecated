#ifndef LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_TCC 1

#include <cassert>
#include <lawa/flensforlawa.h>

namespace flens {

template <typename T>
template <FunctionSide Side>
RefinementMatrix<T,R,CDF>::RefinementMatrix(const BSpline<T,Side,R,CDF> &spline)
      : band(rsqrt2<T>() * spline.mask())
{
}

template <typename T>
template <FunctionSide Side>
RefinementMatrix<T,R,CDF>::RefinementMatrix(const Wavelet<T,Side,R,CDF> &wavelet)
      : band(rsqrt2<T>() * wavelet.mask())
{
}


//------------------------------------------------------------------------------

/*
template <typename X, typename Y>
void
mv(cxxblas::Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,R,CDF> &A,
   const DenseVector<X> &x, typename X::ElementType beta, DenseVector<Y> &y)
{
    assert(0);
    typedef typename X::ElementType T;

    assert(alpha==1.);
    assert(x.engine().stride()==1);

    const DenseVector<Array<T> > &a = A.band;
    int lx = x.length();
    int la = a.length();
    
    if (transA==cxxblas::NoTrans) {
        if (beta==0) {
            y.engine().resize(2*lx,x.firstIndex()-a.firstIndex()) || y.engine().fill();
        } else {
            assert(y.length()==2*lx);
            y.engine().changeIndexBase(x.firstIndex()-a.firstIndex());
        }
        const T *xp = x.engine().data();
        for (int k=0; k<lx; ++k, ++xp) {
            int mMin = a.firstIndex() + 2*k;
            const T *abegin = a.engine().data();
                  T *ybegin = y.engine().data();
            cxxblas::axpy(la, *xp, abegin, 1, ybegin+mMin, 1);
        }
    } else { // (transA==cxxblas::Trans)
        if (beta==0) {
            y.engine().resize(lx/2,x.firstIndex()+a.firstIndex()) || y.engine().fill();
        } else {
            assert(y.length()==lx/2);
            y.engine().changeIndexBase(x.firstIndex()+a.firstIndex());
        }
        T *iter = y.engine().data();
        for (int m=0; m<lx/2; ++m, ++iter) {
            int kMin = 2*m;
            const T *abegin = a.engine().data();
            T dotValue;
            cxxblas::dot(la, abegin, 1, x.engine().data()+kMin, 1, dotValue);
            *iter += dotValue;
        }
    }
}
*/

} // namespace flens

#endif // LAWA_CONSTRUCTIONS_REALLINE_REFINEMENTMATRIX_TCC
