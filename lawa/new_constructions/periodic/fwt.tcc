#ifndef LAWA_CONSTRUCTIONS_PERIODIC_FWT_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_FWT_TCC 1

#include <cassert>
#include <lawa/aux/integer.h>

namespace lawa {

template <typename X, typename Y>
void
decompose(const DenseVector<X> &x, 
          const Basis<typename X::ElementType,Dual,Periodic,CDF> &basis_, int j,
          DenseVector<Y> &y)
{
    assert(j>=basis_.j0);
    typedef typename X::ElementType T;

    DenseVector<Array<T> > sf, w;
    sf = transpose(basis_.mra_.M0_) * x;
    w  = transpose(basis_.M1_) * x;
    concatenate(sf,w,y,y.firstIndex());
}

template <typename X, typename Y>
void
reconstruct(const DenseVector<X> &x, 
            const Basis<typename X::ElementType,Primal,Periodic,CDF> &basis, int j,
            DenseVector<Y> &y)
{
    assert(j>=basis.j0);
    assert(x.length()%2==0);
    typedef typename X::ElementType T;
    
    int middle = x.firstIndex() + x.length()/2 - 1;
    y.engine().resize(x.range()) || y.engine().fill();
    y  = basis.mra.M0*x(_(x.firstIndex(),middle));
    y += basis.M1*x(_(middle+1,x.lastIndex()));
}

template <typename X, typename Y>
void
fwt(const DenseVector<X> &x, 
    const Basis<typename X::ElementType,Dual,Periodic,CDF> &basis_, int j,
    DenseVector<Y> &y)
{
    assert(j>=basis_.j0);

    typedef typename DenseVector<X>::ElementType T;
    int n = x.length();
    y.engine().resize(n,x.firstIndex()) || y.engine().fill();
    
    DenseVector<Array<T> > tmp = x;
    n /= 2;
    int from = x.firstIndex(), middle = from+n-1, to = x.lastIndex();
    for (int l=j; l>=basis_.j0; --l) {
        n /=2;
        typename DenseVector<Y>::View yview = y(_(from,to));
        decompose(tmp, basis_, l, yview);
        to = middle;
        middle = from + n - 1;
        tmp = y(_(from,to));
    }
}

template <typename X, typename Y>
void
ifwt(const DenseVector<X> &x, 
     const Basis<typename X::ElementType,Primal,Periodic,CDF> &basis, int j,
     DenseVector<Y> &y)
{
    assert(j>=basis.j0);
    typedef typename X::ElementType T;
    
    int n = x.length() / pow2i<T>(j-basis.j0);
    y = x;
    for (int l=basis.j0; l<=j; ++l) {
        Range<Integer> r = _(y.firstIndex(),y.firstIndex()+n-1);
        typename DenseVector<Y>::View yview = y(r);
        DenseVector<Array<T> > z = y(r);
        reconstruct(z, basis, l, yview);
        n *= 2;
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_FWT_TCC
