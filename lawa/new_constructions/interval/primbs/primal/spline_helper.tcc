#ifndef LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_SPLINE_HELPER_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_SPLINE_HELPER_TCC 1

namespace lawa {
    
template <typename T>
T
w(Integer i, int d, const DenseVector<Array<T> > &knots, T x)
{
    assert(1<=i);
    assert(i<=knots.length()-d+1);
    
    if (x<=knots(i)) {
        return 0.0;
    } else if (x>=knots(i+d-1)) {
        return 1.0;
    } else {
        return (x-knots(i)) / (knots(i+d-1)-knots(i));
    }
}

template <typename T>
GeMatrix<FullStorage<T,cxxblas::ColMajor> >
insertKnot(int d, DenseVector<Array<T> > &knots, T x)
{
    assert(knots.length()-d-1>=1);
    
    GeMatrix<FullStorage<T,cxxblas::ColMajor> > ret(knots.length()-d, 
                                                    knots.length()-d-1);
    for (Integer i=ret.firstCol(); i<=ret.lastCol(); ++i) {
        ret(i,i) = w(i,d+1,knots,x);
        ret(i+1,i) = 1-w(i+1,d+1,knots,x);
    }
    std::list<T> temp;
    for (Integer i=knots.firstIndex(); i<=knots.lastIndex(); ++i) {
        temp.push_back(knots(i));
    }
    temp.push_back(x);
    temp.sort();
    knots.engine().resize(knots.length()+1);
    typename std::list<T>::const_iterator it=temp.begin();
    for (int i=1; it!=temp.end(); ++it, ++i) {
        knots(i) = *it;
    }
    return ret;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_SPLINE_HELPER_TCC
