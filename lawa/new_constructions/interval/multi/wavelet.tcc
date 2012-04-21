#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC 1

#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
Wavelet<T,Orthogonal,Interval,Multi>::Wavelet(const Basis<T,Orthogonal,Interval,Multi> &_basis)
    : basis(_basis), d(_basis.d), vanishingMoments(_basis.d)
{
}
    
template <typename T>
Wavelet<T,Orthogonal,Interval,Multi>::~Wavelet()
{
}

template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::operator()(T x, int j, Integer k, unsigned short deriv) const
{
    k -= 1;
    
    // left boundary
    if (k<basis._numLeftParts) {
        return pow2ih<T>(2*j*deriv+j) * basis._leftEvaluator[k](pow2i<T>(j)*x, deriv);
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        Integer shift = iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
        return pow2ih<T>(2*j*deriv+j) * 
        basis._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
    
    // right part
    int type  = (int)(k - (basis.cardJ(j)-1 - basis._numRightParts + 1));
    Integer shift = iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
    return pow2ih<T>(2*j*deriv+j) * basis._rightEvaluator[type](pow2i<T>(j)*x-shift, deriv);
}
    
template <typename T>
Support<T>
Wavelet<T,Orthogonal,Interval,Multi>::support(int j, Integer k) const
{
    k -= 1;
    
    // left boundary
    if (k<basis._numLeftParts) {
        return pow2i<T>(-j) * basis._leftSupport[k];
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        int type = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        Integer shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
        return pow2i<T>(-j) * (basis._innerSupport[type]+shift);
    }
    
    // right part
    int type  = (int)(k - (basis.cardJ(j) -1 - basis._numRightParts + 1));
    Integer shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
    return pow2i<T>(-j) * (basis._rightSupport[type]+shift);
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Orthogonal,Interval,Multi>::singularSupport(int j, Integer k) const
{
    k -= 1;    
    // left boundary
    if (k<basis._numLeftParts) {
        return pow2i<T>(-j) * basis._leftSingularSupport[k];
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        Integer shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
        DenseVector<Array<T> > result = basis._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
    
    // right part
    int type  = (int)(k - (basis.cardJ(j)-1 - basis._numRightParts + 1));
    Integer shift = iceil<T>((k+1. - basis._numLeftParts)/basis._numInnerParts);
    DenseVector<Array<T> > result = basis._rightSingularSupport[type];
    result += shift;
    return pow2i<T>(-j) * result;
}
    
template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}
    
} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC
