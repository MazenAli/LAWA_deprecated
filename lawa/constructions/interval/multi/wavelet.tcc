#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC 1

#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
Wavelet<T,Orthogonal,Interval,Multi>::Wavelet(const Basis<T,Orthogonal,Interval,Multi> &_basis)
    : basis(_basis), d(_basis.d), vanishingMoments(_basis.d)
{
    switch (d) {
        case 2:
            initialticsize = pow2i<T>(-3);
            break;

        case 3:
            initialticsize = pow2i<T>(-4);
            break;

        case 4:
            initialticsize = pow2i<T>(-4);
            break;

        default: std::cerr << "BSpline<T,Orthogonal,Interval,Multi> not yet realized"
                    " for d = " << d << ". Stopping." << std::endl;
                    exit(-1);
    }
}
    
template <typename T>
Wavelet<T,Orthogonal,Interval,Multi>::~Wavelet()
{
}

template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    k -= 1;
    // left boundary
    if (k<basis._numLeftParts) {
        //std::cerr << " Left boundary: k = " << k << std::endl;
        return pow2ih<T>(2*j*deriv+j) * basis._leftEvaluator[k](pow2i<T>(j)*x, deriv);
    }
    
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = lawa::iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
        //std::cerr << "  k = " << k << " : type = " << type << ", shift = " << shift << std::endl;
        return pow2ih<T>(2*j*deriv+j) * 
        basis._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
    
    // right part
    int type  = (int)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
    return pow2ih<T>(2*j*deriv+j) * basis._rightEvaluator[type](pow2i<T>(j)*x-shift, deriv);
}
    
template <typename T>
Support<T>
Wavelet<T,Orthogonal,Interval,Multi>::support(int j, long k) const
{
    k -= 1;
    
    // left boundary
    if (k<basis._numLeftParts) {
        return pow2i<T>(-j) * basis._leftSupport[k];
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        int type = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
        return pow2i<T>(-j) * (basis._innerSupport[type]+shift);
    }
    
    // right part
    int type  = (int)(k - (basis.cardJ(j) -1 - basis._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
    return pow2i<T>(-j) * (basis._rightSupport[type]+shift);
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Orthogonal,Interval,Multi>::singularSupport(int j, long k) const
{
    k -= 1;    
    // left boundary
    if (k<basis._numLeftParts) {
        return pow2i<T>(-j) * basis._leftSingularSupport[k];
    }
    
    // inner part
    if (k<basis.cardJL()+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = iceil<T>((k+1.-basis._numLeftParts)/basis._numInnerParts);
        DenseVector<Array<T> > result = basis._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
    
    // right part
    int type  = (int)(k - (basis.cardJ(j)-1 - basis._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1. - basis._numLeftParts)/basis._numInnerParts);
    DenseVector<Array<T> > result = basis._rightSingularSupport[type];
    result += shift;
    return pow2i<T>(-j) * result;
}
    
template <typename T>
T
Wavelet<T,Orthogonal,Interval,Multi>::tic(int j) const
{
    //todo: Critical! Yields totally wrong results if too small, increases run-time when too large!
    //return pow2i<T>(-(j+4));
    return initialticsize*pow2i<T>(-j);
}

template <typename T>
DenseVector<Array<long double> > *
Wavelet<T,Orthogonal,Interval,Multi>::getRefinement(int j, long k,
                                                    int &refinement_j, long &refinement_k_first) const
{
    k -= 1;
    refinement_j = j + basis._addRefinementLevel;
    // left boundary
    if (k<basis._numLeftParts) {
        refinement_k_first = basis._leftOffsets[k];
        return &(basis._leftRefCoeffs[k]);
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = (long)lawa::iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
        //refinement_k_first = pow2i<long>(basis._addRefinementLevel)*shift+basis._innerOffsets[type];
        refinement_k_first = basis._shiftFactor*shift+basis._innerOffsets[type];
        return &(basis._innerRefCoeffs[type]);
    }
    // right part
    int type  = (int)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    long shift = (long)pow2i<long>(j)-1;
    //refinement_k_first = pow2i<long>(basis._addRefinementLevel)*shift+basis._rightOffsets[type];
    refinement_k_first = basis._shiftFactor*shift+basis._rightOffsets[type];
    return &(basis._rightRefCoeffs[type]);
}

template <typename T>
void
Wavelet<T,Orthogonal,Interval,Multi>::getRefinementNeighbors(int j, long k, int &refinement_j,
                                                             long &refinement_k_first,
                                                             long &refinement_k_last) const
{
    k -= 1;
    refinement_j = j + basis._addRefinementLevel;
    // left boundary
    if (k<basis._numLeftParts) {
        refinement_k_first = basis._leftOffsets[k];
        refinement_k_last  = refinement_k_first + basis._leftRefCoeffs[k].lastIndex();
        return;
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        int type  = (int)((k-basis._numLeftParts) % basis._numInnerParts);
        long shift = (long)lawa::iceil<T>((k+1.-basis._numLeftParts)/(double)basis._numInnerParts);
        refinement_k_first = basis._shiftFactor*shift+basis._innerOffsets[type];
        refinement_k_last  = refinement_k_first + basis._innerRefCoeffs[type].lastIndex();
        return;
    }
    // right part
    int type  = (int)(k+1 - (basis.cardJ(j) - basis._numRightParts + 1));
    long shift = (long)pow2i<long>(j)-1;
    refinement_k_first = basis._shiftFactor*shift+basis._rightOffsets[type];
    refinement_k_last  = refinement_k_first + basis._rightRefCoeffs[type].lastIndex();
    return;
}

template <typename T>
void
Wavelet<T,Orthogonal,Interval,Multi>::getRefinedNeighbors(int refinement_j, long refinement_k,
                                                          int &j, long &k_first, long &k_last) const
{
    j = refinement_j - basis._addRefinementLevel+1;
    Support<T> supp = basis.refinementbasis.mra.phi.support(refinement_j,refinement_k);
    T a = supp.l1, b = supp.l2;

    if (a==0.L) {
        k_first = 1;
        k_last = basis._numLeftParts + basis._numInnerParts;
        k_last  = std::min(k_last, (long)basis.rangeJR(j).lastIndex());
        return;
    }
    if (0<a && b<1.L) {
        long k_tilde = (long)std::floor(pow2i<T>(j)*a);
        k_tilde += 1;
        T tmp = k_tilde*pow2i<T>(-j);
        if (a<tmp && tmp<b) {
            k_first  = (k_tilde-2)*basis._numInnerParts;
            k_last   = (k_tilde+2)*basis._numInnerParts;
        }
        else {
            k_first  = (k_tilde-2)*basis._numInnerParts;
            k_last   = (k_tilde+1)*basis._numInnerParts;
        }
        long k_first_min = basis.rangeJL(j).firstIndex();
        long k_last_max  = basis.rangeJR(j).lastIndex();
        if (k_first < k_first_min) k_first = k_first_min;
        if (k_last  > k_last_max)  k_last  = k_last_max;
        return;
    }
    k_last   = basis.rangeJ(j).lastIndex();
    k_first  = k_last - (basis._numRightParts + basis._numInnerParts) + 1;
    k_first  = std::max(1L, k_first);


    return;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_WAVELET_TCC
