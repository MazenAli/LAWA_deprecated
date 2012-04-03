#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BSPLINE_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BSPLINE_TCC 1

#include <cassert>
#include <iostream>
#include <lawa/math/math.h>

namespace lawa {

template <typename T>
BSpline<T,Orthogonal,Interval,Multi>::BSpline(const MRA<T,Orthogonal,Interval,Multi> &_mra)
    : mra(_mra), d(_mra.d), initialticsize(pow2i<T>(-3))
{
    switch (d) {
        case 2:
            initialticsize = pow2i<T>(-2);
            break;

        case 3:
            initialticsize = pow2i<T>(-3);
            break;

        case 4:
            initialticsize = pow2i<T>(-2);
            break;

        default: std::cerr << "BSpline<T,Orthogonal,Interval,Multi> not yet realized"
                    " for d = " << d << ". Stopping." << std::endl;
                    exit(-1);
    }

}
    
template <typename T>
BSpline<T,Orthogonal,Interval,Multi>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2ih<T>(2*j*deriv+j) * mra._leftEvaluator[k](pow2i<T>(j)*x, deriv);
    }
    
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
        return pow2ih<T>(2*j*deriv+j) * 
               mra._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
    
    // right part
    int type  = (int)(k+1 - (mra.cardI(j) - mra._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1. - mra._numLeftParts)/mra._numInnerParts);
    return pow2ih<T>(2*j*deriv+j) * mra._rightEvaluator[type](pow2i<T>(j)*x-shift, deriv);
}
    
template <typename T>
Support<T>
BSpline<T,Orthogonal,Interval,Multi>::support(int j, long k) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2i<T>(-j) * mra._leftSupport[k];
    }
    
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
        return pow2i<T>(-j) * (mra._innerSupport[type]+shift);
    }
    
    // right part
    int type  = (int)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
    return pow2i<T>(-j) * (mra._rightSupport[type]+shift);
}

template <typename T>
DenseVector<Array<T> >
BSpline<T,Orthogonal,Interval,Multi>::singularSupport(int j, long k) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        return pow2i<T>(-j) * mra._leftSingularSupport[k];
    }
    
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        long shift = iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
        DenseVector<Array<T> > result = mra._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
    
    // right part
    int type  = (int)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    long shift = pow2i<long>(j)-1;
    //long shift = iceil<T>((k+1. - mra._numLeftParts)/mra._numInnerParts);
    DenseVector<Array<T> > result = mra._rightSingularSupport[type];
    result += shift;
    return pow2i<T>(-j) * result;
}

template <typename T>
T
BSpline<T,Orthogonal,Interval,Multi>::tic(int j) const
{
    //return pow2i<T>(-(j+3));
    return initialticsize*pow2i<T>(-j);
}

template <typename T>
DenseVector<Array<long double> > *
BSpline<T,Orthogonal,Interval,Multi>::getRefinement(int j, long k, long &shift, long &offset) const
{
    // left boundary
    if (k<mra._numLeftParts) {
        offset = mra.leftOffsets[k];
        shift  = 0L;
        return &(mra.leftRefCoeffs[k]);
    }
    // inner part
    if (k<mra.cardIL()+mra.cardII(j)) {
        int type  = (int)((k-mra._numLeftParts) % mra._numInnerParts);
        shift = (long)iceil<T>((k+1.-mra._numLeftParts)/mra._numInnerParts);
        offset = mra.innerOffsets[type];
        return &(mra.innerRefCoeffs[type]);
    }
    // right part
    int type  = (int)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    shift = pow2i<long>(j)-1;
    offset = mra.rightOffsets[type];
    return &(mra.rightRefCoeffs[type]);
}


} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_BSPLINE_TCC
