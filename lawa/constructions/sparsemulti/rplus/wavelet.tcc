#include <cassert>
#include <iostream>

namespace lawa {

template <typename T>
Wavelet<T,Primal,RPlus,SparseMulti>::Wavelet(const Basis<T,Primal,RPlus,SparseMulti> &_basis)
    : basis(_basis), d(_basis.d), vanishingMoments(_basis.d)
{
    switch (d) {
            case 4: _numSplines = 4;
                    _max_support = Support<T>(-2.,2.);
                    break;

            default: std::cerr << "Wavelet<T,Primal,R,SparseMulti> not yet realized"
                                          " for d = " << d << ". Stopping." << std::endl;
                     exit(1);
    }
}
    
template <typename T>
Wavelet<T,Primal,RPlus,SparseMulti>::~Wavelet()
{
}

template <typename T>
T
Wavelet<T,Primal,RPlus,SparseMulti>::operator()(T x, int j, long k, unsigned short deriv) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2ih<T>(2*j*deriv+j) * basis._leftScalingFactors(k) *
                   basis._leftEvaluator[k](pow2i<T>(j)*x, deriv);
        }
        k-=(basis._numLeftParts-1);
        int type  = (int)(k % basis._numInnerParts);
        long shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
        //std::cerr << "type = " << type << ", shift = " << shift << " "  << std::endl;
        return pow2ih<T>(2*j*deriv+j) * basis._innerScalingFactors(type) *
               basis._innerEvaluator[type](pow2i<T>(j)*x-shift,deriv);
    }
}
    
template <typename T>
Support<T>
Wavelet<T,Primal,RPlus,SparseMulti>::support(int j, long k) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2i<T>(-j) * basis._leftSupport[k];
        }

        k-=(basis._numLeftParts-1);
        // inner part
        int type  = (int)(k % basis._numInnerParts);
        long shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
        return pow2i<T>(-j) * (basis._innerSupport[type]+shift);
    }
}

template <typename T>
Support<T>
Wavelet<T,Primal,RPlus,SparseMulti>::max_support() const
{
    return _max_support;
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Primal,RPlus,SparseMulti>::singularSupport(int j, long k) const
{
    if (d==4) {
        k-=1;
        // left boundary
        if (k<basis._numLeftParts) {
            return pow2i<T>(-j) * basis._leftSingularSupport[k];
        }
        k-=(basis._numLeftParts-1);
        // inner part
        int type  = (int)(k % basis._numInnerParts);
        long shift = 2*std::ceil((T(k) / T(basis._numInnerParts)));
        DenseVector<Array<T> > result = basis._innerSingularSupport[type];
        result += shift;
        return pow2i<T>(-j) * result;
    }
}
    
template <typename T>
T
Wavelet<T,Primal,RPlus,SparseMulti>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}
    
} // namespace lawa
