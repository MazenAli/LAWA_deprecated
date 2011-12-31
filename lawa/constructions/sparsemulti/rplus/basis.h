#ifndef LAWA_CONSTRUCTIONS_SPARSEMULTI_RPLUS_BASIS_H
#define LAWA_CONSTRUCTIONS_SPARSEMULTI_RPLUS_BASIS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Primal,RPlus,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = RPlus;
        static const Construction Cons = SparseMulti;
        
        typedef BasisFunction<T,Primal,RPlus,SparseMulti> BasisFunctionType;
        typedef BSpline<T,Primal,RPlus,SparseMulti>       BSplineType;
        typedef Wavelet<T,Primal,RPlus,SparseMulti>       WaveletType;

        Basis(const int d, const int j=-1);
    
        virtual
        ~Basis();
    
        int
        level() const;
    
        void
        setLevel(const int j) const;
    
        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();
    
        const BasisFunctionType &
        generator(XType xtype) const;

        //--- cardinalities of left index set.
        long
        cardJL(const int j=-1) const;

        //--- ranges of left index set.
        const flens::Range<long>
        rangeJL(const int j=-1) const;
    
        MRA<T,Primal,RPlus,SparseMulti> mra;
    
        const int d;
        const int j0;          // minimal used(!) level.
    
    private:
        DenseVector<Array<int> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
        
        mutable int _j;                // the current level.
    
        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        friend class Wavelet<T,Primal,RPlus,SparseMulti>;

        unsigned int _numLeftParts,
                     _numInnerParts;
        Evaluator *_leftEvaluator, 
                  *_innerEvaluator;
        Support<T> *_leftSupport, 
                   *_innerSupport;
        DenseVector<Array<T> > *_leftSingularSupport, 
                               *_innerSingularSupport;
        DenseVector<Array<T> > _leftScalingFactors, _innerScalingFactors;
        
    public:
        Wavelet<T,Primal,RPlus,SparseMulti> psi;
};

} // namespace lawa

#include <lawa/constructions/sparsemulti/rplus/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_SPARSEMULTI_RPLUS_BASIS_H
