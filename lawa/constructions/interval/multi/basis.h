/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008-2011 Sebastian Kestler, Kristina Steih,
                         Alexander Stippler, Mario Rometsch.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Orthogonal,Interval,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Multi;
        
        typedef Basis<T,Orthogonal,Interval,MultiRefinement> RefinementBasis;
        typedef BasisFunction<T,Orthogonal,Interval,Multi>   BasisFunctionType;
        typedef BSpline<T,Orthogonal,Interval,Multi>         BSplineType;
        typedef Wavelet<T,Orthogonal,Interval,Multi>         WaveletType;

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

        //--- cardinalities of whole, left, inner, right index set.
        int
        cardJ(const int j) const;
        
        int
        cardJL(const int j=-1) const;

        int
        cardJI(const int j) const;

        int
        cardJR(const int j=-1) const;
    
        //--- ranges of whole, left, inner, right index set.
        const flens::Range<int>
        rangeJ(const int j) const;
    
        const flens::Range<int>
        rangeJL(const int j=-1) const;

        const flens::Range<int>
        rangeJI(const int j) const;

        const flens::Range<int>
        rangeJR(const int j=-1) const;
    
        MRA<T,Orthogonal,Interval,Multi> mra;

        const int d;
        const int j0;          // minimal used(!) level.

        Wavelet<T,Orthogonal,Interval,Multi> psi;

        Basis<T,Orthogonal,Interval,MultiRefinement> refinementbasis;

    
    private:
        DenseVector<Array<int> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
                                       // bc(1) = 1 -> Dirichlet BC right.
        
        mutable int _j;                // the current level.
    
        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        friend class Wavelet<T,Orthogonal,Interval,Multi>;

        unsigned int _numLeftParts, 
                     _numInnerParts, 
                     _numRightParts;
        Evaluator *_leftEvaluator, 
                  *_innerEvaluator, 
                  *_rightEvaluator;
        Support<T> *_leftSupport, 
                   *_innerSupport, 
                   *_rightSupport;
        DenseVector<Array<T> > *_leftSingularSupport, 
                               *_innerSingularSupport, 
                               *_rightSingularSupport;
        
        // Refinement coefficients for representation in B-splines, only d \in \{2,3\} so far,
        // because d=4 for multiwavelets requires cubic Hermite B-splines which are not used
        // in the biorthogonal construction
        DenseVector<Array<long double> > *_leftRefCoeffs,
                                         *_innerRefCoeffs,
                                         *_rightRefCoeffs;

        long *_leftOffsets,
             *_innerOffsets,
             *_rightOffsets;

        int  _addRefinementLevel;

};

} // namespace lawa

#include <lawa/constructions/interval/multi/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_H
