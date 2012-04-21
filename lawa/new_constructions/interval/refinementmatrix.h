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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_REFINEMENTMATRIX_H
#define LAWA_CONSTRUCTIONS_INTERVAL_REFINEMENTMATRIX_H 1

#include <lawa/aux/integer.h>
#include <lawa/flensforlawa.h>
#include <lawa/constructions/refinementmatrix.h>

namespace flens {

using namespace lawa;
using namespace cxxblas;

template <typename T, Construction Cons>
class RefinementMatrix<T, Interval, Cons>
    : public Matrix<RefinementMatrix<T, Interval, Cons> >
{
    public:
        typedef T ElementType;
                
        RefinementMatrix();
        
        RefinementMatrix(int nLeft, int nRight,
                         const GeMatrix<FullStorage<T, ColMajor> > &A,
                         int _min_j0, int cons_j);

        const typename DenseVector<Array<T> >::ConstView
        operator()(int j, const Underscore<Integer> &u, Integer col) const;
        
        Range<Integer>
        rows() const;

        Range<Integer>
        cols() const;

        Integer
        numRows() const;

        Integer
        numCols() const;

        Integer
        firstRow() const;

        Integer
        lastRow() const;

        Integer
        firstCol() const;

        Integer
        lastCol() const;

        int
        level() const;

        void
        setLevel(int j) const;

        DenseVector<Array<DenseVector<Array<T> > > > left, right;
        DenseVector<Array<T> > leftband, rightband;
        DenseVector<Array<Integer> > lengths;
        int min_j0;

    private:
        void
        _extractMasks(const GeMatrix<FullStorage<T,ColMajor> > &A);        
        
        int _cons_j;
        mutable int _j;
        Integer _firstRow, _firstCol, _lastRow, _lastCol;
        mutable Integer _additionalRows, _additionalCols;
};

template <typename T, Construction Cons>
struct TypeInfo<RefinementMatrix<T,Interval,Cons> >
{
    typedef RefinementMatrix<T,Interval,Cons> Impl;
    typedef T                                 ElementType;
    
};

template <typename X, Construction Cons, typename Y>
void
mv(Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,Interval,Cons> &A,
   const DenseVector<X> &x, typename X::ElementType beta, DenseVector<Y> &y);

} // namespace flens

#endif // LAWA_CONSTRUCTIONS_REFINEMENTMATRIX_H
