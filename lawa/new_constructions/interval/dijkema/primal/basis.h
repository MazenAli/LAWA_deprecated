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

#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H 1

#include <lawa/aux/integer.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/interval/dijkema/dual/mra.h>
#include <lawa/constructions/interval/primal/bspline.h>
#include <lawa/constructions/interval/dijkema/primal/mra.h>
#include <lawa/constructions/interval/primal/wavelet.h>

namespace lawa {
    
template <typename _T>
class Basis<_T,Primal,Interval,Dijkema>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = Interval;
        static const Construction Cons = Dijkema;

        typedef BasisFunction<T,Primal,Interval,Dijkema> BasisFunctionType;
        typedef BSpline<T,Primal,Interval,Dijkema> BSplineType;
        typedef Wavelet<T,Primal,Interval,Dijkema> WaveletType;

        Basis(int _d, int _d_, int j=-1);
        
        int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();

        const BasisFunctionType &
        generator(XType xtype) const;

        // cardinalities of whole, left, inner, right index sets (primal).
        Integer
        cardJ(int j) const;

        Integer
        cardJL(int j=-1) const;

        Integer
        cardJI(int j) const;

        Integer
        cardJR(int j=-1) const;

        // ranges of whole, left, inner, right index sets (primal).
        const Range<Integer>
        rangeJ(int j) const;

        const Range<Integer>
        rangeJL(int j=-1) const;

        const Range<Integer>
        rangeJI(int j) const;

        const Range<Integer>
        rangeJR(int j=-1) const;

        MRA<T,Primal,Interval,Dijkema> mra;
        MRA<T,Dual,Interval,Dijkema>  mra_;

        RefinementMatrix<T,Interval,Dijkema> M1;

        const int d, d_, mu;   // mu = mu(d) = d&1.
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level.

    private:
        DenseVector<Array<int> > _bc;    // the boundary conditions
                                           // bc(0) = 1 -> Dirichlet BC left.
                                           // bc(1) = 1 -> Dirichlet BC right.

        mutable int _j;                // the current level.

    public:
        Wavelet<T,Primal,Interval,Dijkema> psi;
};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_BASIS_H
