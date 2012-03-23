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

#ifndef LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BASIS_H
#define LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BASIS_H 1

#include <lawa/aux/integer.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/periodic/refinementmatrix.h>
#include <lawa/constructions/periodic/primal/mra.h>
#include <lawa/constructions/periodic/dual/mra.h>
#include <lawa/constructions/periodic/dual/wavelet.h>

namespace lawa {

template <typename _T>
class Basis<_T,Dual,Periodic,CDF>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Dual;
        static const DomainType Domain = Periodic;
        static const Construction Cons = CDF;

        typedef BasisFunction<T,Dual,Periodic,CDF> BasisFunctionType;
        typedef BSpline<T,Dual,Periodic,CDF> BSplineType;
        typedef Wavelet<T,Dual,Periodic,CDF> WaveletType;

        Basis(unsigned int _d, unsigned int _d_, int j=0);

        int
        level() const;

        void
        setLevel(int j) const;

        const BasisFunctionType &
        generator(XType xtype) const;

        int
        cardJ_(int j) const;

        Range<Integer>
        rangeJ_(int j) const;

        const unsigned int d, d_; 
        const int j0;
        MRA<T,Primal,Periodic,CDF> mra;
        MRA<T,Dual,Periodic,CDF> mra_;
        Wavelet<T,Dual,Periodic,CDF> psi_;
        RefinementMatrix<T,Periodic,CDF> M1_;
        
    private:
        mutable int _j;
};

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BASIS_H
