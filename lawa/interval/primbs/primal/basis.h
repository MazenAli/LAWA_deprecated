/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_INTERVAL_PRIMBS_PRIMAL_BASIS_H
#define LAWA_INTERVAL_PRIMBS_PRIMAL_BASIS_H 1

#include <lawa/basis.h>
#include <lawa/interval/primbs/dual/mra.h>
#include <lawa/interval/primbs/primal/mra.h>

namespace lawa {
    
template <typename T>
class Basis<T,Primal,Interval,Primbs>
{
    public:
        Basis(int _d, int _d_, int j=-1);

        int
        level() const;

        void
        setLevel(int j) const;

        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();

        // cardinalities of whole, left, inner, right index sets (primal).
        int
        cardJ(int j) const;

        int
        cardJL(int j=-1) const;

        int
        cardJI(int j) const;

        int
        cardJR(int j=-1) const;

        // ranges of whole, left, inner, right index sets (primal).
        const Range<int>
        rangeJ(int j) const;

        const Range<int>
        rangeJL(int j=-1) const;

        const Range<int>
        rangeJI(int j) const;

        const Range<int>
        rangeJR(int j=-1) const;

        MRA<T,Primal,Interval,Primbs> mra;
        MRA<T,Dual,Interval,Primbs>  mra_;

        RefinementMatrix<T,Interval,Primbs> M1;

        const int d, d_, mu;   // mu = mu(d) = d&1.
        const int min_j0;      // minimal allowed(!) level;
        const int j0;          // minimal used(!) level.

    private:
        void
        _calcM1();
            
        DenseVector<Array<int> > _bc;    // the boundary conditions
                                           // bc(0) = 1 -> Dirichlet BC left.
                                           // bc(1) = 1 -> Dirichlet BC right.

        mutable int _j;                // the current level.

};

} // namespace lawa

#include <lawa/interval/primbs/primal/basis.tcc>

#endif // LAWA_INTERVAL_PRIMBS_PRIMAL_BASIS_H