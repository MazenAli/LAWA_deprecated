/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef LAWA_PRECONDITIONERS_SPACETIMEPRECONDITIONERS_RIGHTNORMPRECONDITIONER2D_C_PG_H
#define LAWA_PRECONDITIONERS_SPACETIMEPRECONDITIONERS_RIGHTNORMPRECONDITIONER2D_C_PG_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integrals.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
class RightNormPreconditioner2D_c_PG
{
    typedef typename TrialBasis::FirstBasisType FirstBasis_Trial;
    typedef typename TestBasis::FirstBasisType FirstBasis_Test;
    typedef typename TrialBasis::SecondBasisType SecondBasis_Trial;
    typedef typename TestBasis::SecondBasisType SecondBasis_Test;

    public:
        RightNormPreconditioner2D_c_PG(const TrialBasis& trialbasis, const TestBasis& testbasis, T s=2.);  //s=2: A: H^1 -> H^{-1}

        T
        operator()(XType xtype1, int j1, int k1,
                   XType xtype2, int j2, int k2) const;

        T
        operator()(const Index2D &index) const;

    private:
        //const Basis2D &_basis;
        T              _s;  // scaling for certain classes of integral operators
        Integral<Gauss,FirstBasis_Test,FirstBasis_Trial>   _integral_t;
        Integral<Gauss,SecondBasis_Test,SecondBasis_Trial> _integral_x;
};

}   // namespace lawa

#include <lawa/preconditioners/spacetimepreconditioners/rightnormpreconditioner2d_c_pg.tcc>

#endif // LAWA_PRECONDITIONERS_SPACETIMEPRECONDITIONERS_RIGHTNORMPRECONDITIONER2D_C_PG_H
