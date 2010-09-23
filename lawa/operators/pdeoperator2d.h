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


#ifndef LAWA_OPERATORS_PDEOPERATOR2D_H
#define LAWA_OPERATORS_PDEOPERATOR2D_H 1

#include <lawa/enum.h>
#include <lawa/integrals.h>
#include <lawa/adaptive/index.h>
#include <lawa/operators/rieszoperator1d.h>
#include <lawa/operators/weaklaplaceoperator1d.h>


namespace lawa {

template <typename T, typename Basis2D>
class PDEOperator2D {
	typedef typename Basis2D::FirstBasisType  Basis_x;
	typedef typename Basis2D::SecondBasisType Basis_y;


public:
	const Basis2D &basis;
	WeakLaplaceOperator1D<T,Basis_x> diffusion_x;
	RieszOperator1D<T,Basis_x>       reaction_x;
	WeakLaplaceOperator1D<T,Basis_y> diffusion_y;
	RieszOperator1D<T,Basis_y>       reaction_y;
	T c;

    PDEOperator2D(const Basis2D &basis, T c);

    T getc() const;

    const Basis2D& getBasis() const;

    T
    operator()(XType row_xtype_x, int j1_x, int k1_x,
               XType row_xtype_y, int j1_y, int k1_y,
               XType col_xtype_x, int j2_x, int k2_x,
               XType col_xtpye_y, int j2_y, int k2_y) const;

    T
    operator()(const Index2D &row_index, const Index2D &col_index) const;
};

}    // namespace lawa

#include <lawa/operators/pdeoperator2d.tcc>

#endif // LAWA_ADAPTIVE_PDEOPERATOR2D_H

