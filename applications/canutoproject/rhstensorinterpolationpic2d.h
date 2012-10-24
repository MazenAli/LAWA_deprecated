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

#ifndef APPLICATIONS_CANUTOPROJECT_RHSTENSORINTERPOLATIONPIC2D_H
#define APPLICATIONS_CANUTOPROJECT_RHSTENSORINTERPOLATIONPIC2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/algorithms/multitreeoperations.h>

namespace lawa {

template <typename Basis, typename LocalOperator>
struct RHSTensorInterpolationPic2D
{
    typedef typename LocalOperator::T T;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff2d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::iterator       coeff2d_it;

    RHSTensorInterpolationPic2D(const Basis &_basis, LocalOperator &_localOp,
                                const Coefficients<Lexicographical,T,Index2D> &_u);

    void
    initializeRHS(const Coefficients<Lexicographical,T,Index2D> &_f);

    T
    operator()(const Index2D &index);

    T
    operator()(XType xtype_x, int j_x, long k_x,
               XType xtype_y, int j_y, long k_y);


    const Basis                                     &basis;
    LocalOperator                                   &localOp;
    const Coefficients<Lexicographical,T,Index2D>   &u;
    Coefficients<Lexicographical,T,Index2D>         f;
};

}   // namespace lawa

#include <applications/canutoproject/rhstensorinterpolationpic2d.tcc>

#endif  // APPLICATIONS_CANUTOPROJECT_RHSTENSORINTERPOLATIONPIC2D_H
