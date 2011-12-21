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

#ifndef APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR2D_H
#define APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR2D_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <applications/new_eval_scheme/source/localoperator.h>

namespace lawa {

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
struct LocalOperator2D {

    typedef typename Basis::T T;

    typedef IndexSet<Index2D>::const_iterator                                  const_set2d_it;

    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator   const_coeff2d_it;

    typedef typename TreeCoefficients1D<T>::const_by_level_it                  const_by_level_it;
    typedef typename TreeCoefficients1D<T>::by_level_it                        by_level_it;

    typedef AlignedCoefficients<T,Index2D,Index1D,Index1D>                     alignedCoefficients;

    LocalOperator2D(const Basis &_basis, const LocalOperator1 &_localoperator1,
                    const LocalOperator2 &_localoperator2);

    void
    evalAA(const Coefficients<Lexicographical,T,Index2D> &v,
           Coefficients<Lexicographical,T,Index2D> &AAv);

    void
    evalUI(const Coefficients<Lexicographical,T,Index2D> &v,
           Coefficients<Lexicographical,T,Index2D> &UIv);

    void
    evalIA(const Coefficients<Lexicographical,T,Index2D> &UIv,
           Coefficients<Lexicographical,T,Index2D> &AAv);

    const Basis                   &basis;
    const LocalOperator1          &localoperator1;
    const LocalOperator2          &localoperator2;
};

}   // namespace lawa

#include <applications/new_eval_scheme/source/localoperator2d.tcc>

#endif  // APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR2D_H
