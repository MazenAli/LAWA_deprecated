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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_UNIDIRECTIONALOPERATOR_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_UNIDIRECTIONALOPERATOR_H 1

#include <lawa/flensforlawa.h>
#include <lawa/settings/enum.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/operators/localoperators/localoperator1d.h>

namespace lawa {

template <typename Index, CoordinateDirection CoordX, typename LocalOperator1D,
                          CoordinateDirection NotCoordX, typename NotCoordXIndex>
struct UniDirectionalOperator
{
    typedef typename LocalOperator1D::T T;

    typedef typename LocalOperator1D::TrialWaveletBasis                    TrialBasis_x1;
    typedef typename LocalOperator1D::TestWaveletBasis                     TestBasis_x1;

    typedef typename TreeCoefficients1D<T>::const_by_level_it              const_by_level_it;
    typedef typename TreeCoefficients1D<T>::by_level_it                    by_level_it;

    typedef AlignedCoefficients<T,Index,NotCoordXIndex,Index1D,NotCoordX>  NotCoordXAlignedCoefficients;

    UniDirectionalOperator(LocalOperator1D &_localOperator1D);

    LocalOperator1D &localOperator1D;

};

}   // namespace lawa

#endif  // LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_UNIDIRECTIONALOPERATOR_H
