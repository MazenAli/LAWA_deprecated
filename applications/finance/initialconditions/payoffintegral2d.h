/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Sebastian Kestler, Kristina Steih, Mario Rometsch, Alexander Stippler.

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

#ifndef APPLICATIONS_FINANCE_INITIALCONDITIONS_EUROPEANBASKET2D_H
#define APPLICATIONS_FINANCE_INITIALCONDITIONS_EUROPEANBASKET2D_H 1

#include <lawa/constructions/constructions.h>
#include <lawa/integrals/integrals.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <applications/finance/initialconditions/adaptivepayoffquadrature2d.h>
#include <applications/finance/options/options.h>
#include <applications/finance/processes/processes.h>

namespace lawa {

template <OptionTypenD OType, ProcessType2D PType, typename Basis>
struct PayoffIntegral2D
{
    typedef typename Basis::T T;
    typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

    PayoffIntegral2D(const Option2D<T,OType> &_option,
                     ProcessParameters2D<T,PType> &_processparameters,
                     const Basis &_basis, const T _left1=0., const T _right1=1.,
                     const T _left2=0., const T _right2=1.);

    T
    operator()(Index2D &lambda) const;

    T
    payoff(T x1, T x2) const;

    T
    integrand(T x1, T x2) const;

    T
    operator()(Index2D &index2d);

    const Option2D<T,OType>         &option;
    ProcessParameters2D<T,PType>    &processparameters;
    const Basis                     &basis;
    AdaptivePayoffQuadrature2D<OType, PayoffIntegral2D<OType,PType,Basis> > adapquad;
    const T                         left1, right1;
    const T                         left2, right2;
    const T                         RightmLeft1, SqrtRightmLeft1;
    const T                         RightmLeft2, SqrtRightmLeft2;

    mutable int j1, deriv1,
                j2, deriv2;
    mutable long k1, k2;
    mutable XType e1, e2;
};

}   // namespace lawa

#include <applications/finance/initialconditions/payoffintegral2d.tcc>

#endif  // APPLICATIONS_FINANCE_INITIALCONDITIONS_BASKETPUT2D_H
