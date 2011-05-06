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


#ifndef LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDL2SCALARPRODCUT1D_H
#define LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDL2SCALARPRODCUT1D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/functiontypes/exponentialweightfunction1d.h>
#include <lawa/integrals/integral.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename Basis>
class WeightedL2ScalarProduct1D {

    public:

        const Basis                     &basis;
        ExponentialWeightFunction1D<T>  exponentialweightfunction;

        WeightedL2ScalarProduct1D(const Basis& _basis, const T _eta=0.,
                                  T R1=0., T R2=1., int order=10);

        T
        operator()(XType xtype1, int j1, int k1,
                   XType xtype2, int j2, int k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;

    private:

        const T                        eta;
        Function<T>                    weight;
        IntegralF<Gauss, Basis, Basis> integral_w;
        Integral<Gauss, Basis, Basis>  integral;


};

}   //namespace lawa

#include <lawa/operators/pdeoperators1d/weightedL2scalarproduct1d.tcc>

#endif  // LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDL2SCALARPRODCUT1D_H
