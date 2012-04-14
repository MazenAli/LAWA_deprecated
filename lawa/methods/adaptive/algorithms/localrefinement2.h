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

#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_LOCALREFINEMENT2_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_LOCALREFINEMENT2_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename PrimalBasis>
class LocalRefinement2
{
    typedef typename PrimalBasis::T T;

    typedef typename PrimalBasis::RefinementBasis                              RefinementBasis;

    typedef flens::DenseVector<flens::Array<long double> >                     DenseVectorLD;

    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    typedef typename CoefficientsByLevel<T>::const_it                          const_coeffbylevel_it;

    public:


        LocalRefinement2(const PrimalBasis &_basis);

        const PrimalBasis     &basis;
        const RefinementBasis &refinementbasis;

        // Computes the common local refinement of a scaling coefficient vector (including
        // multiscaling!!) and a multilevel wavelet coefficient vector (including multiwavelets!!).
        void
        reconstruct(const Coefficients<Lexicographical,T,Index1D> &u, int j_scaling,
                    Coefficients<Lexicographical,T,Index1D> &u_loc_single) const;

        // Computes the common local refinement of a b-spline vector (not multiscaling!!) and a
        // wavelet (including multi-wavelets!!).
        void
        reconstruct(const CoefficientsByLevel<T> &u_bspline, int j_bspline,
                    const CoefficientsByLevel<T> &u_wavelet, int j_wavelet,
                    CoefficientsByLevel<T> &u_loc_single,    int &refinement_j) const;

        //Computes the local refinement of a multiscaling representation
        void
        reconstructMultiScaling(const CoefficientsByLevel<T> &u_scaling, const int &j,
                                CoefficientsByLevel<T> &u_loc_single, int &refinement_j) const;

        //Computes the local refinement of a B-spline (not multiscaling!!)
        void
        reconstructBSpline(const int &j, const long &k, T coeff,
                           CoefficientsByLevel<T> &u_loc_single, int &refinement_j) const;

        //Computes the local refinement of a Wavelet (also multiwavelet!!)
        void
        reconstructWavelet(const int &j, const long &k, T coeff,
                           CoefficientsByLevel<T> &u_loc_single, int &refinement_j) const;
};


}   // namespace lawa

#include <lawa/methods/adaptive/algorithms/localrefinement2.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_LOCALREFINEMENT2_H
