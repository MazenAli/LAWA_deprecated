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

#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_LOCALREFINEMENT_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_LOCALREFINEMENT_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename PrimalBasis>
class LocalRefinement
{
    typedef typename PrimalBasis::T T;
    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    typedef typename CoefficientsByLevel<T>::const_it                          const_coeffbylevel_it;

    public:


        LocalRefinement(const PrimalBasis &_basis, bool withDirichletBC);

        /*
         * Expects a vector u_multi_j of wavelets _and_ scaling functions on level j.
         * Computes the _local_ representation in terms of scaling functions.
         * _Adds_ the coefficients to a vector _u_loc_single_jP1 .
         */
        void
        reconstruct(const Coefficients<Lexicographical,T,Index1D> &u_multi_j,
                    Coefficients<Lexicographical,T,Index1D> &u_loc_single_jP1) const;

        void
        reconstruct(const Coefficients<Lexicographical,T,Index1D> &u_multi_j, int j,
                    Coefficients<Lexicographical,T,Index1D> &u_loc_single_jP1) const;

        void
        reconstruct(const CoefficientsByLevel<T> &u_scaling, const CoefficientsByLevel<T> &u_wavelet,
                    int j, CoefficientsByLevel<T> &u_loc_single_jP1) const;

        /*
         * Expects an coefficient entry (wavelets or scaling function) and computes the _local_ refinement.
         * Result is _added_ to u_loc_single_jP1
         */
        void
        reconstruct(const short &j, const long &k, const XType &xtype, T coeff,
                    Coefficients<Lexicographical,T,Index1D> &u_loc_single) const;

        void
        reconstruct(const short &j, const long &k, const XType &xtype, T coeff,
                    CoefficientsByLevel<T> &u_loc_single) const;

        /*
         * Expects a vector <Phi_{j+1},v> (u_loc_single) with scaling functions on level j+1.
         * Computes the _local_ representation
         *  | |  M_{j,0}^T |                   |
         *  | |  M_{j,1}^T | < \Phi_{j+1},v>   |_{\Lambda}
         *  and _adds_ it to the result vector
         */
        void
        decompose_(Coefficients<Lexicographical,T,Index1D> &u_loc_single,
                   const IndexSet<Index1D> &Lambda,
                   Coefficients<Lexicographical,T,Index1D> &u_loc_single_jM1,
                   Coefficients<Lexicographical,T,Index1D> &u_multi) const;

        void
        decompose_(CoefficientsByLevel<T>  &u_loc_single, int j,
                   CoefficientsByLevel<T>  &u_scaling,
                   CoefficientsByLevel<T>  &u_wavelet) const;

        T
        decompose_(Coefficients<Lexicographical,T,Index1D> &u_loc_single,
                   const short &j, const long &k, const XType &xtype) const;

        T
        decompose_(const CoefficientsByLevel<T> &u_loc_single,
                   const short &j, const long &k, const XType &xtype) const;



        /*
         * "Graphical test" if refinement is correct.
         */
        void
        test_reconstruct(int j, long k, XType xtype) const;

        const PrimalBasis &basis;

    private:
        int numLeftBoundaryScalings, numRightBoundaryScalings;
        int numLeftBoundaryWavelets, numRightBoundaryWavelets;
        int inner_scaling_offset, inner_wavelet_offset;
        DenseVector<Array<DenseVector<Array<T> > > > leftM0, rightM0;
        DenseVector<Array<T> > innerM0;
        DenseVector<Array<DenseVector<Array<T> > > > leftM1, rightM1;
        DenseVector<Array<T> > innerM1;
};


}   // namespace lawa

#include <lawa/methods/adaptive/algorithms/localrefinement.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_LOCALREFINEMENT_H
