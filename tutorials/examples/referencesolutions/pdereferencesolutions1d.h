/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_PDEREFERENCESOLUTIONS1D_H
#define TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_PDEREFERENCESOLUTIONS1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

/*
 * Reference solutions u and corresponding righthand sides for second order PDEs
 * with constant coefficients:
 *       - diffusion * u'' + convection * u' + reaction * u = f
 */

template<typename T>
struct PDEReferenceSolutions1D
{
    static int nr;

    static T diffusion, convection, reaction;

    static DomainType domain;

    static DenseVector<Array<T> > sing_pts;

    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > deltas;

    static void
    setExample(int _nr, T _diffusion, T _convection, T _reaction, DomainType domain);

    static T
    exact(T x, int deriv);

    static T
    exact(T x);

    static T
    d_exact(T x);

    static T
    rhs(T x);

    static T
    H1norm();
};


} // namespace lawa

#include <tutorials/examples/referencesolutions/pdereferencesolutions1d.tcc>


#endif // TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_PDEREFERENCESOLUTIONS1D_H
