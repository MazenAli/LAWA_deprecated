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

#ifndef TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_REFERENCESOLUTIONS1D_H
#define TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_REFERENCESOLUTIONS1D_H 1

#include <iostream>
#include <lawa/settings/enum.h>
#include <lawa/operators/pdeoperators1d/helmholtzoperator1d.h>
#include <lawa/operators/pdeoperators2d/helmholtzoperator2d.h>

namespace lawa {


template <typename T, typename Basis, typename BilinearForm>
struct ReferenceSolution1D
{
};

template<typename T, typename Basis>
struct ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >
{
    static int nr;
    static T c;
    static DomainType domain;

    static T
    exact(T x, int deriv);

    static DenseVector<Array<T> >
    sing_pts;

    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
    deltas;

    static void
    setExample(int _nr, const HelmholtzOperator1D<T,Basis> &a, DomainType domain);

    static T
    exact(T x);

    static T
    d_exact(T x);

    static T
    rhs(T x);

    static T
    smoothened_rhs(T x);

    static T
    H1norm();
};


} // namespace lawa

#include <tutorials/examples/referencesolutions/referencesolutions1d.tcc>


#endif // TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_H
