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

#ifndef LAWA_ADAPTIVE_HELMHOLTZ_H
#define LAWA_ADAPTIVE_HELMHOLTZ_H 1

#include <lawa/flensforlawa.h>
#include <lawa/adaptive/problem.h>
#include <lawa/adaptive/waveletindex.h>

#include <string>

namespace lawa {

template <typename T,Construction Cons>
struct Helmholtz
    : Problem<T,Cons>
{
    Helmholtz(const Basis<T,Primal,Interval,Cons> &_basis,
              const Basis<T,Dual,Interval,Cons> &_basis_);

    T
    preconditioner(const WaveletIndex<T,Cons> &l1) const;

    T
    bilinearForm(const WaveletIndex<T,Cons> &l1,
                 const WaveletIndex<T,Cons> &l2) const;
    
    T
    rhs(const WaveletIndex<T,Cons> &l1) const;

    bool
    compression(const WaveletIndex<T,Cons> &l1,
                const WaveletIndex<T,Cons> &l2, int J) const;

    bool compression(int j1, int k1, XType xtype1,
                     int j2, int k2, XType xtype2, int J) const;

    T
    alpha(int j) const;

    void
    rescale(DenseVector<Array<T> > &x, T s) const;

    void
    rescale(Coefficient<Lexicographical,T,Cons> &coeff, T s) const;

    void
    rescale(Coefficient<AbsoluteValue,T,Cons> &coeff, T s) const;

    void
    exactSolution(DenseVector<Array<T> > &coefficient, int Level);

    Coefficient<Lexicographical,T,Cons>
    exactSolution(int J) const;

    Coefficient<Lexicographical,T,Cons>
    exactSolution(T eps) const;

    void
    output(const Coefficient<Lexicographical,T,Cons> &u_exact,
           const Coefficient<Lexicographical,T,Cons> &U_lambda,
           const std::string &mode, T delta=0);

    std::string
    name() const;

    using Problem<T,Cons>::P;
    using Problem<T,Cons>::A;
    using Problem<T,Cons>::basis;
    using Problem<T,Cons>::basis_;
    using Problem<T,Cons>::Lambda;
    using Problem<T,Cons>::LambdaCheck;
    using Problem<T,Cons>::f;
    using Problem<T,Cons>::CA;
    using Problem<T,Cons>::cA;
    using Problem<T,Cons>::kappa;
    using Problem<T,Cons>::s;
    using Problem<T,Cons>::t;
    using Problem<T,Cons>::Error_Estimate_factor;
};

} // namespace lawa

#include <lawa/adaptive/helmholtz.tcc>

#endif // LAWA_ADAPTIVE_HELMHOLTZ_H
