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

#ifndef LAWA_ADAPTIVE_PROBLEM_H
#define LAWA_ADAPTIVE_PROBLEM_H 1

#include <lawa/adaptive/coefficients.h>
#include <lawa/operators/operators.h>

namespace lawa {

template <typename T, typename Basis, typename BilinearForm>
class Problem1d
{
};

/*  Not to be used!! Structure not yet clear...
 *
template <typename T, typename Basis>
class Problem1d<T,Basis,HelmholtzOperator1d<T,Basis> >
{
private:
    const Basis &basis;
    const int d,d_;
    Coefficients<Lexicographical,T,Index1d> u, f;
    ReferenceSolution1d<T,Basis,HelmholtzOperator1d<T,Basis> > refsol;

public:
    Problem1d(const Basis &basis, int example_nr, const HelmholtzOperator1d<T,Basis> &a, DomainType domain);

    Coefficients<Lexicographical,T,Index1d>
    RHS(const IndexSet<Index1d> Lambda);

    Coefficients<Lexicographical,T,Index1d>
    RHS(T tol);
};

*/

}  //namespace lawa

//#include <adaptive/problem.tcc>

#endif // LAWA_ADAPTIVE_PROBLEM_H
