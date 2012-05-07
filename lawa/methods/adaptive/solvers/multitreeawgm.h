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

#ifndef  LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_H
#define  LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_H 1

#include <map>
#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/algorithms/algorithms.h>

namespace lawa {

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
struct MultiTreeAWGM {

    typedef typename LocalOperator::T T;

    typedef typename IndexSet<Index>::const_iterator                          const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index>::iterator          coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator    const_coeff_it;

    MultiTreeAWGM(const Basis &_basis, LocalOperator &_A, RHS &_F, Preconditioner &_Prec);

    void
    setParameters(T _alpha, T _gamma);

    void
    cg_solve(Coefficients<Lexicographical,T,Index> &u, T _eps, const char *filename, int NumOfIterations=100, T EnergyNorm=0.);

    const Basis      &basis;
    LocalOperator    &A;
    RHS              &F;
    Preconditioner   &Prec;
    T                alpha, gamma;
    T                eps;
    size_t           hashMapSize;
};


}    //namespace lawa

#include <lawa/methods/adaptive/solvers/multitreeawgm.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_H

