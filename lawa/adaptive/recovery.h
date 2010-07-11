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

#ifndef LAWA_ADAPTIVE_RECOVERY_H
#define LAWA_ADAPTIVE_RECOVERY_H 1

#include <lawa/enum.h>
#include <lawa/adaptive/waveletcoefficient.h>

namespace lawa {

template <typename T, RecoveryType R>
struct Nonlinearity
{
};

template <typename T>
struct Nonlinearity<T,BU>
{
    Nonlinearity(T (*_C_F)(T), T _t, T _r, T _s, T (*_F)(T x), T (*_dF)(T x));
    
    T t, r, s;
    T (*C_F)(T);
    T (*F)(T x);
    T (*dF)(T x);
};

template <typename T>
struct Nonlinearity<T,CDD>
{
    Nonlinearity(T _t, T _gamma, T (*_F)(T x), T (*_dF)(T x));
    
    T t, gamma;
    T (*F)(T x);
    T (*dF)(T x);
};

template <typename T,RecoveryType R,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Recovery(const Problem<T,Cons> &P,
         const Nonlinearity<T,R> &F,
         const Coefficient<Lexicographical,T,Cons> &u, T eps, int J);

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Recovery_APPLY(const StiffnessMatrix<T,Cons> &A,
               const Nonlinearity<T,BU> &F,
               const Coefficient<Lexicographical,T,Cons> &u, T eps,
               int J=std::numeric_limits<int>::max());

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Recovery_APPLY(const StiffnessMatrix<T,Cons> &A,
               const Nonlinearity<T,CDD> &F,
               const Coefficient<Lexicographical,T,Cons> &u, T eps,
               int J=std::numeric_limits<int>::max());

//----------------------------------------------------------------------------

template <typename T,RecoveryType R,Construction Cons>
IndexSet<T,Cons>
prediction(const Nonlinearity<T,R> &F,
           const Coefficient<Lexicographical,T,Cons> &u, T delta);

template <typename T,Construction Cons>
IndexSet<T,Cons>
prediction(const IndexSet<T,Cons> Lambda,
           const Nonlinearity<T,BU> &F,
           const Coefficient<Lexicographical,T,Cons> &u, T eps);

template <typename T,Construction Cons>
IndexSet<T,Cons>
prediction(const Nonlinearity<T,CDD> &F,
           const Coefficient<Lexicographical,T,Cons> &u, T eps,
           int J=std::numeric_limits<int>::max());

template <typename T,Construction Cons>
IndexSet<T,Cons>
local_prediction(const IndexSet<T,Cons> &Lambda);

template <typename T,Construction Cons>
IndexSet<T,Cons>
local_prediction(const Coefficient<Lexicographical,T,Cons> &u);

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
localReconstruction(const Coefficient<Lexicographical,T,Cons> &v,
                    const IndexSet<T,Cons> &Gamma);

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
localReconstruction_(const Coefficient<Lexicographical,T,Cons> &v,
                     const IndexSet<T,Cons> &Gamma);

template <typename T,RecoveryType R,Construction Cons>
Coefficient<Lexicographical,T,Cons>
quasiInterpolation(const Nonlinearity<T,R> &F,
                   const Coefficient<Lexicographical,T,Cons> &v);

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
localDecomposition(const Coefficient<Lexicographical,T,Cons> &v);

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
localDecomposition_(const Coefficient<Lexicographical,T,Cons> &v);

} // namespace lawa

#include <lawa/adaptive/recovery.tcc>

#endif // LAWA_ADAPTIVE_RECOVERY_H 
