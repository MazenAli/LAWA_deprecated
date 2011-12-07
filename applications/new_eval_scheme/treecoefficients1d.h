/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_TREEFCOEFFICIENTS1D_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_TREEFCOEFFICIENTS1D_H 1

#include <ext/hash_map>
#include <tr1/unordered_map>

#include <iostream>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

template <typename T>
struct CoefficientsByLevel : public std::tr1::unordered_map<long, T>
{

    typedef typename std::tr1::unordered_map<long, T>::const_iterator const_it;
    typedef typename std::tr1::unordered_map<long, T>::value_type     val_type;

    CoefficientsByLevel<T>&
    operator=(const CoefficientsByLevel<T> &_coeff);

    CoefficientsByLevel<T>&
    operator+=(const CoefficientsByLevel<T> &_coeff);
};

template <typename T>
std::ostream& operator<<(std::ostream &s, const CoefficientsByLevel<T> &_treecoeff);

template <typename T>
struct TreeCoefficients1D
{
    typedef typename CoefficientsByLevel<T>::value_type                         val_type;
    typedef typename CoefficientsByLevel<T>::const_iterator                     const_by_level_it;
    typedef typename CoefficientsByLevel<T>::iterator                           by_level_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator    const_coeff1d_it;

    CoefficientsByLevel<T> bylevel[JMAX+1];

    TreeCoefficients1D(void);

    TreeCoefficients1D(unsigned int n);

    TreeCoefficients1D<T>&
    operator=(const Coefficients<Lexicographical,T,Index1D> &_coeff);

    TreeCoefficients1D<T>&
    operator-=(const Coefficients<Lexicographical,T,Index1D> &_coeff);

    const CoefficientsByLevel<T>&
    operator[](short j)  const;

    CoefficientsByLevel<T>&
    operator[](short j);

    T
    norm(T factor);

};

template <typename T>
std::ostream& operator<<(std::ostream &s, const TreeCoefficients1D<T> &_treecoeff);


}   // namespace lawa

#include <applications/new_eval_scheme/treecoefficients1d.tcc>

#endif  //  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_TREEFCOEFFICIENTS1D_H
