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

#ifndef APPLICATIONS_NEWEVALSCHEME_COMPOUNDLOCALOPERATOR_H
#define APPLICATIONS_NEWEVALSCHEME_COMPOUNDLOCALOPERATOR_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/algorithms/localrefinement2.h>

namespace lawa {

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator=SecondLocalOperator,
          typename FourthLocalOperator=SecondLocalOperator>
class CompoundLocalOperator {

    public:
        typedef typename FirstLocalOperator::T T;

        CompoundLocalOperator(const FirstLocalOperator  &_firstLocalOp,
                              const SecondLocalOperator &_secondLocalOp);

        CompoundLocalOperator(const FirstLocalOperator  &_firstLocalOp,
                              const SecondLocalOperator &_secondLocalOp,
                              const ThirdLocalOperator  &_thirdLocalOp);

        void
        eval(const Coefficients<Lexicographical,T,Index> &v,
             Coefficients<Lexicographical,T,Index> &auxiliary,
             Coefficients<Lexicographical,T,Index> &Av);

        int                               numOfLocalOp;
        const FirstLocalOperator          &firstLocalOp;
        const SecondLocalOperator         &secondLocalOp;
        const ThirdLocalOperator          &thirdLocalOp;
        const FourthLocalOperator         &fourthLocalOp;


};

}   // namespace lawa

#include <applications/new_eval_scheme/source/compoundlocaloperator.tcc>

#endif // APPLICATIONS_NEWEVALSCHEME_COMPOUNDLOCALOPERATOR_H
