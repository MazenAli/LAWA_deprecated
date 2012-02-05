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

#ifndef  APPLICATIONS_NEWEVALSCHEME_MULTITREEOPERATIONS_H
#define  APPLICATIONS_NEWEVALSCHEME_MULTITREEOPERATIONS_H 1

#include <boost/functional/hash.hpp>

#include <lawa/settings/enum.h>
#include <iostream>

namespace lawa {

template <typename Index, typename Basis>
void
extendMultiTree2(const Basis &basis, const Index &index2d, const int offset, IndexSet<Index> &Lambda);


}

#include <applications/new_eval_scheme/source/multitreeoperations.tcc>

#endif // APPLICATIONS_NEWEVALSCHEME_MULTITREEOPERATIONS_H
