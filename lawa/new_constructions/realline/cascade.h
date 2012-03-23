/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008-2011 Sebastian Kestler, Kristina Steih,
                         Alexander Stippler, Mario Rometsch.

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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_CASCADE_H
#define LAWA_CONSTRUCTIONS_REALLINE_CASCADE_H 1

#include <lawa/constructions/bspline.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template <typename X, typename Y>
void
evalAtDyadicGrid_Cascade(const DenseVector<X> &sf, int J,
                         DenseVector<Y> &scaling);

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_CASCADE_H

