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

#include <lawa/flensforlawa.h>

#include <lawa/bspline.h>
#include <lawa/enum.h>
#include <lawa/function.h>
#include <lawa/integrals.h>
#include <lawa/interval/interval.h>
#include <lawa/math/math.h>
#include <lawa/mra.h>
#include <lawa/param.h>
#include <lawa/periodic/periodic.h>
#include <lawa/quadrature.h>
#include <lawa/realline/realline.h>
#include <lawa/support.h>
#include <lawa/wavelet.h>
