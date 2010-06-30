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

#include <lawa/interval/dijkema/dijkema.h>
#include <lawa/interval/dku/dku.h>
#include <lawa/interval/fwt.h>
#include <lawa/interval/primal/primal.h>
#include <lawa/interval/primbs/primbs.h>
#include <lawa/interval/initial_stable_completion.h>
#include <lawa/interval/evaluate.h>
#include <lawa/interval/refinementmatrix.h>
