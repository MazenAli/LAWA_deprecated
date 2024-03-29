/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2014  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

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

#ifndef LAWA_PRECONDITIONERS_PRECONDITIONERS1D_DIAGONALMATRIXPRECONDITIONER1D_H
#define LAWA_PRECONDITIONERS_PRECONDITIONERS1D_DIAGONALMATRIXPRECONDITIONER1D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/settings/enum.h>

namespace lawa {

template <typename T, typename Basis, typename BilinearForm>
class DiagonalMatrixPreconditioner1D
{
    public:
        DiagonalMatrixPreconditioner1D(const BilinearForm &a);

        T
        operator()(XType xtype1, int j1, long k1) const;

        T
        operator()(const Index1D &index) const;

    private:
        const BilinearForm &_a;
};

}   // namespace lawa

#include <lawa/preconditioners/preconditioners1d/diagonalmatrixpreconditioner1d.tcc>

#endif // LAWA_PRECONDITIONERS_PRECONDITIONERS1D_DIAGONALMATRIXPRECONDITIONER1D_H

