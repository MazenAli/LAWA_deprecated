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

#ifndef LAWA_BOX_TENSORBASIS_H
#define LAWA_BOX_TENSORBASIS_H 1

namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
struct TensorBasis
{
    typedef FirstBasis FirstBasisType;
    typedef SecondBasis SecondBasisType;
    
    TensorBasis(const FirstBasis &_basis1, const SecondBasis &_basis2);

    const FirstBasis &first;
    const SecondBasis &second;
    
    int
    dim(const int J_x, const int J_y) const;
    
    int 
    Jx_max(const int J_x, const int J_y, const int jy) const;
    
    int 
    Jy_max(const int J_x, const int J_y, const int jx) const;
};
    
    
} // namespace lawa

#include <lawa/box/tensorbasis.tcc>


#endif // LAWA_BOX_TENSORBASIS_H