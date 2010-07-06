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

#ifndef BOXINDEX_H
#define BOXINDEX_H 1

namespace lawa {
  
template<typename Basis>    
class BoxIndex
{
    private:
        const Basis& basis;
        const int J_x, J_y;
        
        const int offsetIx;
        const int offsetIy;
        const int offsetJx;
        const int offsetJy;
    
    public:
        BoxIndex(const Basis& _basis, const int _J_x, const int _J_y);
        
        int
        operator()(bool XisSpline, int jx, int kx,
                   bool YisSpline, int jy, int ky) const;
};
    
    
} // namespace lawa

#include <lawa/box/boxindex.tcc>

#endif // BOXINDEX_H