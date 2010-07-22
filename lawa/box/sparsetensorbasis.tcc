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
 
namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
SparseTensorBasis<FirstBasis, SecondBasis>::SparseTensorBasis(const FirstBasis &_basis1, 
                                                              const SecondBasis &_basis2)  
    : first(_basis1), second(_basis2)
{
}

template<typename FirstBasis, typename SecondBasis>
int
SparseTensorBasis<FirstBasis, SecondBasis>::dim(const int J_x, const int J_y) const
{
    int d = first.mra.cardI(first.j0) * second.mra.cardI(J2_max(J_x, J_y, first.j0 - 1));
    for(int jx = first.j0; jx <= J1_max(J_x, J_y, second.j0-1) - 1; ++jx){
        d += first.cardJ(jx)*second.mra.cardI(second.j0);
        for(int jy = second.j0; jy <= J2_max(J_x,J_y, jx) - 1; ++jy){
            d += first.cardJ(jx) * second.cardJ(jy);
        }
    }
    return d;
}

template<typename FirstBasis, typename SecondBasis>
int 
SparseTensorBasis<FirstBasis, SecondBasis>::J1_max(const int J_x, const int J_y, const int jy) const
{
    return std::min(J_x, std::max(J_x + second.j0 - 2, J_y + first.j0 - 2) - jy + 1);
} 

template<typename FirstBasis, typename SecondBasis>
int 
SparseTensorBasis<FirstBasis, SecondBasis>::J2_max(const int J_x, const int J_y, const int jx) const
{
    return std::min(J_y, std::max(J_x + second.j0 - 2, J_y + first.j0 - 2)- jx  + 1);
}
    
} // namespace lawa