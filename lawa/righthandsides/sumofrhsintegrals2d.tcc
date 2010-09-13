/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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


namespace lawa {



template<typename T, typename RHS2D>
SumOfTwoRHSIntegrals2D<T, RHS2D>::SumOfTwoRHSIntegrals2D(const RHS2D &_rhs1, const RHS2D &_rhs2)
    : rhs1(_rhs1), rhs2(_rhs2)
{
}

template<typename T, typename RHS2D>
T
SumOfTwoRHSIntegrals2D<T, RHS2D>::operator()(XType xtype_x, int j_x, int k_x,
                                 XType xtype_y, int j_y, int k_y) const
{
	return rhs1(xtype_x, j_x, k_x, xtype_y, j_y, k_y)
		  +rhs2(xtype_x, j_x, k_x, xtype_y, j_y, k_y);
}

template<typename T, typename RHS2D>
T
SumOfTwoRHSIntegrals2D<T, RHS2D>::operator()(const Index2D &index) const
{
    return rhs1(index.index1.xtype, index.index1.j, index.index1.k,
                index.index2.xtype, index.index2.j, index.index2.k)
          +rhs2(index.index1.xtype, index.index1.j, index.index1.k,
                index.index2.xtype, index.index2.j, index.index2.k);
}



template<typename T, typename RHS2D>
SumOfThreeRHSIntegrals2D<T, RHS2D>::SumOfThreeRHSIntegrals2D(const RHS2D &_rhs1, const RHS2D &_rhs2,
	                                                         const RHS2D &_rhs3)
    : rhs1(_rhs1), rhs2(_rhs2), rhs3(_rhs3)
{
}

template<typename T, typename RHS2D>
T
SumOfThreeRHSIntegrals2D<T, RHS2D>::operator()(XType xtype_x, int j_x, int k_x,
                                 XType xtype_y, int j_y, int k_y) const
{
	return rhs1(xtype_x, j_x, k_x, xtype_y, j_y, k_y)
		  +rhs2(xtype_x, j_x, k_x, xtype_y, j_y, k_y)
		  +rhs3(xtype_x, j_x, k_x, xtype_y, j_y, k_y);
}

template<typename T, typename RHS2D>
T
SumOfThreeRHSIntegrals2D<T, RHS2D>::operator()(const Index2D &index) const
{
    return rhs1(index.index1.xtype, index.index1.j, index.index1.k,
                index.index2.xtype, index.index2.j, index.index2.k)
          +rhs2(index.index1.xtype, index.index1.j, index.index1.k,
                index.index2.xtype, index.index2.j, index.index2.k)
          +rhs3(index.index1.xtype, index.index1.j, index.index1.k,
                index.index2.xtype, index.index2.j, index.index2.k);
}


}	//namespace lawa