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

namespace lawa {

template<typename T, typename Basis2D>
SeparableRHS2D<T, Basis2D>::SeparableRHS2D(const Basis2D& _basis, const SeparableFunction2D<T>& _F, int order)
    : basis(_basis), F(_F), phi_x(_basis.first.mra), 
      phi_y(_basis.second.mra), psi_x(_basis.first), psi_y(_basis.second),
      integral_sff_x(phi_x, _F.F_x), integral_sff_y(phi_y, _F.F_y),
      integral_wf_x(psi_x, _F.F_x), integral_wf_y(psi_y, _F.F_y)
{
	integral_sff_x.quadrature.setOrder(order);
	integral_sff_y.quadrature.setOrder(order);
	integral_wf_x.quadrature.setOrder(order);
	integral_wf_y.quadrature.setOrder(order);
}

template<typename T, typename Basis2D>
T
SeparableRHS2D<T, Basis2D>::operator()(XType xtype_x, int j_x, int k_x,
                                     XType xtype_y, int j_y, int k_y) const
{
    T val_x = 0;
    T val_y = 0;
    if(xtype_x == XBSpline){
        val_x = integral_sff_x(j_x, k_x);
    }
    else{
        val_x = integral_wf_x(j_x, k_x);
    }
    
    if(xtype_y == XBSpline){
        val_y = integral_sff_y(j_y, k_y);
    }
    else{
        val_y = integral_wf_y(j_y, k_y);
    }
    
    return val_x * val_y;
}

template<typename T, typename Basis2D>
T
SeparableRHS2D<T, Basis2D>::operator()(const Index2D &index) const
{
	return SeparableRHS2D<T, Basis2D>::operator()(index.index1.xtype, index.index1.j, index.index1.k,
			                                      index.index2.xtype, index.index2.j, index.index2.k);
}

//==========================================================================================//

template<typename T, typename RHS2D>
SumOfRHS2D<T, RHS2D>::SumOfRHS2D(const RHS2D &_rhs1, const RHS2D &_rhs2)
    : rhs1(_rhs1), rhs2(_rhs2)
{
}

template<typename T, typename RHS2D>
T
SumOfRHS2D<T, RHS2D>::operator()(XType xtype_x, int j_x, int k_x,
                                 XType xtype_y, int j_y, int k_y) const
{
	return rhs1(xtype_x, j_x, k_x, xtype_y, j_y, k_y)
		  +rhs2(xtype_x, j_x, k_x, xtype_y, j_y, k_y);
}

template<typename T, typename RHS2D>
T
SumOfRHS2D<T, RHS2D>::operator()(const Index2D &index) const
{
	return rhs1(index.index1.xtype, index.index1.j, index.index1.k,
				index.index2.xtype, index.index2.j, index.index2.k)
		  +rhs2(index.index1.xtype, index.index1.j, index.index1.k,
				index.index2.xtype, index.index2.j, index.index2.k);
}

template<typename T, typename RHS2D>
Coefficients<Lexicographical,T,Index2D>
SumOfRHS2D<T, RHS2D>::operator()(const IndexSet<Index2D> &Lambda) const
{
	typedef typename IndexSet<Index2D>::iterator const_set_it;
	Coefficients<Lexicographical,T,Index2D> ret(Lambda.d,Lambda.d_);
	for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
		T tmp = SumOfRHS2D<T, RHS2D>::operator()(*lambda);
		ret[*lambda] = tmp;
	}
	return ret;
}

} // namespace lawa
