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

template <typename T, typename Basis, typename Compression, typename Preconditioner>
TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, Compression, Preconditioner>::TensorMatrix2D(const Basis &_basis,
																					   const Preconditioner &_p, Compression &_c)
	: basis(_basis), p(_p), c(_c),
	  dd_x(basis.first), id_x(basis.first), dd_y(basis.second), id_y(basis.second),
	  c_x(basis.first), c_y(basis.second),
	  data_dd_x(dd_x, prec1d, c_x, 4*1024, 1024), data_id_x(id_x, prec1d, c_x, 4*1024, 1024),
	  data_dd_y(dd_y, prec1d, c_y, 4*1024, 1024), data_id_y(id_y, prec1d, c_y, 4*1024, 1024)
{
}

template <typename T, typename Basis, typename Compression, typename Preconditioner>
T
TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, Compression, Preconditioner>::operator()(const Index2D &row_index,
																								 const Index2D &col_index)
{
	typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
	T prec = 1.;
	const_coeff_it it_P_end       = P_data.end();
	const_coeff_it it_row_index   = P_data.find(row_index);
	if (it_row_index != it_P_end) {
	    prec *= (*it_row_index).second;
	}
	else {
	    T tmp = p(row_index);
	    P_data[row_index] = tmp;
	    prec *= tmp;
	}
	it_P_end       = P_data.end();
	const_coeff_it it_col_index   = P_data.find(col_index);
	if (it_col_index != it_P_end) {
		prec *= (*it_col_index).second;
	}
	else {
		T tmp = p(col_index);
		P_data[col_index] = tmp;
		prec *= tmp;
	}
	return prec * (
		  data_dd_x(row_index.index1,col_index.index1) * data_id_y(row_index.index2,col_index.index2) +
		  data_id_x(row_index.index1,col_index.index1) * data_dd_y(row_index.index2,col_index.index2) +
		  data_id_x(row_index.index1,col_index.index1) * data_id_y(row_index.index2,col_index.index2) );
}

template <typename T, typename Basis, typename Compression, typename Preconditioner>
void
TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, Compression, Preconditioner>::clear()
{
	data_dd_x.clear();
	data_id_x.clear();
	data_dd_y.clear();
	data_id_y.clear();
}


}	//namespace lawa
