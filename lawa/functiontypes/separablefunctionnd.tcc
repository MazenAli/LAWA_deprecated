/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

template<typename T>
SeparableFunction2D<T>::SeparableFunction2D(Function<T> _F_x, Function<T> _F_y)
    : F_x(_F_x), F_y(_F_y)
{
}

template<typename T>
SeparableFunction2D<T>::SeparableFunction2D(T (*_f_x)(T),
        const DenseVector<Array<T> > &_singularPts_x,
        T (*_f_y)(T), const DenseVector<Array<T> > &_singularPts_y)
    : F_x(Function<T>(_f_x, _singularPts_x)), F_y(Function<T>(_f_y, _singularPts_y))
{
}

template<typename T>
T
SeparableFunction2D<T>::operator()(T x, T y) const
{
    return F_x(x) * F_y(y);
}



template<typename T>
SeparableFunction3D<T>::SeparableFunction3D(Function<T> _F_x, Function<T> _F_y, Function<T> _F_z)
    : F_x(_F_x), F_y(_F_y), F_z(_F_z)
{
}

template<typename T>
SeparableFunction3D<T>::SeparableFunction3D(T (*_f_x)(T), const DenseVector<Array<T> > &_singularPts_x,
											T (*_f_y)(T), const DenseVector<Array<T> > &_singularPts_y,
											T (*_f_z)(T), const DenseVector<Array<T> > &_singularPts_z)
    : F_x(Function<T>(_f_x, _singularPts_x)), F_y(Function<T>(_f_y, _singularPts_y)), F_z(Function<T>(_f_z, _singularPts_z))
{
}

template<typename T>
T
SeparableFunction3D<T>::operator()(T x, T y, T z) const
{
    return F_x(x) * F_y(y) * F_z(z);
}


} //  namespace lawa
