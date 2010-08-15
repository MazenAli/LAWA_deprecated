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


#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H 1

#include <adaptive/index.h>
#include <adaptive/bilinearform.h>

namespace lawa {

template <typename T, typename Index, typename Basis, typename BilinearForm>
class Preconditioner
{
public:
	T
	operator()(const Index &index) const;

	T
	rescale(const Index &index) const;
};

template <typename T, typename Basis>
class Preconditioner<T,Index1d,Basis,HelmholtzOperator1d<T,Basis> >
{
public:
	T
	operator()(const Index1d &index) const;

	T
	rescale(const Index1d &index) const;

};

}

#include <adaptive/preconditioner.tcc>


#endif /* PRECONDITIONER_H_ */
