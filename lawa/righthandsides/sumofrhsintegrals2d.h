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

#ifndef LAWA_RIGHTHANDSIDES_SUMOFRHSINTEGRALS2D_H
#define LAWA_RIGHTHANDSIDES_SUMOFRHSINTEGRALS2D_H 1

namespace lawa {

template <typename T, typename RHS2D>
class SumOfTwoRHSIntegrals2D
{
private:
    const RHS2D &rhs1;
    const RHS2D &rhs2;

public:
    SumOfTwoRHSIntegrals2D(const RHS2D &rhs1, const RHS2D &rhs2);

    T
    operator()(XType xtype_x, int j_x, int k_x,
               XType xtpye_y, int j_y, int k_y) const;

    T
    operator()(const Index2D &index) const;
};

template <typename T, typename RHS2D>
class SumOfThreeRHSIntegrals2D
{
private:
    const RHS2D &rhs1;
    const RHS2D &rhs2;
    const RHS2D &rhs3;

public:
    SumOfThreeRHSIntegrals2D(const RHS2D &rhs1, const RHS2D &rhs2, const RHS2D &rhs3);

    T
    operator()(XType xtype_x, int j_x, int k_x,
               XType xtpye_y, int j_y, int k_y) const;

    T
    operator()(const Index2D &index) const;
};



}	//namespace lawa

#include <lawa/righthandsides/sumofrhsintegrals2d.tcc>

#endif  //LAWA_RIGHTHANDSIDES_SUMOFRHSINTEGRALS2D_H
