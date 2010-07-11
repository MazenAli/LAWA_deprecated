/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

template <typename T,Construction Cons>
T
evaluate(const Coefficient<Uniform,T,Cons> &coeff, T x, int deriv)
{
    const Basis<T,Primal,Interval,Cons> &basis = coeff.basis;
    int j0=coeff.j0, J=coeff.J;
    T val = 0;
    for (int k=basis.rangeI(j0).firstIndex(); k<=basis.rangeI(j0).lastIndex(); ++k) {
        val += coeff(k) * basis.phi(x,j0,k,deriv);
    }

    for (int j=j0; j<=J-1; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            val += coeff(basis.rangeI(j).lastIndex() + k) * basis.psi(x, j, k, deriv);
        }
    }
    return val;
}

template <typename T,Construction Cons>
T
evaluate(const Coefficient<Lexicographical,T,Cons> &coeff, T x, int deriv)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    const Basis<T,Primal,Interval,Cons> &basis = coeff.basis;
    BSpline<T,Primal,Interval,Primbs> phi(basis.mra,deriv);
    Wavelet<T,Primal,Interval,Cons> psi(basis,deriv);
    
    T val = 0;
    if (coeff.size()) {
        for (const_it it = coeff.begin(); it != coeff.end(); it++) {
            if ((*it).first.xtype==XBSpline) {
                val += (*it).second * phi(x, (*it).first.j, (*it).first.k);
            } else if ((*it).first.xtype==XWavelet) {
                val += (*it).second * psi(x, (*it).first.j, (*it).first.k);
            } else {
                assert(0);
            }
        }
    }
    return val;
}

template <typename T,Construction Cons>
DenseVector<Array<T> >
evaluate(const Coefficient<Lexicographical,T,Cons> &coeff,
         const DenseVector<Array<T> > &x, int deriv)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    const Basis<T,Primal,Interval,Cons> &basis = coeff.basis;
    BSpline<T,Primal,Interval,Primbs> phi(basis.mra,deriv);
    Wavelet<T,Primal,Interval,Cons> psi(basis,deriv);
    
    DenseVector<Array<T> > val(x.range());
    if (coeff.size()) {
        for (const_it it = coeff.begin(); it != coeff.end(); it++) {
            if ((*it).first.xtype==XBSpline) {
                val += (*it).second * phi(x, (*it).first.j, (*it).first.k);
            } else if ((*it).first.xtype==XWavelet) {
                val += (*it).second * psi(x, (*it).first.j, (*it).first.k);
            } else {
                assert(0);
            }
        }
    }
    return val;
}

template <typename T,Construction Cons>
T
evaluate(const Coefficient<AbsoluteValue,T,Cons> &coeff, T x, int deriv)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    const Basis<T,Primal,Interval,Cons> &basis = coeff.basis;
    
    T val = 0;
    if (coeff.size()) {
        for (const_it it = coeff.begin(); it != coeff.end(); it++) {
            if ((*it).second.xtype==XBSpline) {
                val += (*it).first * basis.phi(x, (*it).second.j, (*it).second.k, deriv);
            } else if ((*it).second.xtype==XWavelet) {
                val += (*it).first * basis.psi(x, (*it).second.j, (*it).second.k, deriv);
            } else {
                assert(0);
            }
        }
    }
    return val;
}

} // namespace lawa