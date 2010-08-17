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

template <typename T, typename Basis>
void
getSingularPoints(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff, DenseVector<Array<T> > &sing_pts)
{
	typedef typename Basis::BSplineType PrimalSpline;
	typedef typename Basis::WaveletType PrimalWavelet;
	PrimalSpline phi(basis.mra);
	PrimalWavelet psi(basis);

	typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;
	std::list<T> temp;
	for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
		DenseVector<Array<T> > phi_singpts;
		if ((*it).first.xtype == XBSpline) phi_singpts = phi.singularSupport((*it).first.j, (*it).first.k);
		else							   phi_singpts = psi.singularSupport((*it).first.j, (*it).first.k);

		for (int i = phi_singpts.firstIndex(); i <= phi_singpts.lastIndex(); ++i) {
			temp.push_back(phi_singpts(i));
		}
	}
	temp.sort(); temp.unique();
	sing_pts.engine().resize((int)temp.size());
	int i = 1;
	for (typename std::list<T>::const_iterator it = temp.begin(); it != temp.end(); ++it ) {
	    sing_pts(i) = *it; ++i;
	}
}

template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
	 const Preconditioner &P, T (*u)(T), const char* filename)
{
	typedef typename Basis::BSplineType PrimalSpline;
	typedef typename Basis::WaveletType PrimalWavelet;
	typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;

	PrimalSpline phi(basis.mra);
	PrimalWavelet psi(basis);

	std::stringstream PlotFileName;
	PlotFileName << filename << ".dat";
	std::ofstream plotfile(PlotFileName.str().c_str());

	DenseVector<Array<T> > sing_pts;
	getSingularPoints(basis, coeff, sing_pts);

	for (int i=sing_pts.firstIndex(); i<=sing_pts.lastIndex(); ++i) {
		T x = sing_pts(i);
		T appr = 0.0;
		T exact= u(x);
		for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
			int j = (*it).first.j, k = (*it).first.k;
			T coeff = (*it).second, prec = P((*it).first);
			if ((*it).first.xtype == XBSpline) {
				appr  += prec * coeff * phi(x,j,k);
			}
			else {
				appr  += prec * coeff * psi(x,j,k);
			}
		}
		plotfile << x << " " << exact << " " << appr  << " " << fabs(exact-appr) << std::endl;
	}
	plotfile.close();
}


}  // namespace lawa
