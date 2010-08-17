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

template <typename T>
void
plotCoeff(const Coefficients<AbsoluteValue,T,Index1D > &coeff, const Basis<T,Primal,R,CDF> &basis, const char* filename)
{
    typedef typename Coefficients<AbsoluteValue,T,Index1D>::const_iterator const_it;
    if (coeff.size() == 0) {
        return;
    }

    std::stringstream gpsFilename;
    gpsFilename << filename << ".gps";
    std::ofstream gps(gpsFilename.str().c_str());

    gps << "reset" << std::endl;
    gps << "set terminal postscript eps enh color; set output '" << filename << ".eps'" << std::endl;
    gps << "set palette color; set colorbox vertical" << std::endl;

    BSpline<T,Primal,R> phi(basis.mra);
    Wavelet<T,Primal,R> psi(basis);

    T maxCoeffSca = -1.0;
    T maxCoeffWav = -1.0;
    int j = (*coeff.begin()).second.j, k = (*coeff.begin()).second.k;
    int j0 = j;
    int J  = j;
    T a_sca = 5000.0, a_wav = 5000.0;
    T b_sca = -5000.0, b_wav = -5000.0;
    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
    	j = (*it).second.j; k = (*it).second.k;
    	j0 = std::min(j0, j);
        J  = std::max(J, j);
        if ((*it).second.xtype == XBSpline) {
        	maxCoeffSca = std::max(maxCoeffSca, fabs((*it).first));
        	a_sca = std::min(a_sca, phi.support(j,k).l1);
            b_sca = std::max(b_sca, phi.support(j,k).l2);
        }
        else {
        	maxCoeffWav = std::max(maxCoeffWav, fabs((*it).first));
        	a_wav = std::min(a_wav, psi.support(j,k).l1);
        	b_wav = std::max(b_wav, psi.support(j,k).l2);
        }
    }
    T maxCoeff = std::max(maxCoeffWav,maxCoeffSca);
    T l1_sca = phi.support(0,0).l1, l2_sca = phi.support(0,0).l2;
    T l1_wav = psi.support(0,0).l1, l2_wav = psi.support(0,0).l2;

    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        T lineWidth = 0.1;
        T ctr, fromX, toX, fromY, toY;
        T color = 0.0;

        if ((*it).second.xtype==XBSpline) {
        	int k1 = ceil(pow2i<T>((*it).second.j)*a_sca - l1_sca), k2 = floor(pow2i<T>((*it).second.j)*b_sca - l2_sca);
        	int N = k2 - k1 + 1;
            fromX = a_sca + ((*it).second.k-k1)*(b_sca-a_sca)/(T)N;
            toX   = a_sca + ((*it).second.k-k1+1)*(b_sca-a_sca)/(T)N;

            fromY = (*it).second.j-1.5;
            toY   = (*it).second.j-0.5;
            color = fabs((*it).first) / maxCoeffSca;
        }

        else {
        	long int k1 = ceil(pow2i<T>((*it).second.j)*a_wav - l1_wav), k2 = floor(pow2i<T>((*it).second.j)*b_wav - l2_wav);
        	long int N = k2 - k1 + 1;
        	fromX = a_wav + ((*it).second.k-k1)*(b_wav-a_wav)/(T)N;
        	toX   = fromX + std::max((b_wav-a_wav)/N,0.01);  //was 0.05

        	fromY = (*it).second.j-0.5;
            toY   = (*it).second.j+0.5;
            color = fabs((*it).first) / maxCoeffWav;

        }

        gps << "set object rectangle from " << fromX << ", " << fromY
            << " to " << toX << "," << toY << " fc rgb";
        if (color > 0.5) gps << " 'black' ";
        else if ((0.5 >= color) && (color > 0.25)) gps << " 'purple' ";
        else if ((0.25 >= color) && (color > 0.125)) gps << " 'magenta' ";
        else if ((0.125 >= color) && (color > 0.0625)) gps << " 'red' ";
        else if ((0.0625 >= color) && (color > 0.03125)) gps << " 'orangered' ";
        else if ((0.03125 >= color) && (color > 0.015625)) gps << " 'orange' ";
        else if ((0.015625 >= color) && (color > 0.0078125)) gps << " 'yellow' ";
        else gps << " 'grey' ";
        gps << " linewidth " << 0.1 << " fillstyle solid" << std::endl;

    }

    gps << "set xrange["<< std::min(a_sca, a_wav) <<":"<< std::max(b_sca, b_wav) <<"]" << std::endl;
    gps << "set yrange[" << j0-1.5 << ":" << J+0.5 << "]" << std::endl;
    //gps << "set xtics 1" << endl;
    gps << "set ytics ('" << j0 << "' " << j0-1;
    gps << ", '" << j0 << "' " << j0;
    for (int j = j0+1; j <= J; ++j) {
        gps << ", '" << j << "' " << j;
    }
    gps << ")" << std::endl;
    gps << "plot " << j0-1.5 << " with lines linecolor rgb 'black' notitle" << std::endl;
    gps << "reset; set terminal pop" << std::endl;
    gps.close();
}

}  // namespace lawa
