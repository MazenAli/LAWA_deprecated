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

#include <iostream>
#include <vector>
#include <lawa/lawa.h>

typedef double T;
using namespace lawa;
using namespace std;

typedef Basis<T,Primal,Interval,Dijkema>                                     Basis1D;

typedef flens::DenseVector<flens::Array<T> >                      DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > FullMatT;

//Operator definitions
typedef WeightedIdentityOperator1D<T, Basis1D>                    WeightedPDEOp1D;

//typedef WeightedSobolevNormPreconditioner1D<T,Basis1D>            NormPreconditioner1D;
//typedef WeightedSobolevMidPointPreconditioner1D<T,Basis1D>        MidPointPreconditioner1D;


/// Non-constant weight function
const T singPt = 0.4;

T
weight_f(T x)
{
    return (x<singPt) ? 1. : -1.;
}



int main (int argc, char *argv[]) {
    if (argc != 5) {
        cout << "usage " << argv[0] << " d d_ j0 J" << endl;
        exit(1);
    }

    int d=atoi(argv[1]), d_=atoi(argv[2]);
    int j0=atoi(argv[3]);
    int J=atoi(argv[4]);
    cout.precision(8);

    Basis1D                  basis(d,d_,j0);
    DenseVectorT             weight_singPts(1);
    weight_singPts = singPt;
    Function<T>              weightFct(weight_f, weight_singPts);
    WeightedPDEOp1D          a(basis, weightFct, 4, 0., 1.);

    int k_middle = (basis.rangeJ(j0).firstIndex()+basis.rangeJ(j0).lastIndex())/2;
    T ref_val = a(XWavelet,j0,k_middle,XWavelet,j0,k_middle);
    std::cout << "Support of reference wavelet: " << basis.generator(XWavelet).support(j0,k_middle)
              << std::endl;


    std::stringstream filename;
    filename << "decay_matrix_entries_" << j0 << "_" << d << "_" << d_  << ".dat";
    std::ofstream file(filename.str().c_str());
    for (int j=j0; j<=J; ++j) {
        T largest_entry_per_level = 0.;
        T largest_entry_per_level_wo_singPt = 0.;
        cout << "Computing level j=" << j << endl;
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            T tmp = a(XWavelet,j0,k_middle,XWavelet,j,k);
            largest_entry_per_level = std::max(largest_entry_per_level,fabs(tmp));
            if (   (basis.generator(XWavelet).support(j,k).l2<=singPt)
                || (basis.generator(XWavelet).support(j,k).l1>=singPt) ) {
                largest_entry_per_level_wo_singPt = std::max(largest_entry_per_level_wo_singPt,
                                                             fabs(tmp));
            }
            else {
                cout << "(" << j << ", " << k << "): contains " << singPt << endl;
            }
        }
        file << j-j0 << " " << largest_entry_per_level
                     << " " << largest_entry_per_level_wo_singPt << endl;
    }
    file.close();
    return 0;
}
