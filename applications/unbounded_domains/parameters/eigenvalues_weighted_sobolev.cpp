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

typedef Basis<T,Primal,R,CDF>                                     Basis1D;
typedef Wavelet<T,Primal,R,CDF>                                   WaveletR;

typedef flens::DenseVector<flens::Array<T> >                      DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > FullMatT;

//Operator definitions
typedef WeightedIdentityOperator1D<T, Basis1D>                    WeightedPDEOp1D;
typedef NoCompression<T, Index1D, Basis1D>                        Compression;
typedef WeightedSobolevNormPreconditioner1D<T,Basis1D>            NormPreconditioner1D;
typedef WeightedSobolevMidPointPreconditioner1D<T,Basis1D>        MidPointPreconditioner1D;

//MapMatrix definition
typedef MapMatrix<T,Index1D,WeightedPDEOp1D,
                  Compression,NormPreconditioner1D>               MA_norm_prec;
typedef MapMatrix<T,Index1D,WeightedPDEOp1D,
                  Compression,MidPointPreconditioner1D>           MA_midpoint_prec;


template <typename T>
void
computeEV(flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A, T &cB, T &CB);

template <typename T>
void
computeSV(flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A, T &cB, T &CB);

template<typename T>
IndexSet<Index1D>
LambdaForEigenvalues(int jmin, int jmax, const Basis1D &basis, T radius);

typedef IndexSet<Index1D>::const_iterator const_set_it;



/// weight function
const T eta=2.;
T
weight_f(T x)
{
    return exp(-2*eta*fabs(x));
}

/// Non-constant weight function
T
smoothweight_f(T x)
{
    if (fabs(x)>=1.) {
        return exp(-2*eta*fabs(x));
    }
    else {
        return exp(2*eta*(1./16.)*(-5-15*x*x+5*x*x*x*x-x*x*x*x*x*x));
    }
}





int main (int argc, char *argv[]) {
    if (argc != 6) {
        cout << "usage " << argv[0] << " d d_ radius jmin max_level" << endl;
        exit(1);
    }

    int d=atoi(argv[1]), d_=atoi(argv[2]);
    T   radius    =atof(argv[3]);
    int jmin =atoi(argv[4]);
    int max_level =atoi(argv[5]);
    cout.precision(8);

    ofstream plotfile("test.txt");
    for (T x=-5.; x<=5.; x+=0.1) {
        plotfile << x << " " << weight_f(x) << " " << smoothweight_f(x) << endl;
    }
    plotfile.close();

    Basis1D                  basis(d,d_,jmin);
    DenseVectorT             weight_singPts(3);
    weight_singPts = -1.,0.,1.;
    Function<T>              weightFct(smoothweight_f, weight_singPts);
    WeightedPDEOp1D          Bil(basis, weightFct, 8, 0., 1.);
    //NormPreconditioner1D     NormP(basis,weight,0);
    MidPointPreconditioner1D MidPointP(basis,weightFct,0);
    Compression              compression(basis);
    //MA_norm_prec             A_norm_prec(Bil,NormP,compression);
    MA_midpoint_prec         A_midpoint_prec(Bil,MidPointP,compression);

    cout << "Chosen parameters: d=" << d << ", d_=" << d_  << ", radius="<< radius
         << ", jmin=" << jmin << ", max_level=" << max_level << endl;



    T cB_norm, CB_norm, cB_midpoint, CB_midpoint;
    FullMatT A_dense;
    std::stringstream filename;
    filename << "eigenvalues_weighted_sobolev_" << jmin << "_" << d << "_" << d_  << ".dat";
    std::ofstream file(filename.str().c_str());
    for (int jmax=jmin; jmax<=max_level; jmax+=1) {
        for (T r=1.; r<=radius; r+=1.) {
            IndexSet<Index1D> Lambda = LambdaForEigenvalues(jmin, jmax, basis, r);
            int N = Lambda.size();
            cout << "Size of Lambda: " << N << endl;
            /*
            cout << "Norm preconditioning ..." << endl;
            flens::SparseGeMatrix<CRS<T,CRS_General> > A_norm_prec_sparse(N,N);
            toFlensSparseMatrix(A_norm_prec, Lambda, Lambda, A_norm_prec_sparse);
            spy(A_norm_prec_sparse,"A_norm_prec");
            densify(cxxblas::NoTrans,A_norm_prec_sparse,A_dense);
            cB_norm=0., CB_norm=0.;
            computeEV(A_dense, cB_norm, CB_norm);
            cout             << " " << jmax << " " << r
                             << " " << cB_norm << " " << " " << CB_norm << endl;
            */

            cout << "Midpoint preconditioning..." << endl;
            flens::SparseGeMatrix<CRS<T,CRS_General> > A_midpoint_prec_sparse(N,N);
            toFlensSparseMatrix(A_midpoint_prec, Lambda, Lambda, A_midpoint_prec_sparse);
            densify(cxxblas::NoTrans,A_midpoint_prec_sparse,A_dense);
            spy(A_midpoint_prec_sparse,"A_midpoint_prec");
            cB_midpoint=0., CB_midpoint=0.;
            computeEV(A_dense, cB_midpoint, CB_midpoint);
            cout             << " " << jmax << " " << r
                             << " " << cB_midpoint << " " << " " << CB_midpoint << endl << endl;

            file << r << " " << jmax-jmin << " " << cB_norm << " " << CB_norm << " "
                 << cB_midpoint << " " << CB_midpoint << endl;
        }
        file << endl;
    }
    file.close();
    return 0;
}

template <typename T>
void
computeEV(flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A, T &cB, T &CB) {
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > U(A.numRows(),A.numRows()),
                                                               V(A.numCols(),A.numCols());
    int N = A.numRows();
    DenseVector<Array<T> > wr(N), wi(N);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > vl,vr;
    ev(false, false, A, wr, wi, vl, vr);
    cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
    for (int i=1; i<=wr.lastIndex(); ++i) {
        cB = std::min(cB,wr(i));
        CB = std::max(CB,wr(i));
    }
}

template <typename T>
void
computeSV(flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &A, T &cB, T &CB) {
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > U(A.numRows(),A.numRows()),
                                                               V(A.numCols(),A.numCols());
    DenseVector<Array<T> > s(A.numCols());
    int iterations = svd(A,s,U,V);
    CB = s(s.firstIndex());
    cB = s(s.lastIndex());
}

template<typename T>
IndexSet<Index1D>
LambdaForEigenvalues(int jmin, int jmax, const Basis1D &basis, T radius)
{
    const BSpline<T,Primal,R,CDF> phi = basis.mra.phi;
    const Wavelet<T,Primal,R,CDF> psi = basis.psi;
    IndexSet<Index1D> Lambda;
    int k_left, k_right;

    for (int j=jmin; j<=jmax; ++j) {
        k_left = std::floor(-pow2i<T>(j)*radius-psi.support(0,0).l2);
        k_right = std::ceil(pow2i<T>(j)*radius-psi.support(0,0).l1);
        for (int k=k_left; k<=k_right; ++k) {
            Lambda.insert(Index1D(j,k,XWavelet));
        }
    }

    k_left  = int(std::floor(-pow2i<T>(jmin)*radius-phi.support(0,0).l2));
    k_right = int(std::ceil(  pow2i<T>(jmin)*radius-phi.support(0,0).l1));
    for (int k=k_left; k<=k_right; ++k) {
        Lambda.insert(Index1D(jmin,k,XBSpline));
    }

    return Lambda;
}
