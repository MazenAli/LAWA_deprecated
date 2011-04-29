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
#include <lawa/lawa.h>
#include <applications/unbounded_domains/referencesolutions/referencesolutions.h>

typedef double T;
using namespace lawa;
using namespace std;

typedef Basis<T,Primal,R,CDF> Basis1D;

typedef flens::DenseVector<flens::Array<T> >                      DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > FullMatT;

//Operator definitions
typedef PDENonConstCoeffOperator1D<T,Basis1D>                     PDEOperator1D;
typedef NoCompression<T,Index1D,Basis1D>                          Compression1D;
typedef WeightedSobolevNormPreconditioner1D<T,Basis1D>            Preconditioner1D;

//Righthandsides definitions
typedef RHSWithPeaks1D<T,Basis1D >                                RhsIntegral1D;
typedef RHS<T,Index1D, RhsIntegral1D, Preconditioner1D>           Rhs;

//MapMatrix definition
//typedef MapMatrixWithZeros<T,Index1D,HelmholtzBilinearForm1D,Compression1D,Preconditioner1D> MA;
typedef MapMatrix<T,Index1D,PDEOperator1D,
                  Compression1D,Preconditioner1D>                 MA;

//Algorithm definition
typedef S_ADWAV<T,Index1D, Basis1D, MA, Rhs> S_Adwav;

const T eta = 2.;

T
weight_f(T x) {
    return exp(-2*eta*fabs(x));
}

T
sqrt_weight_f(T x) {
    return sqrt(weight_f(x));
}

T
zero_f(T x) {
    return 0.;
}

T
u_f(T x) {
    return 100*std::max(exp(x-0.1)-1,0.);
}

T
rhs_f(T x) {
    return weight_f(x) * u_f(x);
}


template <typename T, typename Basis, typename Preconditioner>
void
plot_w(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
       const Preconditioner &P, T (*u)(T), T (*w)(T), T a, T b, T h, const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;

    std::ofstream plotfile(filename);
    for (T x=a; x<=b; x+=h) {
        T appr=0., d_appr = 0.0;
        T exact= u(x);
        for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int j = (*it).first.j, k = (*it).first.k;
            T coeff = (*it).second, prec = P((*it).first);
            appr   += prec * coeff * basis.generator((*it).first.xtype)(x,j,k,0);
        }
        plotfile << x << " " << exact << " " << exact*w(x)
                      << " " << appr  << " " << appr*w(x) << std::endl;
    }
    plotfile.close();
}


int main (int argc, char *argv[]) {
    if (argc != 5) {
        cout << "usage " << argv[0] << " d d_ max_its jmin" << endl; exit(1);
    }
    T c = 1.;
    T contraction = 0.125;
    T threshTol = 0.1, cgTol = 0.001*threshTol, resTol=1e-4;

    int d=atoi(argv[1]), d_=atoi(argv[2]);
    int NumOfIterations=atoi(argv[3]);
    int jmin=atoi(argv[4]);

    cout << "Initializing S-ADWAV, jmin = " << jmin << endl;

    Basis1D basis(d,d_,jmin);

    DenseVectorT             weight_singPts(1); weight_singPts = 0.;
    DenseVectorT             zero_singPts;
    Function<T>              weight(weight_f,weight_singPts);
    Function<T>              zero(zero_f, zero_singPts);
    PDEOperator1D            Bil(basis, weight, zero, zero, 8);

    Preconditioner1D         P(basis,weight,0);

    Compression1D            Compr(basis);

    MA A(Bil,P,Compr);

    DenseVectorT             rhs_singPts(1); rhs_singPts = 0.1;
    Function<T>              rhs(rhs_f,rhs_singPts);
    FullMatT                 rhs_deltas;
    RhsIntegral1D            rhsintegral1d(basis, rhs, rhs_deltas, 16);
    Rhs                      F(rhsintegral1d,P);

    IndexSet<Index1D> InitialLambda;
    InitialLambda.insert(Index1D(jmin,0,XBSpline));

    S_Adwav s_adwav(basis, A, F, contraction, threshTol, cgTol, resTol, NumOfIterations, 2, 1e-2);
    cout << "... finished." << endl;

    Timer time;
    time.start();
    s_adwav.solve_cg(InitialLambda, std::sqrt(558.6000383630328));
    time.stop();
    cout << "S-ADWAV required " << time.elapsed() << " seconds real time" << endl;

    stringstream plot_filename;
    plot_filename << "s-adwav-realline-weightedL2-1d-plot_" << d << "_" << d_
                  << "_" << jmin << ".dat";
    cout << "Plot of solution started." << endl;
    plot_w<T, Basis1D, Preconditioner1D>(basis, s_adwav.solutions[NumOfIterations-1], P, u_f,
                                         sqrt_weight_f, -30., 30., pow2i<T>(-5),
                                         plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;
    return 0;
}

