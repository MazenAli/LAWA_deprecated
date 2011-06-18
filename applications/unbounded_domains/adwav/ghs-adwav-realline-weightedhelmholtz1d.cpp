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
#include <applications/unbounded_domains/parameters/parameters.h>

typedef double T;
using namespace lawa;
using namespace std;

typedef Basis<T,Primal,R,CDF>                                           Basis1D;
typedef Wavelet<T,Primal,R,CDF>                                         WaveletR;

typedef IndexSet<Index1D>::const_iterator const_set_it;

//Operator definitions
typedef WeightedHelmholtzOperator1D<T,Basis1D>                          WeightedHelmholtzOp1D;
typedef WeightedSobolevMidPointPreconditioner1D<T,Basis1D>              WeightedPreconditioner1D;

typedef CompressionWeightedPDE1D<T,Basis1D>                             WeightedCompression1D;

//Righthandsides definitions
typedef RHSWithPeaks1D<T,Basis1D>                                       RhsIntegral1D;
typedef RHS<T,Index1D, RhsIntegral1D, WeightedPreconditioner1D>         Rhs1D;


//Matrix definitions
typedef MapMatrix<T,Index1D,WeightedHelmholtzOp1D,
                  WeightedCompression1D, WeightedPreconditioner1D>      MA;

//APPLY definitions
typedef Parameters<T, Basis1D, WeightedHelmholtzOp1D,
                   WeightedPreconditioner1D>                            Parameters1D;
typedef SYM_WEIGHTED_APPLY_1D<T,Basis1D,Parameters1D,MA>                APPLY1D;

//Algorithm definition
typedef GHS_ADWAV1D<T,Basis1D,APPLY1D,Rhs1D>                            GHS_Adwav;


template <typename T>
Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const Basis1D &basis, const RhsIntegral1D &rhsintegral1d_smooth,
                    const WeightedPreconditioner1D &P);


ExponentialWeightFunction1D<T> expweight;

const int nr = 1;
const T   eta = 2.;

T
u_f(T x) {
    if (nr==1)  return x*x;
    else        return 0.;
}

T
du_f(T x) {
    if (nr==1)  return 2.*x;
    else        return 0.;
}

T
ddu_f(T x) {
    if (nr==1)  return 2.;
    else        return 0.;
}

T
rhs_f(T x) {
    return (-ddu_f(x) - expweight.dalpha(x)*du_f(x) + u_f(x))*expweight.weight(x);
}

T
H1norm() {
    if (nr==1) return sqrt( 0.04608930968270510 + 0.2238051764944967);
}

int main (int argc, char *argv[]) {
    cout.precision(8);
    if (argc != 4) {
        cout << "usage " << argv[0] << " d d_ max_its" << endl; exit(1);
    }
    int d=atoi(argv[1]), d_=atoi(argv[2]);
    int NumOfIterations=atoi(argv[3]);

    int jmin=1;

    expweight.setEta(eta);
    expweight.setSingularPoints();
    Function<T> expweightFct(expweight.weight,expweight.singularPoints);

    T eps=1e-5;
    Basis1D                     basis(d,d_,jmin);
    WeightedHelmholtzOp1D       op(basis,1.,expweightFct,4);
    WeightedPreconditioner1D    P(basis,expweightFct,1);
    WeightedCompression1D       Compr(basis);
    Parameters1D                parameters(basis,op,true,basis.j0);

    DenseVector<Array<T> >      rhs_singularPoints = expweight.singularPoints;
    cout << rhs_singularPoints << endl;
    GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > deltas;
    Function<T>                 rhsFct(rhs_f,rhs_singularPoints);
    RhsIntegral1D               rhsintegral1d(basis, rhsFct, deltas, 42);
    Coefficients<Lexicographical,T,Index1D> f;
    f = initializeRHSVector<T>(basis, rhsintegral1d, P);
    Rhs1D                       F(rhsintegral1d,P,f);

    MA                          A(op,P,Compr);

    APPLY1D                     Apply(parameters,basis,A);

    GHS_Adwav                   ghs_adwav(basis,Apply,F);



    cout << "ADWAV started." << endl;
    ghs_adwav.SOLVE(f.norm(2.),eps,NumOfIterations,H1norm());
    cout << "ADWAV finished." << endl;

    stringstream plot_filename;
    plot_filename << "ghs-adwav-realline-weightedhelmholtz1d-plot_"
                  << "_" << d << "_" << d_ << "_" << jmin << ".dat";
    cout << "Plot of solution started." << endl;
    w_plot<T, Basis1D, WeightedPreconditioner1D>(basis, ghs_adwav.solutions[NumOfIterations-1],
                                                 P, u_f, expweight.sqrtweight, -20., 20.,
                                                 pow2i<T>(-5), plot_filename.str().c_str());
    cout << "Plot of solution finished." << endl;



    return 0;
}


template <typename T>
Coefficients<Lexicographical,T,Index1D>
initializeRHSVector(const Basis1D &basis, const RhsIntegral1D &rhsintegral1d_smooth,
                    const WeightedPreconditioner1D &P)
{
    Coefficients<Lexicographical,T,Index1D> f_smooth;
    IndexSet<Index1D> Lambda;
    for (int k=-60; k<=60; ++k) {
        Lambda.insert(Index1D(basis.j0,k,XBSpline));
    }
    for (int j=basis.j0; j<=6; ++j) {
        for (int k=-60*pow2i<T>(j-basis.j0); k<=60*pow2i<T>(j-basis.j0); ++k) {
            Lambda.insert(Index1D(j,k,XWavelet));
        }
    }

    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        f_smooth[*it] = P(*it)*rhsintegral1d_smooth(*it);
        //cout << *it << ": " << basis.generator((*it).xtype).support((*it).j,(*it).k)
        //     << " " << P(*it) << " " << rhsintegral1d_smooth(*it) << endl;
    }

    return f_smooth;
}
