#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;

///  Typedefs for problem components:
///     Multiwavelet basis over an interval
typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;

typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;


int main(int argc, char*argv[])
{
    cout.precision(8);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d  = atoi(argv[1]);
    int j0 = atoi(argv[2]);
    int J  = atoi(argv[3]);
    bool sparsetree = true;
    T x_sing = 1./3.;

    /// Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, j0);     // For L2_orthonormal and special MW bases
    //PrimalBasis basis(d, d, j0);     // For biorthogonal wavelet bases
    if (d>1) basis.enforceBoundaryCondition<DirichletBC>();

    Coefficients<Lexicographical,T,Index1D> tree;

    Index1D lambda(J,1,XWavelet);
    for (int k=basis.rangeJ(J).firstIndex(); k<=basis.rangeJ(J).lastIndex(); ++k) {
        Support<T> supp = basis.psi.support(J,k);
        if (supp.l1 < x_sing && x_sing < supp.l2) {
            lambda.k = k;
            break;
        }
    }
    cout << "lambda = " << lambda << endl;

    completeMultiTree(basis, lambda, tree, sparsetree);

    cout << "Singular point: " << x_sing << endl;
    Support<T> covered_interval(0.,1.);
    for (int j=J; j>=j0; --j) {
        int k_min=basis.rangeJ(j).lastIndex(), k_max=basis.rangeJ(j).firstIndex();
        T x_min = 1., x_max = 0.;
        cout << "************ Wavelet Level " << j << " *************" << endl;
        for (const_coeff1d_it it=tree.begin(); it!=tree.end(); ++it) {
            if ((*it).first.j==j && (*it).first.xtype==XWavelet) {
                int k = (*it).first.k;
                Support<T> supp = basis.psi.support(j,(*it).first.k);
                k_min = std::min(k_min,k);
                k_max = std::max(k_max,k);
                x_min = std::min(x_min,supp.l1);
                x_max = std::max(x_max,supp.l2);
                cout << "   " << (*it).first << " " << basis.psi.support(j,k) << endl;
            }
        }
        cout << "Index range: [" << k_min << ", " << k_max << "], covered interval: ["
             << x_min << ", " << x_max << "]" << endl;
        if (!sparsetree) {
            if (j==J) {
                covered_interval.l1 = x_min; covered_interval.l2 = x_max;
            }
            else {
                for (int k=k_min-1; k>=std::max(k_min-6,basis.rangeJ(j).firstIndex()); --k) {
                    if (overlap(covered_interval,basis.psi.support(j,k))>0) {
                        cout << "Error: Index (" << j << ", " << k << ") is missing: " << endl;
                        cout << covered_interval << " " << basis.psi.support(j,k) << endl;
                    }
                }
                for (int k=k_max+1; k<=std::min(k_max+6,basis.rangeJ(j).lastIndex()); ++k) {
                    if (overlap(covered_interval,basis.psi.support(j,k))>0) {
                        cout << "Error: Index (" << j << ", " << k << "is missing: " << endl;
                        cout << covered_interval << " " << basis.psi.support(j,k) << endl;
                    }
                }
                covered_interval.l1 = x_min; covered_interval.l2 = x_max;
            }
        }
        cout << "***************************************" << endl << endl;

        if (j==j0) {
            cout << "************ Scaling Level " << j << " *************" << endl;
            for (const_coeff1d_it it=tree.begin(); it!=tree.end(); ++it) {
                if ((*it).first.xtype==XBSpline) {
                    int k = (*it).first.k;
                    Support<T> supp = basis.mra.phi.support(j,(*it).first.k);
                    k_min = std::min(k_min,k);
                    k_max = std::max(k_max,k);
                    x_min = std::min(x_min,supp.l1);
                    x_max = std::max(x_max,supp.l2);
                    cout << "   " << (*it).first << " " << basis.mra.phi.support(j,k) << endl;
                }
            }
        }
    }

}
