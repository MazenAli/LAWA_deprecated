#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/localrefinement.h>

#define DERIV 0

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;

typedef Integral<Gauss, PrimalBasis, PrimalBasis>                   Integral_Psi_Psi;

int main (int argc, char *argv[]) {
    cout.precision(16);
    if (argc!=6) {
        cout << "Usage: " << argv[0] << " d d_ j0 J bc" << endl;
        return 0;
    }
    int d  = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    int j0 = atoi(argv[3]);
    int J  = atoi(argv[4]);
    bool withDirichletBC = atoi(argv[5]);

    PrimalBasis basis(d,d_,j0);
    if (withDirichletBC) basis.enforceBoundaryCondition<DirichletBC>();

    LocalRefinement<PrimalBasis> LocalRefine(basis,withDirichletBC);
    Integral_Psi_Psi integral(basis,basis);

    Coefficients<Lexicographical,T,Index1D> u, Au, Au2;

    for (long k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        u[Index1D(j0,k,XBSpline)] = (T)rand() / RAND_MAX;
        Au[Index1D(j0,k,XBSpline)] = 0.;
    }

    for (int j1=j0; j1<=J; ++j1) {
        for (long k=basis.rangeJ(j1).firstIndex(); k<=basis.rangeJ(j1).lastIndex(); ++k) {
            u[Index1D(j1,k,XWavelet)] = (T)rand() / RAND_MAX;
            Au[Index1D(j1,k,XWavelet)] = 0.;
        }
    }
    cout << "u = " << u << endl;
    Au2 = Au;

    for (coeff1d_it row=Au.begin(); row!=Au.end(); ++row) {
        T val = 0.;
        for (const_coeff1d_it col=u.begin(); col!=u.end(); ++col) {
            val += integral((*row).first.j,(*row).first.k,(*row).first.xtype,DERIV,
                            (*col).first.j,(*col).first.k,(*col).first.xtype,DERIV) * (*col).second;
        }
        (*row).second = val;
    }
    //cout << "Au =  " << Au << endl;

    Coefficients<Lexicographical,T,Index1D> u_loc_single, Au_loc_single, help;
    for (int j1=j0; j1<=J; ++j1) {
        LocalRefine.reconstruct(u, j1, u_loc_single);
        LocalRefine.reconstruct(Au2, j1, Au_loc_single);
    }

    for (coeff1d_it row=Au2.begin(); row!=Au2.end(); ++row) {
        T val1 = 0.;
        for (const_coeff1d_it col=u_loc_single.begin(); col!=u_loc_single.end(); ++col) {
            val1 += integral((*row).first.j,(*row).first.k,(*row).first.xtype,DERIV,
                             (*col).first.j,(*col).first.k,(*col).first.xtype,DERIV) * (*col).second;
        }
        (*row).second = val1;
    }
    cout << "error =  " << Au-Au2 << endl;

    T max_rel_error = 0.;
    T max_abs_error = 0.;
    for (const_coeff1d_it it=Au.begin(); it!=Au.end(); ++it) {
        max_rel_error = std::max(max_rel_error, fabs(((*it).second-Au2[(*it).first])/(*it).second) );
        max_abs_error = std::max(max_abs_error, fabs(((*it).second-Au2[(*it).first])) );
    }
    cout << "Maximum relative error: " << max_rel_error << endl;
    cout << "Maximum absolute error: " << max_abs_error << endl;
    for (const_coeff1d_it it=Au.begin(); it!=Au.end(); ++it) {
        T rel_val = fabs(((*it).second-Au2[(*it).first])/(*it).second);
        if (fabs(rel_val-max_rel_error)==0) {
            cout << "  Relative error: " << (*it).first << " " << (*it).second << " " <<  Au2[(*it).first] << " " << rel_val << endl;
        }
        T abs_val = fabs(((*it).second-Au2[(*it).first]));
        if (fabs(abs_val-max_abs_error)==0) {
            cout << "  Absolute error: " << (*it).first << " " << (*it).second << " " <<  Au2[(*it).first] << " " << abs_val << endl;
        }
    }


    cout << " *******************************************" << endl << endl;


    for (coeff1d_it row=Au_loc_single.begin(); row!=Au_loc_single.end(); ++row) {
        T val2 = 0.;
        for (const_coeff1d_it col=u_loc_single.begin(); col!=u_loc_single.end(); ++col) {
            val2 += integral((*row).first.j,(*row).first.k,(*row).first.xtype,DERIV,
                             (*col).first.j,(*col).first.k,(*col).first.xtype,DERIV) * (*col).second;
        }
        (*row).second = val2;
    }

    Coefficients<Lexicographical,T,Index1D> Au3;
    for (int j1=J; j1>j0; --j1) {
        IndexSet<Index1D> Lambda;
        for (long k=basis.mra.rangeI(j1-1).firstIndex(); k<=basis.mra.rangeI(j1-1).lastIndex(); ++k) {
            Lambda.insert(Index1D(j1-1,k,XBSpline));
        }
        for (long k=basis.rangeJ(j1-1).firstIndex(); k<=basis.rangeJ(j1-1).lastIndex(); ++k) {
            Lambda.insert(Index1D(j1-1,k,XWavelet));
        }
        LocalRefine.decompose_(Au_loc_single, supp(Au2), help, Au3);
        Au_loc_single = help;
    }
    Au3 += help;

    T max_rel_error2 = 0.;
    T max_abs_error2 = 0.;
    for (const_coeff1d_it it=Au.begin(); it!=Au.end(); ++it) {
        max_rel_error2 = std::max(max_rel_error2, fabs(((*it).second-Au3[(*it).first])/(*it).second) );
        max_abs_error2 = std::max(max_abs_error2, fabs(((*it).second-Au3[(*it).first])) );
    }
    cout << "Maximum relative error: " << max_rel_error2 << endl;
    cout << "Maximum absolute error: " << max_abs_error2 << endl;
    for (const_coeff1d_it it=Au.begin(); it!=Au.end(); ++it) {
        T rel_val = fabs(((*it).second-Au3[(*it).first])/(*it).second);
        if (fabs(rel_val-max_rel_error2)==0) {
            cout << "  Relative error: " << (*it).first << " " << (*it).second << " " <<  Au3[(*it).first] << " " << rel_val << endl;
        }
        T abs_val = fabs(((*it).second-Au3[(*it).first]));
        if (fabs(abs_val-max_abs_error2)==0) {
            cout << "  Absolute error: " << (*it).first << " " << (*it).second << " " <<  Au3[(*it).first] << " " << abs_val << endl;
        }
    }


    return 0;
}


/*
 *
    ofstream file("test.dat");
    for (T x=basis.psi.support(j_row,k_row).l1; x<=basis.psi.support(j_row,k_row).l2; x+=0.00001) {
        T refine = 0.;
        for (const_coeff1d_it col=u_loc_single.begin(); col!=u_loc_single.end(); ++col) {
            refine += (*col).second * basis.mra.phi(x,(*col).first.j,(*col).first.k,0);
        }
        file << x << " " << basis.psi(x,j_row,k_row,0) << " " << refine << endl;
    }
 *
 */
