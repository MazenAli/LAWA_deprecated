#include <iostream>
#include <lawa/lawa.h>
#include <applications/unbounded_domains/parameters/parameters.h>

using namespace std;
using namespace lawa;

typedef double T;

/// Thresh bound
const T thresh = 1e-15;


// Basis definitions
const DomainType   domain = R;
const Construction construction = CDF;
typedef Basis<T,Primal,domain,construction>     Basis1D;

/// Operator definitions
typedef WeightedHelmholtzOperator1D<T, Basis1D>             WeightedPDEOp1D;
typedef H1NormPreconditioner1D<T,Basis1D>                   Preconditioner1D;
typedef NoCompression<T,Index1D,Basis1D>                    Compression1D;

/// Matrix definition
typedef MapMatrix<T,Index1D,WeightedPDEOp1D,Compression1D,Preconditioner1D> MA;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

/// Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff_it;

T
coeff(T x) {
    return 2+std::cos(x);

}

int main (int argc, char *argv[]) {
    cout.precision(10);
    if (argc != 6) {
        cout << "usage " << argv[0] << " d d_ j0 k0 J" << endl; exit(1);
    }
    int d =atoi(argv[1]);
    int d_=atoi(argv[2]);
    int j0=atoi(argv[3]);
    int k =atoi(argv[4]);
    int J =atoi(argv[5]);

    Basis1D basis(d,d_,j0);
    Coefficients<Lexicographical,T,Index1D> entries;

    DenseVectorT singPts(40001);
    int pos = 1;
    for (int i=-20000; i<=20000; ++i,++pos) {
        singPts(pos) = i;
    }

    Function<T>         coeffFct(coeff,singPts);
    WeightedPDEOp1D     op(basis,0.,coeffFct,20);
    Preconditioner1D    P(basis);
    Compression1D       Compr(basis);
    MA                  A(op,P,Compr);

    Index1D col(j0,k,XWavelet);
    Support<T> supp_col=basis.generator(col.xtype).support(col.j,col.k);
    // DenseVectorT singsupp_col=basis.generator(col.xtype).singularSupport(col.j,col.k);
    DenseVectorT singsupp_col=basis.psi.optim_singularSupport(col.j,col.k);
    int s_tilde_level = J+1;
    int s_tilde_singsupp = J+1;
    IndexSet<Index1D> LambdaTilde = lambdaTilde1d_WeightedPDE(col, basis, s_tilde_level, j0, j0+J,
                                                              s_tilde_singsupp);
    for (const_set_it row=LambdaTilde.begin(); row!=LambdaTilde.end(); ++row) {
        //cout << col << ", " << *row  << ": "  << A(*row,col) << endl;
        entries[(*row)] = A(*row,col);
        /*
        cout << col  << " " << basis.generator(col.xtype).singularSupport(col.j,col.k) << " "
             << *row << " " << basis.generator((*row).xtype).singularSupport((*row).j,(*row).k)
             << " " << entries[*row] << endl << endl;
        */

    }

    ofstream file("decaymatrixentries.txt");
    for (int s_tilde_level=0; s_tilde_level<=J; ++s_tilde_level) {
        T max_entry_level=0., max_entry_singsupp=0.;
        for (const_coeff_it it=entries.begin(); it!=entries.end(); ++it) {
            if ((*it).first.j-j0!=s_tilde_level) continue;
            Support<T> supp_row=
                    basis.generator((*it).first.xtype).support((*it).first.j,(*it).first.k);
            if (overlap(supp_col,supp_row)>0) {
                if (distance(singsupp_col,supp_row)>=0) {
                    max_entry_singsupp=std::max(max_entry_singsupp,fabs((*it).second));
                    //cout << s_tilde_level << " " << max_entry_singsupp << endl;
                }

                max_entry_level=std::max(max_entry_level,fabs((*it).second));

            }
        }
        file << s_tilde_level << " " << max_entry_level << " " << max_entry_singsupp << endl;
    }
    file.close();

    return 0;
}

