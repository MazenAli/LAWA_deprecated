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
typedef HelmholtzOperator1D<T, Basis1D>                     PDEOp;
typedef WeightedHelmholtzOperator1D<T, Basis1D>             WeightedPDEOp;
typedef WeightedSobolevMidPointPreconditioner1D<T,Basis1D>  WeightedPreconditioner;
typedef NoCompression<T,Index1D,Basis1D>                    Compression;

/// Matrix definition
typedef MapMatrix<T,Index1D,PDEOp,Compression,WeightedPreconditioner>         MA;
typedef MapMatrix<T,Index1D,WeightedPDEOp,Compression,WeightedPreconditioner> WeightedMA;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

/// Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff_it;

// Weight definition
const T eta=2;

T
weight(T x) {
    return 2.+sin(x);
/*
    if (fabs(x)>=1.) {
        return exp(-2*eta*fabs(x));
    }
    else {
        return exp(2*eta*(1./16.)*(-5-15*x*x+5*x*x*x*x-x*x*x*x*x*x));
    }
*/
}

T
prec_weight(T x) {
    return 1.;
    //return weight(x);
}

int main (int argc, char *argv[]) {
    if (argc != 6) {
        cout << "usage " << argv[0] << " d d_ j0 k0 J" << endl; exit(1);
    }
    int d =atoi(argv[1]);
    int d_=atoi(argv[2]);
    int j0=atoi(argv[3]);
    int k =atoi(argv[4]);
    int J =atoi(argv[5]);

    Basis1D basis(d,d,j0);
    Coefficients<Lexicographical,T,Index1D> entries;

    cout << basis.mra.phi.support(0,0) << endl;

    ofstream weight_file("weight_function.txt");
    for (T x=-200; x<=200.; x+=0.1) {
        weight_file << x << " " << weight(x) << endl;
    }

    DenseVectorT singPts(2);
    singPts = -1., 1.;
    Function<T> weightFct(weight,singPts);
    Function<T> prec_weightFct(prec_weight,singPts);
    PDEOp                   op(basis,1.);
    WeightedPDEOp           weighted_op(basis,1.,weightFct,100);
    WeightedPreconditioner  P(basis,prec_weightFct,1);
    Compression             Compr(basis);
    MA                      A(op,P,Compr);
    WeightedMA              weighted_A(weighted_op,P,Compr);

    Index1D col(j0,k,XWavelet);
    Support<T> supp_col=basis.generator(col.xtype).support(col.j,col.k);
    // DenseVectorT singsupp_col=basis.generator(col.xtype).singularSupport(col.j,col.k);
    DenseVectorT singsupp_col=basis.psi.optim_singularSupport(col.j,col.k);
    int s_tilde_level = J+1;
    int s_tilde_singsupp = J+1;
    IndexSet<Index1D> LambdaTilde = lambdaTilde1d_WeightedPDE(col, basis, s_tilde_level, j0, j0+J,
                                                              s_tilde_singsupp);
    for (const_set_it row=LambdaTilde.begin(); row!=LambdaTilde.end(); ++row) {
        //cout << col << ", " << *row  << ": " << weighted_A(*row,col) << " " << A(*row,col) << endl;
        entries[(*row)] = weighted_A(*row,col);
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
                else {
                    max_entry_level=std::max(max_entry_level,fabs((*it).second));
                }
            }
        }
        file << s_tilde_level << " " << max_entry_level << " " << max_entry_singsupp << endl;
    }
    file.close();

    ofstream file2("matrixentries_translation.txt");
    for (int m=-50; m<=50; ++m) {
        file2 << m << " " << weighted_A(Index1D(j0,m,XBSpline),Index1D(j0,m,XBSpline)) << endl;
    }
    file2.close();

    return 0;
}

