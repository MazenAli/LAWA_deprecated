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
typedef HelmholtzOperator1D<T, Basis1D>                     PDEOp1D;
//typedef H1NormPreconditioner1D<T,Basis1D>                 Preconditioner1D;
typedef NoPreconditioner<T,Index1D>                         Preconditioner1D;
typedef NoCompression<T,Index1D,Basis1D>                    Compression1D;

/// Matrix definition
typedef MapMatrix<T,Index1D,WeightedPDEOp1D,Compression1D,Preconditioner1D> MA_Weighted;
typedef MapMatrix<T,Index1D,PDEOp1D,Compression1D,Preconditioner1D>         MA_Const;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

/// Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff_it;

T
coeff(T x) {
    return 2.+log(200000+x);
    //return 2.+exp(-0.01*x*x);
}



IndexSet<Index1D>
computeLambda(const Basis1D &basis, T a, T b, int J_plus);

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
    Coefficients<Lexicographical,T,Index1D> entries_weighted;
    Coefficients<Lexicographical,T,Index1D> entries_const;


    DenseVectorT singPts1(32001);
    int pos = 1;
    for (int i=-16000; i<=16000; ++i,++pos) {
        singPts1(pos) = (T)i;
    }
    DenseVectorT singPts2(16001);
    pos = 1;
    for (int i=-8000; i<=8000; ++i,++pos) {
        singPts2(pos) = (T)i;
    }

    Function<T>         coeffFct1(coeff,singPts1);
    Function<T>         coeffFct2(coeff,singPts2);
    WeightedPDEOp1D     weighted_op1(basis,1.,coeffFct1,20);
    WeightedPDEOp1D     weighted_op2(basis,1.,coeffFct2,20);
    PDEOp1D             op(basis,1.);
    //Preconditioner1D    P(basis);
    Preconditioner1D    P;
    Compression1D       Compr(basis);
    MA_Weighted         A_weighted(weighted_op1,P,Compr);
    MA_Const            A_const(op,P,Compr);

    IndexSet<Index1D> LambdaTilde;

    ofstream file("tmp.txt");
    for (T x=-20.; x<=20.; x+=0.01) {
        file << x << " " << coeff(x) << endl;
    }
    file.close();



    Index1D col(j0,k,XWavelet);
    Support<T> supp_col=basis.generator(col.xtype).support(col.j,col.k);
    cout << supp_col << endl;

    DenseVectorT singsupp_col=basis.psi.optim_singularSupport(col.j,col.k);
    int s_tilde_level = J+1;
    int s_tilde_singsupp = J+1;

    LambdaTilde = computeLambda(basis,supp_col.l1,supp_col.l2,basis.j0+J);

    for (const_set_it row=LambdaTilde.begin(); row!=LambdaTilde.end(); ++row) {
        entries_weighted[(*row)] = A_weighted(*row,col);
        Support<T> supp_row=basis.generator((*row).xtype).support((*row).j,(*row).k);
        //cout << col << ", " << *row  << ": " << supp_col << " " << supp_row << " " <<  A_weighted(*row,col) << endl;

    }
    cout << endl << endl;

    LambdaTilde = lambdaTilde1d_PDE(col, basis, s_tilde_level, j0, j0+J,false);
    for (const_set_it row=LambdaTilde.begin(); row!=LambdaTilde.end(); ++row) {
        entries_const[(*row)] = A_const(*row,col);
    }

    ofstream file1("decaymatrixentries_const.txt");
    for (int s_tilde_level=0; s_tilde_level<=J; ++s_tilde_level) {
        T max_entry_level=0., max_entry_singsupp=0.;
        for (const_coeff_it it=entries_const.begin(); it!=entries_const.end(); ++it) {
            if ((*it).first.j-j0!=s_tilde_level) continue;
            Support<T> supp_row=
                    basis.generator((*it).first.xtype).support((*it).first.j,(*it).first.k);
            if (overlap(supp_col,supp_row)>0) {

                max_entry_level=std::max(max_entry_level,fabs((*it).second));

            }
        }
        file1 << s_tilde_level << " " << max_entry_level << endl;
    }
    file1.close();


    ofstream file2("decaymatrixentries_weighted.txt");
    for (int s_tilde_level=0; s_tilde_level<=J; ++s_tilde_level) {
        T max_entry_level=0., max_entry_singsupp=0.;
        for (const_coeff_it it=entries_weighted.begin(); it!=entries_weighted.end(); ++it) {
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
        file2 << s_tilde_level << " " << max_entry_level << " " << max_entry_singsupp << endl;
    }
    file2.close();

    return 0;
}

IndexSet<Index1D>
computeLambda(const Basis1D &basis, T a, T b, int J_plus)
{
    IndexSet<Index1D> ret;
    T l1, l2;
    l1 = basis.mra.phi.support(0,0).l1, l2 = basis.mra.phi.support(0,0).l2;
    int k_left =  std::floor(float(pow2i<T>(basis.j0)*a-l2));
    int k_right = std::ceil(float(pow2i<T>(basis.j0)*b-l1));
    for (int k=k_left; k<=k_right; ++k) {
        ret.insert(Index1D(basis.j0,k,XBSpline));
    }

    l1 = basis.psi.support(0,0).l1, l2 = basis.psi.support(0,0).l2;
    for (int j=basis.j0; j<=J_plus; ++j) {
        int k_left =  std::floor(float(pow2i<T>(j)*a-l2));
        int k_right = std::ceil(float(pow2i<T>(j)*b-l1));
        for (int k=k_left; k<=k_right; ++k) {
            ret.insert(Index1D(j,k,XWavelet));
        }
    }
    return ret;
}
