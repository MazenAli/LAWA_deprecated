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
const T leftbound=-1.;   //for realline constructions, we only consider wavelets with support intersecting [left,right]
const T rightbound=1.;
typedef Basis<T,Primal,domain,construction>     Basis1D;

/// Operator definitions
typedef HelmholtzOperator1D<T, Basis1D>                     PDEOp;
typedef H1NormPreconditioner1D<T,Basis1D>                   Preconditioner;
typedef WeightedHelmholtzOperator1D<T, Basis1D>             WeightedPDEOp;
typedef WeightedSobolevMidPointPreconditioner1D<T,Basis1D>  WeightedPreconditioner;
typedef NoCompression<T,Index1D,Basis1D>                    Compression;

/// Matrix definition
typedef MapMatrix<T,Index1D,PDEOp,Compression,Preconditioner>                 MA;
typedef MapMatrix<T,Index1D,WeightedPDEOp,Compression,WeightedPreconditioner> WeightedMA;

//APPLY definitions
typedef Parameters<T, Basis1D, PDEOp, Preconditioner>                   Params;
typedef Parameters<T, Basis1D, WeightedPDEOp, WeightedPreconditioner>   WeightedParams;
typedef SYM_APPLY_1D<T,Index1D,Basis1D,Params,MA>                       APPLY;
typedef SYM_WEIGHTED_APPLY_1D<T,Basis1D,WeightedParams,WeightedMA>      WeightedAPPLY;


// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

/// Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff_it;

// Weight definition
const T eta=2;

T
weight(T x) {
    //return 2.+sin(x);

    if (fabs(x)>=1.) {
        return exp(-2*eta*fabs(x));
    }
    else {
        return exp(2*eta*(1./16.)*(-5-15*x*x+5*x*x*x*x-x*x*x*x*x*x));
    }

}

long double
test_weight(long double x) {
    //return 2.+sin(x);

    if (fabs(x)>=1.) {
        return exp(-2*eta*fabs(x));
    }
    else {
        return exp(2*eta*(1./16.)*(-5-15*x*x+5*x*x*x*x-x*x*x*x*x*x));
    }

}

T
prec_weight(T x) {
    return weight(x);
}

int main (int argc, char *argv[]) {
    if (argc != 9) {
        cout << "usage " << argv[0] << " d d_ j0 J set_J set_K_left set_K_right weighted" << endl;
        exit(1);
    }
    int d          =atoi(argv[1]);
    int d_         =atoi(argv[2]);
    int j0         =atoi(argv[3]);   //minimal level for basis
    int J          =atoi(argv[4]);   //For each column index, add J row levels by lamdabTilde routine
    int set_J      =atoi(argv[5]);   //maximum level difference for column index set
    int set_K_left =atoi(argv[6]);   //indcates left translation range for column index set
    int set_K_right=atoi(argv[7]);   //indcates right translation range for column index set
    int weighted   =atoi(argv[8]);   //use a weighted PDE-operator (=1) or not (else).

    Basis1D basis(d,d_,j0);
    Coefficients<Lexicographical,T,Index1D> u, Au;
    Coefficients<AbsoluteValue,T,Index1D> u_abs;

    IndexSet<Index1D> Lambda;
    for (int k=set_K_left; k<=set_K_right; ++k) {
        Lambda.insert(Index1D(j0,k,XBSpline));
    }
    for (int j=0; j<=set_J; ++j) {
        for (int k=pow2i<T>(j)*set_K_left; k<=pow2i<T>(j)*set_K_right; ++k) {
            Lambda.insert(Index1D(j0+j,k,XWavelet));
        }
    }


    cout << "Size of Lambda: " << Lambda.size() << endl;

    int count=1;
    for (const_set_it it=Lambda.begin(); it!=Lambda.end(); ++it,++count) {
        u[*it] = std::pow(T(count),-1.51);
    }

    u_abs = u;
    ofstream file("nterm.txt");
    for (int n=1; n<=int(Lambda.size()/2); ++n) {
        file << n << " " << u_abs.l2bestnterm(n) << endl;
    }
    file.close();

    if (weighted==0) {
        PDEOp           op(basis,1.);
        Preconditioner  P(basis);
        Compression     Compr(basis);
        MA              A(op,P,Compr);

        Params  parameters(basis,op,true,j0);
        APPLY   Apply(parameters,basis,A);

        count=1;
        for (const_set_it col=Lambda.begin(); col!=Lambda.end(); ++col,++count) {
            int s_tilde_level = J+1;
            int s_tilde_singsupp = J+1;
            cout << *col << endl;
            IndexSet<Index1D> LambdaTilde = lambdaTilde1d_WeightedPDE(*col, basis, s_tilde_level,
                                                                      j0, j0+J, s_tilde_singsupp);
            for (const_set_it row=LambdaTilde.begin(); row!=LambdaTilde.end(); ++row) {
                Au[(*row)] += P(*row)*op(*row,*col)*P(*col) * u[*col];
            }
        }

        ofstream apply_conv_file("apply_conv.txt");
        for (int s_tilde_level=0; s_tilde_level<=J; ++s_tilde_level) {
            cout << "Computing Apply for s_tilde_level=" << s_tilde_level << endl;
            Coefficients<Lexicographical,T,Index1D> wk = Apply(u,s_tilde_level);
            Coefficients<Lexicographical,T,Index1D> diff;
            diff = Au-wk;
            apply_conv_file << s_tilde_level << " " << diff.norm(2.) << endl;
            cout << "   " << s_tilde_level << " " << diff.norm(2.) << endl;
        }
        apply_conv_file.close();
    }
    else {
        DenseVectorT singPts(2);
        singPts = -1., 1.;
        Function<T> weightFct(weight,singPts);
        Function<T> prec_weightFct(prec_weight,singPts);
        WeightedPDEOp           op(basis,1.,weightFct,4);
        WeightedPreconditioner  P(basis,prec_weightFct,1);
        Compression             Compr(basis);
        WeightedMA              A(op,P,Compr);

        WeightedParams  parameters(basis,op,true,j0);
        WeightedAPPLY   Apply(parameters,basis,A);

        count=1;
        for (const_set_it col=Lambda.begin(); col!=Lambda.end(); ++col,++count) {
            int s_tilde_level = J+1;
            int s_tilde_singsupp = J+1;
            IndexSet<Index1D> LambdaTilde = lambdaTilde1d_WeightedPDE(*col, basis, s_tilde_level,
                                                                      j0, j0+J, s_tilde_singsupp);
            for (const_set_it row=LambdaTilde.begin(); row!=LambdaTilde.end(); ++row) {
                Au[(*row)] += P(*row)*op(*row,*col)*P(*col) * u[*col];
            }
        }

        ofstream apply_conv_file("apply_conv.txt");
        for (int s_tilde_level=0; s_tilde_level<J-j0; ++s_tilde_level) {
            cout << "Computing Apply for s_tilde_level=" << s_tilde_level << endl;
            Coefficients<Lexicographical,T,Index1D> wk = Apply(u,s_tilde_level);
            Coefficients<Lexicographical,T,Index1D> diff;
            diff = Au-wk;
            apply_conv_file << s_tilde_level << " " << diff.norm(2.) << endl;
            cout << "   " << s_tilde_level << " " << diff.norm(2.) << endl;
        }
        apply_conv_file.close();
    }



    return 0;
}

