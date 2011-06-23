#include <iostream>
#include <lawa/lawa.h>
#include <applications/unbounded_domains/parameters/parameters.h>

using namespace std;
using namespace lawa;

typedef double T;

/// Thresh bound
const T thresh = 1e-15;

/// Basis definitions
typedef Basis<T,Orthogonal,R,Multi>                                 Basis1D;

/// Operator definitions
typedef HelmholtzOperator1D<T, Basis1D>                             HelmholtzOp;
typedef H1NormPreconditioner1D<T,Basis1D>                           Preconditioner;
typedef NoCompression<T,Index1D,Basis1D>                            Compression;

/// Matrix definition
typedef MapMatrix<T,Index1D,HelmholtzOp,Compression,Preconditioner> MA;

//APPLY definitions
typedef Parameters<T, Basis1D, HelmholtzOp, Preconditioner>         Params;
typedef SYM_APPLY_1D<T,Index1D,Basis1D,Params,MA>                   APPLY;



// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

/// Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff_it;


int main (int argc, char *argv[]) {
    if (argc != 7) {
        cout << "usage " << argv[0] << " d j0 J set_J set_K_left set_K_right" << endl;
        exit(1);
    }
    int d          =atoi(argv[1]);
    int j0         =atoi(argv[2]);   //minimal level for basis
    int J          =atoi(argv[3]);   //For each column index, add J row levels by lamdabTilde routine
    int set_J      =atoi(argv[4]);   //maximum level difference for column index set
    int set_K_left =atoi(argv[5]);   //indcates left translation range for column index set
    int set_K_right=atoi(argv[6]);   //indcates right translation range for column index set

    Basis1D basis(d,j0);
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


    HelmholtzOp     helmholtz_op(basis,1.);
    Preconditioner  prec(basis);
    Compression     compr(basis);
    MA              A(helmholtz_op,prec,compr);

    Params  parameters(basis,helmholtz_op,j0);
    APPLY   Apply(parameters,basis,A);

    count=1;
    for (const_set_it col=Lambda.begin(); col!=Lambda.end(); ++col,++count) {
        int s_tilde_level = J-j0;
        IndexSet<Index1D> LambdaTilde = lambdaTilde1d_PDE(*col, basis, s_tilde_level, j0, J, false);
        for (const_set_it row=LambdaTilde.begin(); row!=LambdaTilde.end(); ++row) {
            Au[(*row)] += prec(*row)*helmholtz_op(*row,*col)*prec(*col) * u[*col];
        }
    }

    ofstream apply_conv_file("apply_conv.txt");
    for (int s_tilde_level=0; s_tilde_level<=2*J; ++s_tilde_level) {
        cout << "Computing Apply for s_tilde_level=" << s_tilde_level << endl;
        Coefficients<Lexicographical,T,Index1D> wk = Apply(u,s_tilde_level,J);
        Coefficients<Lexicographical,T,Index1D> diff;
        diff = Au-wk;
        apply_conv_file << s_tilde_level << " " << diff.norm(2.) << endl;
        cout << "   " << s_tilde_level << " " << diff.norm(2.) << endl;
    }
    apply_conv_file.close();




    return 0;
}

