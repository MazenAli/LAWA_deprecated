#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

/// Thresh bound
const T thresh = 1e-30;

// Basis definitions
const FunctionSide functionside = Primal;
const DomainType   domain = R;
const Construction construction = SparseMulti;

const T leftbound=-10.;   //for realline constructions, we only consider wavelets with support intersecting [left,right]
const T rightbound=10.;

typedef Basis<T,functionside,domain,construction>                   Basis1D;

/// Operator definitions
typedef HelmholtzOperator1D<T, Basis1D>                             HelmholtzOp;
typedef IdentityOperator1D<T, Basis1D>                              IdentityOp;

/// Preconditioner definitions
typedef H1NormPreconditioner1D<T, Basis1D>                          Preconditioner;
typedef NoPreconditioner<T, Index1D>                                NoPreconditioner1D;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

/// Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff_it;

template <typename T>
Range<int>
getRange(const Basis<T,Primal,R,CDF> &basis, XType e, int j);

template <typename T>
Range<int>
getRange(const Basis<T,Orthogonal,R,Multi> &basis, XType e, int j);

template <typename T>
Range<int>
getRange(const Basis<T,Primal,R,SparseMulti> &basis, XType e, int j);

template <typename T, Construction Cons>
Range<int>
getRange(const Basis<T,Primal,Interval,Cons> &basis, XType e, int j);

template <typename T, Construction Cons>
Range<int>
getRange(const Basis<T,Primal,Periodic,Cons> &basis, XType e, int j);

int main (int argc, char *argv[]) {
    cout.precision(16);
    if (argc != 5) {
        cout << "usage " << argv[0] << " d d_ j0 J" << endl; exit(1);
    }
    int d =atoi(argv[1]);
    int d_=atoi(argv[2]);
    int j0=atoi(argv[3]);
    int J =atoi(argv[4]);


    //Basis1D basis(d,d_,j0);
    Basis1D         basis(d,j0);
    HelmholtzOp     op(basis,1.);
    //IdentityOp       op(basis);
    Preconditioner  p(basis);
    //NoPreconditioner1D  p;

    T max_abs_error = 0.;
    Index1D max_error_col_index(0,0,XBSpline);
    Index1D max_error_row_index(0,0,XBSpline);


    /*
     * For a certain range of B-spline indices, we set up a coefficient vector where we calculate
     * all entries a(phi_{\lambda},phi_{mu}) and a(phi_{\lambda},psi_{mu}) for indices mu up
     * to a certain maximal level and a certain translation.
     * Afterwards, we calculate the sparsity pattern using lambdaTildePDE up to the maximal level.
     * Then, we perform two tests:
     * -> First, we check the difference between the calculated coefficient vector and the sparsity
     *    pattern, i.e., we check if the sparsity pattern is reliable.
     * -> Second, we check which entries appear in the coefficient vector which do not appear
     *    in the sparsity patter, i.e., we check if the sparsity pattern is efficient.
     */

    //Test B-Spline column indices:
    for (int k_col=getRange(basis,XBSpline,j0).firstIndex();
             k_col<=getRange(basis,XBSpline,j0).lastIndex(); ++k_col) {
        Index1D col_index(j0,k_col,XBSpline);
        Coefficients<Lexicographical,T,Index1D> nonzeros;
        for (int k_row=getRange(basis,XBSpline,j0).firstIndex();
                 k_row<=getRange(basis,XBSpline,j0).lastIndex(); ++k_row) {
            Index1D row_index(j0,k_row,XBSpline);
            T tmp = p(row_index)*op(row_index,col_index)*p(col_index);
            if (fabs(tmp)>thresh) {
                nonzeros[row_index] = tmp;
            }
        }
        for (int j_row=j0; j_row<=j0+J; ++j_row) {
            for (int k_row=getRange(basis,XWavelet,j_row).firstIndex();
                k_row<=getRange(basis,XWavelet,j_row).lastIndex(); ++k_row) {
                Index1D row_index(j_row,k_row,XWavelet);
                T tmp = p(row_index)*op(row_index,col_index)*p(col_index);
                if (fabs(tmp)>thresh) {
                    nonzeros[row_index] = tmp;
                }
            }
        }

        cout << "Current index: " << col_index << endl;
        IndexSet<Index1D> lambdaTilde_nonzeros;
        lambdaTilde_nonzeros=lambdaTilde1d_PDE(col_index, basis, J, j0, j0+J, false);
        //cout << "(XBSpline, " << j0 << ", " << k_col << "): [" << lambdaTilde_nonzeros.size() << "] "
        //     << lambdaTilde_nonzeros << endl;
        //getchar();


        for (const_coeff_it it=nonzeros.begin(); it!=nonzeros.end();++it) {
            if (lambdaTilde_nonzeros.count((*it).first)==0) {
                if (fabs((*it).second)> max_abs_error) {
                    max_abs_error = fabs((*it).second);
                    max_error_col_index = col_index;
                    max_error_row_index = (*it).first;
                }
                /*
                cout << "Error: " << col_index << ": " << (*it).first << ", "
                     << basis.generator(XBSpline).singularSupport(j0,k_col)
                     << " " << basis.generator((*it).first.xtype).singularSupport((*it).first.j,(*it).first.k)
                     << " is missing: ";
                printf("%E",fabs((*it).second));
                cout << endl << endl;
                */
            }
        }
        for (const_set_it it=lambdaTilde_nonzeros.begin(); it!=lambdaTilde_nonzeros.end();++it) {
            if ( (nonzeros.count((*it))==0) &&
                 ( (*it).k >= getRange(basis,(*it).xtype,(*it).j).firstIndex()) &&
                 ( (*it).k <= getRange(basis,(*it).xtype,(*it).j).lastIndex()) )
            {
                T tmp = p(*it)*op(*it,col_index)*p(col_index);
                /*
                cout << "Inefficient: "<< col_index << ": "  << *it << ", "
                     << basis.generator(XBSpline).singularSupport(j0,k_col)
                     << " " << basis.generator((*it).xtype).singularSupport((*it).j,(*it).k)
                     << " is not required: " << tmp << endl << endl;
                */
            }
        }
    }

    //Test Wavelet column indices
    /*
     * For a certain range of wavelets indices, we set up a coefficient vector where we calculate
     * all entries a(psi_{\lambda},phi_{mu}) and a(phi_{\lambda},psi_{mu}) for indices mu up
     * to a certain maximal level and a certain translation.
     * Afterwards, we calculate the sparsity pattern using lambdaTildePDE up to the maximal level.
     * Then, we perform two tests:
     * -> First, we check the difference between the calculated coefficient vector and the sparsity
     *    pattern, i.e., we check if the sparsity pattern is reliable.
     * -> Second, we check which entries appear in the coefficient vector which do not appear
     *    in the sparsity patter, i.e., we check if the sparsity pattern is efficient.
     */
    for (int j_col=j0; j_col<=j0+J; ++j_col) {
        for (int k_col=getRange(basis,XWavelet,j_col).firstIndex();
                 k_col<=getRange(basis,XWavelet,j_col).lastIndex(); ++k_col) {
            Index1D col_index(j_col,k_col,XWavelet);
            Coefficients<Lexicographical,T,Index1D> nonzeros;
            for (int k_row=getRange(basis,XBSpline,j0).firstIndex();
                     k_row<=getRange(basis,XBSpline,j0).lastIndex(); ++k_row) {
                Index1D row_index(j0,k_row,XBSpline);
                T tmp = p(row_index)*op(row_index,col_index)*p(col_index);
                if (fabs(tmp)>thresh) {
                    nonzeros[row_index] = tmp;
                }
            }
            for (int j_row=j0; j_row<=j0+J; ++j_row) {
                for (int k_row=getRange(basis,XWavelet,j_row).firstIndex();
                    k_row<=getRange(basis,XWavelet,j_row).lastIndex(); ++k_row) {
                    Index1D row_index(j_row,k_row,XWavelet);
                    T tmp = p(row_index)*op(row_index,col_index)*p(col_index);
                    if (fabs(tmp)>thresh) {
                        nonzeros[row_index] = tmp;
                    }
                }
            }

            cout << "Current index: " << col_index << endl;
            IndexSet<Index1D> lambdaTilde_nonzeros;
            lambdaTilde_nonzeros=lambdaTilde1d_PDE(col_index, basis, J, j0, j0+J, false);
            /*
            cout << "[" << lambdaTilde_nonzeros.size() << "] " << endl;
            for (const_set_it it=lambdaTilde_nonzeros.begin(); it!=lambdaTilde_nonzeros.end(); ++it) {
                cout << basis.psi.support(j_col,k_col) << " -> "
                     <<  basis.psi.support((*it).j,(*it).k) << endl;
            }
            getchar();
            */

            for (const_coeff_it it=nonzeros.begin(); it!=nonzeros.end();++it) {
                if (lambdaTilde_nonzeros.count((*it).first)==0) {
                    if (fabs((*it).second)> max_abs_error) {
                        max_abs_error = fabs((*it).second);
                        max_error_col_index = col_index;
                        max_error_row_index = (*it).first;
                    }
                    /*
                    cout << "Error: " << col_index << ": " << (*it).first << ", "
                         << basis.generator(XWavelet).singularSupport(j_col,k_col)
                         << " " << basis.generator((*it).first.xtype).singularSupport((*it).first.j,(*it).first.k)
                         << " is missing: ";
                    printf("%E",fabs((*it).second));
                    cout << endl << endl;
                    //getchar();
                    */
                }
            }
            for (const_set_it it=lambdaTilde_nonzeros.begin(); it!=lambdaTilde_nonzeros.end();++it) {
                if ( (nonzeros.count((*it))==0) &&
                     ( (*it).k >= getRange(basis,(*it).xtype,(*it).j).firstIndex()) &&
                     ( (*it).k <= getRange(basis,(*it).xtype,(*it).j).lastIndex()) )
                {
                    T tmp = p(*it)*op(*it,col_index)*p(col_index);
                    /*
                    cout << "Inefficient: "<< col_index << ": "  << *it << ", "
                         << basis.generator(XWavelet).singularSupport(j_col,k_col)
                         << " " << basis.generator((*it).xtype).singularSupport((*it).j,(*it).k)
                         << " is not required: " << tmp << endl << endl;
                    */

                }
            }
        }
    }
    cout << "Maximum absolute error: " ;
    printf("%E", max_abs_error);
    cout << " observed for: " << endl;
    cout << max_error_col_index << ": " << basis.generator(max_error_col_index.xtype).singularSupport(max_error_col_index.j,max_error_col_index.k) << endl;
    cout << max_error_row_index << ": " << basis.generator(max_error_row_index.xtype).singularSupport(max_error_row_index.j,max_error_row_index.k) << endl;

    return 0;
}

template <typename T>
Range<int>
getRange(const Basis<T,Primal,R,CDF> &basis, XType e, int j)
{
    T l1 = basis.generator(e).support(0,0).l1;
    T l2 = basis.generator(e).support(0,0).l2;
    return _(int(pow2i<T>(j)*leftbound-l2),int(pow2i<T>(j)*rightbound-l1));
}

template <typename T>
Range<int>
getRange(const Basis<T,Orthogonal,R,Multi> &basis, XType e, int j)
{
    const BSpline<T,Orthogonal,R,Multi> &phi = basis.mra.phi;
    const Wavelet<T,Orthogonal,R,Multi> &psi = basis.psi;

    T l1=0., l2=1.;
    int numFunc=1;
    if (e==XBSpline) {
        Support<T> supp = phi.max_support();
        l1 = supp.l1;
        l2 = supp.l2;
        numFunc = phi._numSplines;
    }
    else {
        Support<T> supp = psi.max_support();
        l1 = supp.l1;
        l2 = supp.l2;
        numFunc = psi._numSplines;
    }

    return _(int(pow2i<T>(j)*leftbound-numFunc),int(pow2i<T>(j)*rightbound+numFunc));
}

template <typename T>
Range<int>
getRange(const Basis<T,Primal,R,SparseMulti> &basis, XType e, int j)
{
    const BSpline<T,Primal,R,SparseMulti> &phi = basis.mra.phi;
    const Wavelet<T,Primal,R,SparseMulti> &psi = basis.psi;

    T l1=0., l2=1.;
    int numFunc=1;
    if (e==XBSpline) {
        Support<T> supp = phi.max_support();
        l1 = supp.l1;
        l2 = supp.l2;
        numFunc = phi._numSplines;
    }
    else {
        Support<T> supp = psi.max_support();
        l1 = supp.l1;
        l2 = supp.l2;
        numFunc = psi._numSplines;
    }

    return _(int(pow2i<T>(j)*leftbound-numFunc),int(pow2i<T>(j)*rightbound+numFunc));
}

template <typename T, Construction Cons>
Range<int>
getRange(const Basis<T,Primal,Interval,Cons> &basis, XType e, int j)
{
    if (e==XBSpline)    return basis.mra.rangeI(j);
    else                return basis.rangeJ(j);
}

template <typename T, Construction Cons>
Range<int>
getRange(const Basis<T,Primal,Periodic,Cons> &basis, XType e, int j)
{
    if (e==XBSpline)    return basis.mra.rangeI(j);
    else                return basis.rangeJ(j);
}
