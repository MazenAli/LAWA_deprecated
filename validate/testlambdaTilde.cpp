#include <iostream>
#include <lawa/lawa.h>

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
typedef HelmholtzOperator1D<T, Basis1D>         PDEOp;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

/// Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff_it;


template <typename T>
Range<int>
getRange(const Basis<T,Primal,R,CDF> &basis, XType e, int j);

template <typename T, Construction Cons>
Range<int>
getRange(const Basis<T,Primal,Interval,Cons> &basis, XType e, int j);

template <typename T, Construction Cons>
Range<int>
getRange(const Basis<T,Primal,Periodic,Cons> &basis, XType e, int j);

int main (int argc, char *argv[]) {
    if (argc != 5) {
        cout << "usage " << argv[0] << " d d_ j0 J" << endl; exit(1);
    }
    int d =atoi(argv[1]);
    int d_=atoi(argv[2]);
    int j0=atoi(argv[3]);
    int J =atoi(argv[4]);

    Basis1D basis(d,d,j0);
    PDEOp   op(basis,1.);

    T max_abs_error = 0.;

    //Test B-Spline column indices
    for (int k_col=getRange(basis,XBSpline,j0).firstIndex();
             k_col<=getRange(basis,XBSpline,j0).lastIndex(); ++k_col) {
        Index1D col_index(j0,k_col,XBSpline);
        Coefficients<Lexicographical,T,Index1D> nonzeros;
        for (int k_row=getRange(basis,XBSpline,j0).firstIndex();
                 k_row<=getRange(basis,XBSpline,j0).lastIndex(); ++k_row) {
            Index1D row_index(j0,k_row,XBSpline);
            T tmp = op(row_index,col_index);
            if (fabs(tmp)>thresh) {
                nonzeros[row_index] = tmp;
            }
        }
        for (int j_row=j0; j_row<=J; ++j_row) {
            for (int k_row=getRange(basis,XWavelet,j_row).firstIndex();
                k_row<=getRange(basis,XWavelet,j_row).lastIndex(); ++k_row) {
                Index1D row_index(j_row,k_row,XWavelet);
                T tmp = op(row_index,col_index);
                if (fabs(tmp)>thresh) {
                    nonzeros[row_index] = tmp;
                }
            }
        }
        cout << "Current index: " << col_index << endl;
        IndexSet<Index1D> lambdaTilde_nonzeros;
        lambdaTilde_nonzeros=lambdaTilde1d_PDE(col_index, basis, J, j0, J, false);

        for (const_coeff_it it=nonzeros.begin(); it!=nonzeros.end();++it) {
            if (lambdaTilde_nonzeros.count((*it).first)==0) {
                max_abs_error = std::max(fabs((*it).second),max_abs_error);
                //cout << "Error: " << col_index << " " << (*it).first
                //     << " (value=)" << (*it).second << " is missing." << endl;
            }
        }
        for (const_set_it it=lambdaTilde_nonzeros.begin(); it!=lambdaTilde_nonzeros.end();++it) {
            if ( (nonzeros.count((*it))==0) &&
                 ( (*it).k >= getRange(basis,(*it).xtype,(*it).j).firstIndex()) &&
                 ( (*it).k <= getRange(basis,(*it).xtype,(*it).j).lastIndex()) )
            {
                T tmp = op(*it,col_index);
                cout << "Inefficient: " << basis.generator(XBSpline).singularSupport(j0,k_col)
                     << " " << basis.generator((*it).xtype).singularSupport((*it).j,(*it).k)
                     << " is not required: " << tmp << endl;
            }
        }
    }
    //Test Wavelet column indices
    for (int j_col=j0; j_col<=J; ++j_col) {
        for (int k_col=getRange(basis,XWavelet,j_col).firstIndex();
                 k_col<=getRange(basis,XWavelet,j_col).lastIndex(); ++k_col) {
            Index1D col_index(j_col,k_col,XWavelet);
            Coefficients<Lexicographical,T,Index1D> nonzeros;
            for (int k_row=getRange(basis,XBSpline,j0).firstIndex();
                     k_row<=getRange(basis,XBSpline,j0).lastIndex(); ++k_row) {
                Index1D row_index(j0,k_row,XBSpline);
                T tmp = op(row_index,col_index);
                if (fabs(tmp)>thresh) {
                    nonzeros[row_index] = tmp;
                }
            }
            for (int j_row=j0; j_row<=J; ++j_row) {
                for (int k_row=getRange(basis,XWavelet,j_row).firstIndex();
                    k_row<=getRange(basis,XWavelet,j_row).lastIndex(); ++k_row) {
                    Index1D row_index(j_row,k_row,XWavelet);
                    T tmp = op(row_index,col_index);
                    if (fabs(tmp)>thresh) {
                        nonzeros[row_index] = tmp;
                    }
                }
            }
            cout << "Current index: " << col_index << endl;
            IndexSet<Index1D> lambdaTilde_nonzeros;
            lambdaTilde_nonzeros=lambdaTilde1d_PDE(col_index, basis, J, j0, J, false);

            for (const_coeff_it it=nonzeros.begin(); it!=nonzeros.end();++it) {
                if (lambdaTilde_nonzeros.count((*it).first)==0) {
                    max_abs_error = std::max(fabs((*it).second),max_abs_error);
                    //cout << "Error: " << col_index << " " << (*it).first
                    //     << " (value=)" << (*it).second << " is missing." << endl;
                }
            }
            for (const_set_it it=lambdaTilde_nonzeros.begin(); it!=lambdaTilde_nonzeros.end();++it) {
                if ( (nonzeros.count((*it))==0) &&
                     ( (*it).k >= getRange(basis,(*it).xtype,(*it).j).firstIndex()) &&
                     ( (*it).k <= getRange(basis,(*it).xtype,(*it).j).lastIndex()) )
                {
                    T tmp = op(*it,col_index);
                    cout << "Inefficient: " << basis.generator(XWavelet).singularSupport(j_col,k_col)
                         << " " << basis.generator((*it).xtype).singularSupport((*it).j,(*it).k)
                         << " is not required: " << tmp << endl;
                }
            }
        }
    }
    cout << "Maximum absolute error: " << max_abs_error << endl;
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
