#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/source/loc_single_scale_transforms.h>
#include <applications/new_eval_scheme/source/new_eval.h>
#include <applications/new_eval_scheme/source/localoperator.h>
#include <applications/new_eval_scheme/source/localoperator2d.h>
#include <lawa/methods/adaptive/datastructures/alignedindexset.h>
#include <lawa/methods/adaptive/datastructures/alignedcoefficients.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
//typedef HelmholtzOperator1D<T,PrimalBasis>                          BilinearForm;
typedef IdentityOperator1D<T,PrimalBasis>                          BilinearForm;
typedef DiagonalMatrixPreconditioner1D<T,PrimalBasis,BilinearForm>  Preconditioner;

typedef LocalOperator<PrimalBasis,PrimalBasis, BilinearForm, Preconditioner> LocalOp1D;
typedef LocalOperator2D<PrimalBasis, LocalOp1D, LocalOp1D>                   LocalOp2D;

typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;


typedef AlignedCoefficients<T,Index2D,Index1D,Index1D>              alignedCoefficients;

void
getSparseGridIndexSet(const PrimalBasis &basis, Coefficients<Lexicographical,T,Index2D> &coeff, int j);

void
refComputation1(const BilinearForm &Bil, const Coefficients<Lexicographical,T,Index2D> &v,
                Coefficients<Lexicographical,T,Index2D> &Mv);
void
refComputation2(const BilinearForm &Bil, const Coefficients<Lexicographical,T,Index2D> &UIv,
                Coefficients<Lexicographical,T,Index2D> &MMv);

int main (int argc, char *argv[]) {

#ifdef TRONE
    cout << "using tr1." << endl;
#else
    cout << "using gnu_cxx." << endl;
#endif



    cout.precision(20);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);
    int J  = atoi(argv[4]);

    PrimalBasis     basis(d,d_,j0);
    //BilinearForm    Bil(trial_basis,0.);
    BilinearForm    Bil(basis);
    Preconditioner  Prec(Bil);
    int offset=5;
    if (d==2 && d_==2) {
        offset=2;
    }
    LocalOp1D localoperator(basis, false, basis, false, offset, Bil, Prec);
    LocalOp2D localop2d(basis,localoperator,localoperator);

    Timer time;


    ofstream file2("comptimes_mv2d.dat");
    Coefficients<Lexicographical,T,Index2D> v((J+3)*basis.mra.cardI(j0+J+3));
    Coefficients<Lexicographical,T,Index2D> MMv((J+3)*basis.mra.cardI(j0+J+3));
    Coefficients<Lexicographical,T,Index2D> UIv((J+3)*basis.mra.cardI(j0+J+3));

    for (int j=1; j<=J; ++j) {

        getSparseGridIndexSet(basis,v,j);

        for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
            UIv[(*it).first] = 0.;
            MMv[(*it).first] = 0.;
        }

        cout << "New scheme started..." << endl;
        time.start();
        localop2d.evalAA(v,UIv, MMv);
        time.stop();
        T time_evalAA1 = time.elapsed();
        cout << "New scheme finished." << endl;

        /*
        cout << "Reference calculation started..." << endl;
        Coefficients<Lexicographical,T,Index2D> UIv_ref, MMv_ref;
        refComputation1(Bil, v, UIv_ref);
        refComputation2(Bil, UIv_ref, MMv_ref);
        cout << "Reference calculation finished." << endl;

        UIv -= UIv_ref;
        MMv -= MMv_ref;

        cout << "Difference in norm UIv: " << UIv.norm(2.) << endl;
        cout << "Difference in norm MMv: " << MMv.norm(2.) << endl;


        cout << v.size() << " " << time_evalAA1 << endl;
        file2 << v.size() << " " << time_evalAA1 << endl;
        */


    }
    file2.close();

    return 0;
}

void
getSparseGridIndexSet(const PrimalBasis &basis, Coefficients<Lexicographical,T,Index2D> &coeff, int j)
{
    int j0 = basis.j0;
    for (long k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        for (long k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D row(j0,k1,XBSpline);
            Index1D col(j0,k2,XBSpline);
            coeff[Index2D(row,col)] = T(rand()) / T(RAND_MAX);
        }
        for (int i2=0; i2<=j; ++i2) {
            for (long k2=basis.rangeJ(j0+i2).firstIndex(); k2<=basis.rangeJ(j0+i2).lastIndex(); ++k2) {
                Index1D row(j0,k1,XBSpline);
                Index1D col(j0+i2,k2,XWavelet);
                coeff[Index2D(row,col)] = T(rand()) / T(RAND_MAX);
                coeff[Index2D(col,row)] = T(rand()) / T(RAND_MAX);
            }
        }
    }
    for (int i1=0; i1<=j; ++i1) {
        for (long k1=basis.rangeJ(j0+i1).firstIndex(); k1<=basis.rangeJ(j0+i1).lastIndex(); ++k1) {
            for (int i2=0; i1+i2<=j; ++i2) {
                for (long k2=basis.rangeJ(j0+i2).firstIndex(); k2<=basis.rangeJ(j0+i2).lastIndex(); ++k2) {
                    Index1D row(j0+i1,k1,XWavelet);
                    Index1D col(j0+i2,k2,XWavelet);
                    coeff[Index2D(row,col)] = T(rand()) / T(RAND_MAX);
                }
            }
        }
    }
}

void
refComputation1(const BilinearForm &Bil, const Coefficients<Lexicographical,T,Index2D> &v,
                Coefficients<Lexicographical,T,Index2D> &Mv)
{
    IndexSet<Index2D> Lambda = supp(v);
    cout << "   Size of Lambda = " << Lambda.size() << endl;
    for (const_set2d_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).index1;
        Index1D row_y = (*row).index2;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                if ( (row_x.xtype==XBSpline) || ((row_x.xtype==XWavelet && col_x.xtype==XWavelet  && row_x.j<=col_x.j)) ) {
                    val +=   Bil(row_x,col_x) * (*col).second;
                }
            }
        }
        Mv[*row] = val;
    }
}

void
refComputation2(const BilinearForm &Bil, const Coefficients<Lexicographical,T,Index2D> &UIv,
                Coefficients<Lexicographical,T,Index2D> &MMv)
{
    IndexSet<Index2D> Lambda = supp(UIv);
    cout << "   Size of Lambda = " << Lambda.size() << endl;
    for (const_set2d_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).index1;
        Index1D row_y = (*row).index2;
        for (const_coeff2d_it col=UIv.begin(); col!=UIv.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                val +=   Bil(row_y,col_y) * (*col).second;
            }
        }
        MMv[*row] = val;
    }
}
