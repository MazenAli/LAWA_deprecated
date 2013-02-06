/* TEST LOCAL OPERATOR
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

typedef double T;

///  Wavelet basis over an interval
typedef Basis<T, Primal, Interval, Dijkema>                       	IntervalBasis;
typedef Basis<T, Primal, Periodic, CDF>		                        PeriodicBasis;
typedef IntervalBasis::RefinementBasis                              IntervalRefinementBasis;
typedef PeriodicBasis::RefinementBasis                              PeriodicRefinementBasis;

typedef TensorBasis2D<Adaptive,IntervalBasis,IntervalBasis>         IntervalBasis2D;
typedef TensorBasis2D<Adaptive,PeriodicBasis,PeriodicBasis>         PeriodicBasis2D;
typedef TensorBasis2D<Adaptive,IntervalBasis,PeriodicBasis>         IntervalPeriodicBasis2D;
typedef TensorBasis2D<Adaptive,PeriodicBasis,IntervalBasis>         PeriodicIntervalBasis2D;

template<typename Basis>
void
getSparseGridIndexSet(const Basis &basis, IndexSet<Index2D> &Lambda, int j, T gamma=0.);

int main (int argc, char *argv[]) {

#ifdef TRONE
    cout << "using tr1." << endl;
#else
    cout << "using gnu_cxx." << endl;
#endif
    cout.precision(6);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);

    bool useSparseGrid=true;

    /// Basis initialization
    IntervalBasis intervalbasis(d, d, j0);
    intervalbasis.enforceBoundaryCondition<DirichletBC>();
    PeriodicBasis periodicbasis(d, d, j0);


    IntervalBasis2D intervalbasis2d(intervalbasis,intervalbasis);
    PeriodicBasis2D periodicbasis2d(periodicbasis,periodicbasis);
    IntervalPeriodicBasis2D intervalperiodicbasis2d(intervalbasis,periodicbasis);
    PeriodicIntervalBasis2D periodicintervalbasis2d(periodicbasis,intervalbasis);

	
	for(int j = 0; j <= J-j0; ++j){

		// Indizes
		Index1D index1_x_p(j0+j+1,1,XWavelet);
	    Index1D index1_y_p(j0+j+1,4,XWavelet);
	    Index1D index1_x_i(j0+j+1,2,XWavelet);
	    Index1D index1_y_i(j0+j+1,3,XWavelet);

	    Index1D index2_x_p(j0+j+3,17,XWavelet);
	    Index1D index2_y_p(j0+j+3,3,XWavelet);
	    Index1D index2_x_i(j0+j+3,9,XWavelet);
	    Index1D index2_y_i(j0+j+3,11,XWavelet);

	    // Construct initial sets by completing the indizes

		Coefficients<Lexicographical,T,Index2D> v_i, v_p, v_ip, v_pi,
												Cv_i, Cv_p, Cv_ip, Cv_pi;
		IndexSet<Index2D>			Cdiffv_i, Cdiffv_p, Cdiffv_ip, Cdiffv_pi;

		cout << "===== INTERVAL: j = " << j << " ======= " << endl << endl;
		Index2D index1_ii(index1_x_i, index1_y_i);
		completeMultiTree(intervalbasis2d,index1_ii,v_i,0,true);
		cout << "Initial Set v: " << v_i.size() << v_i << endl;

		extendMultiTree(intervalbasis2d,v_i,Cv_i, "standard",false,true);
		cout << "Extended Set Cv: " << Cv_i.size() << Cv_i << endl;
		Cv_i.clear();
		extendMultiTree(intervalbasis2d,v_i,Cv_i, Cdiffv_i, "standard", false,true);
		cout << "Extended Set Diff Cv: " << Cdiffv_i.size() << " (" << Cv_i.size() << ") " << Cdiffv_i << endl << endl;

		cout << "===== PERIODIC: j = " << j << " ======= " << endl << endl;
		Index2D index1_pp(index1_x_p, index1_y_p);
 		completeMultiTree(periodicbasis2d,index1_pp,v_p,0,true);
		cout << "Initial Set v: " << v_p.size() << v_p << endl;

		extendMultiTree(periodicbasis2d,v_p,Cv_p, "standard",false,true);
		cout << "Extended Set Cv: " << Cv_p.size() << Cv_p << endl;
		Cv_p.clear();
		extendMultiTree(periodicbasis2d,v_p,Cv_p, Cdiffv_p, "standard", false,true);
		cout << "Extended Set Diff Cv: " << Cdiffv_p.size() << " (" << Cv_p.size() << ") " << Cdiffv_p << endl << endl;

		cout << "===== INTERVAL x PERIODIC: j = " << j << " ======= " << endl << endl;
		Index2D index1_ip(index1_x_i, index1_y_p);
 		completeMultiTree(intervalperiodicbasis2d,index1_ip,v_ip,0,true);
		cout << "Initial Set v: " << v_ip.size() << v_ip << endl;

		extendMultiTree(intervalperiodicbasis2d,v_ip,Cv_ip, "standard",false,true);
		cout << "Extended Set Cv: " << Cv_ip.size() << Cv_ip << endl;
		Cv_ip.clear();
		extendMultiTree(intervalperiodicbasis2d,v_ip,Cv_ip, Cdiffv_ip, "standard", false,true);
		cout << "Extended Set Diff Cv: " << Cdiffv_ip.size() << " (" << Cv_ip.size() << ") " << Cdiffv_ip << endl << endl;

		cout << "===== PERIODIC X INTERVAL: j = " << j << " ======= " << endl << endl;
		Index2D index1_pi(index1_x_p, index1_y_i);
 		completeMultiTree(periodicintervalbasis2d,index1_pi,v_pi,0,true);
		cout << "Initial Set v: " << v_pi.size() << v_pi << endl;

		extendMultiTree(periodicintervalbasis2d,v_pi,Cv_pi, "standard",false,true);
		cout << "Extended Set Cv: " << Cv_pi.size() << Cv_pi << endl;
		Cv_pi.clear();
		extendMultiTree(periodicintervalbasis2d,v_pi,Cv_pi, Cdiffv_pi, "standard", false,true);
		cout << "Extended Set Diff Cv: " << Cdiffv_pi.size() << " (" << Cv_pi.size() << ") " << Cdiffv_pi << endl << endl;

	}


	return 0;
}


template<typename Basis>
void
getSparseGridIndexSet(const Basis &basis, IndexSet<Index2D> &Lambda, int j, T gamma)
{
    int j0 = basis.j0;
    for (long k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        for (long k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D row(j0,k1,XBSpline);
            Index1D col(j0,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0+i2-1;
            for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
                Lambda.insert(Index2D(col,row));
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0+i1-1;
        for (int i2=1; i1+i2<=j; ++i2) {
            if (T(i1+i2)-gamma*max(i1,i2)>(1-gamma)*j) continue;
            int j2=j0+i2-1;
            for (long k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    Lambda.insert(Index2D(row,col));
                }
            }
        }
    }
    return;
}
