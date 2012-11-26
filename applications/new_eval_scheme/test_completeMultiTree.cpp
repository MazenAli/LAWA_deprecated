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
typedef TensorBasis2D<Adaptive,PeriodicBasis,IntervalBasis>         PeriodicIntervalBasis2D;
typedef TensorBasis2D<Adaptive,IntervalBasis,PeriodicBasis>         IntervalPeriodicBasis2D;


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
    IntervalRefinementBasis &intervalrefinementbasis = intervalbasis.refinementbasis;
    IntervalBasis2D intervalbasis2d(intervalbasis,intervalbasis);
    PeriodicBasis periodicbasis(d, d, j0);
    PeriodicRefinementBasis &periodicrefinementbasis = periodicbasis.refinementbasis;
    PeriodicBasis2D periodicbasis2d(periodicbasis,periodicbasis);

    PeriodicIntervalBasis2D periodicintervalbasis2d(periodicbasis, intervalbasis);
    IntervalPeriodicBasis2D intervalperiodicbasis2d(intervalbasis, periodicbasis);


    cout << "#####################################################################" << endl;
    cout << "########                   1D                             ###########" << endl;
    cout << "#####################################################################" << endl;

	for(int j = 0; j <= J-j0; ++j){
		IndexSet<Index1D> Lambda_Interval, Lambda_Periodic;
		Coefficients<Lexicographical,T,Index1D> Coeffs_Interval, Coeffs_Periodic;


		cout << "===== INTERVAL: j = " << j << " ======= " << endl;

		Index1D index1_i(j0+j+2,16,XWavelet);
		cout << "Adding Index " << index1_i << " with support "
			 << intervalbasis.generator(index1_i.xtype).support(index1_i.j, index1_i.k) << endl;
		completeMultiTree(intervalbasis,index1_i,Coeffs_Interval,true);
		cout << "Extended Coeffs: " << Coeffs_Interval.size() << Coeffs_Interval << endl;

		cout << "===== PERIODIC: j = " << j << " ======= " << endl;

		Index1D index1_p(j0+j+2,15,XWavelet);
		cout << "Adding Index " << index1_p << " with support "
			 << periodicbasis.generator(index1_p.xtype).support(index1_p.j, index1_p.k) << endl;
		completeMultiTree(periodicbasis,index1_p,Coeffs_Periodic,true);
		cout << "Extended Coeffs: " << Coeffs_Periodic.size() << Coeffs_Periodic << endl;

		cout << endl << endl;
	}


	cout << endl << endl;

    cout << "#####################################################################" << endl;
    cout << "########                   2D                             ###########" << endl;
    cout << "#####################################################################" << endl;

	for(int j = 0; j <= J-j0; ++j){
		IndexSet<Index2D> Lambda_Interval, Lambda_Periodic, Lambda_PeriodicInterval, Lambda_IntervalPeriodic;
		Coefficients<Lexicographical,T,Index2D> Coeffs_Interval, Coeffs_Periodic, Coeffs_PeriodicInterval, Coeffs_IntervalPeriodic;
		
		Index1D index1_x_i(j0+j+2,16,XWavelet);
	    Index1D index1_y_i(j0+j+1,2,XWavelet);
		Index1D index1_x_p(j0+j+2,15,XWavelet);
	    Index1D index1_y_p(j0+j+1,1,XWavelet);

		cout << "===== INTERVAL x INTERVAL : j = " << j << " ======= " << endl;
		//getSparseGridIndexSet(intervalbasis,Lambda_Interval,j,0.2);
		//FillWithZeros(Lambda_Interval, Coeffs_Interval);
		
		cout << "Sparse Grid 2D: " << Lambda_Interval.size() << endl;// << Lambda_Interval << endl;		
	    Index2D new_index1_i(index1_x_i,index1_y_i);
		cout << "Adding Index " << new_index1_i << " with support " 
			 << intervalbasis2d.first.generator(index1_x_i.xtype).support(index1_x_i.j, index1_x_i.k) 
			 << " x " << intervalbasis2d.second.generator(index1_y_i.xtype).support(index1_y_i.j, index1_y_i.k)<< endl;
		completeMultiTree(intervalbasis2d,new_index1_i,Coeffs_Interval,0,true);
		cout << "Extended Coeffs: " << Coeffs_Interval.size() << endl << Coeffs_Interval << endl;
		
		cout << "===== PERIODIC x PERIODIC : j = " << j << " ======= " << endl;
		//getSparseGridIndexSet(periodicbasis,Lambda_Periodic,j,0.2);
		//FillWithZeros(Lambda_Periodic, Coeffs_Periodic);
		
		cout << "Sparse Grid 2D: " << Lambda_Periodic.size() << endl;// << Lambda_Periodic << endl;		
	    Index2D new_index1_p(index1_x_p,index1_y_p);
		cout << "Adding Index " << new_index1_p << " with support " 
			 << periodicbasis2d.first.generator(index1_x_p.xtype).support(index1_x_p.j, index1_x_p.k) 
			 << " x " << periodicbasis2d.second.generator(index1_y_p.xtype).support(index1_y_p.j, index1_y_p.k) << endl;
		completeMultiTree(periodicbasis2d,new_index1_p,Coeffs_Periodic,0,true);
		cout << "Extended Coeffs: " << Coeffs_Periodic.size() << endl << Coeffs_Periodic << endl;
		
		cout << "===== PERIODIC x  INTERVAL : j = " << j << " ======= " << endl;
		//getSparseGridIndexSet(periodicbasis,Lambda_Periodic,j,0.2);
		//FillWithZeros(Lambda_Periodic, Coeffs_Periodic);

		cout << "Sparse Grid 2D: " << Lambda_PeriodicInterval.size() << endl;// << Lambda_Periodic << endl;
	    Index2D new_index1_pi(index1_x_p,index1_y_i);
		cout << "Adding Index " << new_index1_pi << " with support "
			 << periodicintervalbasis2d.first.generator(index1_x_p.xtype).support(index1_x_p.j, index1_x_p.k)
			 << " x " << periodicintervalbasis2d.second.generator(index1_y_p.xtype).support(index1_y_p.j, index1_y_p.k) << endl;
		completeMultiTree(periodicintervalbasis2d,new_index1_pi,Coeffs_PeriodicInterval,0,true);
		cout << "Extended Coeffs: " << Coeffs_PeriodicInterval.size() << endl << Coeffs_PeriodicInterval << endl;

		cout << "===== INTERVAL x PERIODIC : j = " << j << " ======= " << endl;
		//getSparseGridIndexSet(periodicbasis,Lambda_Periodic,j,0.2);
		//FillWithZeros(Lambda_Periodic, Coeffs_Periodic);

		cout << "Sparse Grid 2D: " << Lambda_IntervalPeriodic.size() << endl;// << Lambda_Periodic << endl;
	    Index2D new_index1_ip(index1_x_i,index1_y_p);
		cout << "Adding Index " << new_index1_ip << " with support "
			 << intervalperiodicbasis2d.first.generator(index1_x_p.xtype).support(index1_x_p.j, index1_x_p.k)
			 << " x " << intervalperiodicbasis2d.second.generator(index1_y_p.xtype).support(index1_y_p.j, index1_y_p.k) << endl;
		completeMultiTree(intervalperiodicbasis2d,new_index1_ip,Coeffs_IntervalPeriodic,0,true);
		cout << "Extended Coeffs: " << Coeffs_IntervalPeriodic.size() << endl << Coeffs_IntervalPeriodic << endl;

		cout << endl << endl;
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
