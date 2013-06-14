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
//typedef Basis<T, Primal, Interval, Dijkema>                       	IntervalBasis;
typedef Basis<T, Orthogonal, Interval, Multi>                       	IntervalBasis;
typedef Basis<T, Primal, Periodic, CDF>		                        PeriodicBasis;
typedef IntervalBasis::RefinementBasis                              IntervalRefinementBasis;
typedef PeriodicBasis::RefinementBasis                              PeriodicRefinementBasis;

typedef TensorBasis2D<Adaptive,IntervalBasis,IntervalBasis>         IntervalBasis2D;
typedef TensorBasis2D<Adaptive,PeriodicBasis,PeriodicBasis>         PeriodicBasis2D;
typedef TensorBasis2D<Adaptive,PeriodicBasis,IntervalBasis>         PeriodicIntervalBasis2D;
typedef TensorBasis2D<Adaptive,IntervalBasis,PeriodicBasis>         IntervalPeriodicBasis2D;

typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XOne>             XOneAlignedCoefficients;
typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XTwo>             XTwoAlignedCoefficients;

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
    //IntervalBasis intervalbasis(d, d, j0);
    IntervalBasis intervalbasis(d,0);
    intervalbasis.enforceBoundaryCondition<DirichletBC>();
    IntervalRefinementBasis &intervalrefinementbasis = intervalbasis.refinementbasis;
    IntervalBasis2D intervalbasis2d(intervalbasis,intervalbasis);
    PeriodicBasis periodicbasis(d, d, j0);
    PeriodicRefinementBasis &periodicrefinementbasis = periodicbasis.refinementbasis;
    PeriodicBasis2D periodicbasis2d(periodicbasis,periodicbasis);

    PeriodicIntervalBasis2D periodicintervalbasis2d(periodicbasis, intervalbasis);
    IntervalPeriodicBasis2D intervalperiodicbasis2d(intervalbasis, periodicbasis);


    cout << "#####################################################################" << endl;
    cout << "########                   1D                             ###########" << endl;
    cout << "#####################################################################" << endl << endl;

	for(int j = 0; j <= J-j0; ++j){
		IndexSet<Index1D> Lambda_Interval, Lambda_Periodic;
		Coefficients<Lexicographical,T,Index1D> Coeffs_Interval, Coeffs_Periodic;


		cout << "===== INTERVAL: j = " << j << " ======= " << endl;

		Index1D index1_i(j0+j+2,16,XWavelet);
		cout << "Adding Index " << index1_i << " with support "
			 << intervalbasis.generator(index1_i.xtype).support(index1_i.j, index1_i.k) << endl;
		completeMultiTree(intervalbasis,index1_i,Coeffs_Interval,true);
		cout << "Extended Coeffs: " << Coeffs_Interval.size() << Coeffs_Interval << endl << endl;

		cout << "--------    2nd Completion -- return added indizes   ---------------" << endl << endl;

		IndexSet<Index1D> additional_indizes_i;
		Index1D index2_i(j0+j+3,16,XWavelet);
		cout << "Adding Index " << index1_i << " with support "
			 << intervalbasis.generator(index2_i.xtype).support(index2_i.j, index2_i.k) << endl;
		completeMultiTree(intervalbasis,index2_i,Coeffs_Interval,additional_indizes_i,true);
		cout << "Extended Coeffs: " << Coeffs_Interval.size() << Coeffs_Interval << endl;
		cout << "Added Indizes: " << additional_indizes_i.size() << additional_indizes_i << endl << endl;

		cout << "===== PERIODIC: j = " << j << " ======= " << endl;

		Index1D index1_p(j0+j+2,15,XWavelet);
		cout << "Adding Index " << index1_p << " with support "
			 << periodicbasis.generator(index1_p.xtype).support(index1_p.j, index1_p.k) << endl;
		completeMultiTree(periodicbasis,index1_p,Coeffs_Periodic,true);
		cout << "Extended Coeffs: " << Coeffs_Periodic.size() << Coeffs_Periodic << endl << endl;

		cout << "--------    2nd Completion -- return added indizes   ---------------" << endl << endl;

		IndexSet<Index1D> additional_indizes_p;
		Index1D index2_p(j0+j+3,15,XWavelet);
		cout << "Adding Index " << index2_p << " with support "
			 << periodicbasis.generator(index2_p.xtype).support(index2_p.j, index2_p.k) << endl;
		completeMultiTree(periodicbasis,index2_p,Coeffs_Periodic,additional_indizes_p,true);
		cout << "Extended Coeffs: " << Coeffs_Periodic.size() << Coeffs_Periodic << endl;
		cout << "Added Indizes: " << additional_indizes_p.size() << additional_indizes_p << endl;

		cout << endl << endl;


	}


	cout << endl << endl;

    cout << "#####################################################################" << endl;
    cout << "########                   2D                             ###########" << endl;
    cout << "#####################################################################" << endl;

	for(int j = 0; j <= J-j0; ++j){
		IndexSet<Index2D> Lambda_Interval, Lambda_Periodic, Lambda_PeriodicInterval, Lambda_IntervalPeriodic;
		Coefficients<Lexicographical,T,Index2D> Coeffs_Interval, Coeffs_Periodic, Coeffs_PeriodicInterval, Coeffs_IntervalPeriodic;
		
		/*Index1D index1_x_i(j0+j+2,16,XWavelet);
	    Index1D index1_y_i(j0+j+1,2,XWavelet);
		Index1D index1_x_p(j0+j+2,15,XWavelet);
	    Index1D index1_y_p(j0+j+1,1,XWavelet);
	    */

        Index1D index1_x_p(j0+j+4,42,XWavelet);
        Index1D index1_y_p(j0+j+4,4,XWavelet);
        Index1D index1_x_i(j0+j+4,7,XWavelet);
        Index1D index1_y_i(j0+j+4,17,XWavelet);

        Index1D index2_x_p(j0+j+5,30,XWavelet);
        Index1D index2_y_p(j0+j+5,3,XWavelet);
        Index1D index2_x_i(j0+j+5,9,XWavelet);
        Index1D index2_y_i(j0+j+5,11,XWavelet);

		cout << "===== INTERVAL x INTERVAL : j = " << j << " ======= " << endl;
		getSparseGridIndexSet(intervalbasis2d,Lambda_Interval,j,0.2);
		FillWithZeros(Lambda_Interval, Coeffs_Interval);
		
		cout << "Sparse Grid 2D: " << Lambda_Interval.size() << endl;// << Lambda_Interval << endl;		
	    Index2D new_index1_i(index1_x_i,index1_y_i);
		cout << "Adding Index " << new_index1_i << " with support " 
			 << intervalbasis2d.first.generator(index1_x_i.xtype).support(index1_x_i.j, index1_x_i.k) 
			 << " x " << intervalbasis2d.second.generator(index1_y_i.xtype).support(index1_y_i.j, index1_y_i.k)<< endl;
		completeMultiTree(intervalbasis2d,new_index1_i,Coeffs_Interval,0,true);
		cout << "Extended Coeffs: " << Coeffs_Interval.size() << endl << Coeffs_Interval << endl;
		
		/*cout << "X1 Alignment " << endl;
	    XOneAlignedCoefficients x1aligned(6151,193);
	    x1aligned.align(Coeffs_Interval,J+j0);
	    for (XOneAlignedCoefficients::const_map_prindex_it it=x1aligned.map.begin();
	                                                            it!=x1aligned.map.end(); ++it) {
	    	cout << (*it).first << (*it).second << endl;
	    }

		cout << "X2 Alignment " << endl;
	    XTwoAlignedCoefficients x2aligned(6151,193);
	    x2aligned.align(Coeffs_Interval,J+j0);
	    for (XTwoAlignedCoefficients::const_map_prindex_it it=x2aligned.map.begin();
	                                                            it!=x2aligned.map.end(); ++it) {
	    	cout << (*it).first << (*it).second << endl;
	    }*/

		cout << "--------    2nd Completion -- return added indizes   ---------------" << endl << endl;

		IndexSet<Index2D> additional_indizes_ii;
	    Index2D new_index2_i(index2_x_i,index2_y_i);
		cout << "Adding Index " << new_index2_i << " with support "
			 << intervalbasis2d.first.generator(index2_x_i.xtype).support(index2_x_i.j, index2_x_i.k)
			 << " x " << intervalbasis2d.second.generator(index2_y_i.xtype).support(index2_y_i.j, index2_y_i.k)<< endl;
		completeMultiTree(intervalbasis2d,new_index2_i,Coeffs_Interval,additional_indizes_ii,0,true);
		cout << "Extended Coeffs: " << Coeffs_Interval.size() << endl << Coeffs_Interval << endl;
		cout << "Added Indizes: " << additional_indizes_ii.size() << additional_indizes_ii << endl << endl;

		cout << "===== PERIODIC x PERIODIC : j = " << j << " ======= " << endl;
		getSparseGridIndexSet(periodicbasis2d,Lambda_Periodic,j,0.2);
		FillWithZeros(Lambda_Periodic, Coeffs_Periodic);
		
		cout << "Sparse Grid 2D: " << Lambda_Periodic.size() << endl;// << Lambda_Periodic << endl;		
	    Index2D new_index1_p(index1_x_p,index1_y_p);
		cout << "Adding Index " << new_index1_p << " with support " 
			 << periodicbasis2d.first.generator(index1_x_p.xtype).support(index1_x_p.j, index1_x_p.k) 
			 << " x " << periodicbasis2d.second.generator(index1_y_p.xtype).support(index1_y_p.j, index1_y_p.k) << endl;
		completeMultiTree(periodicbasis2d,new_index1_p,Coeffs_Periodic,0,true);
		cout << "Extended Coeffs: " << Coeffs_Periodic.size() << endl << Coeffs_Periodic << endl;
		
		cout << "--------    2nd Completion -- return added indizes   ---------------" << endl << endl;

		IndexSet<Index2D> additional_indizes_pp;
	    Index2D new_index2_p(index2_x_p,index2_y_p);
		cout << "Adding Index " << new_index2_p << " with support "
			 << periodicbasis2d.first.generator(index2_x_p.xtype).support(index2_x_p.j, index2_x_p.k)
			 << " x " << periodicbasis2d.second.generator(index2_y_p.xtype).support(index2_y_p.j, index2_y_p.k) << endl;
		completeMultiTree(periodicbasis2d,new_index2_p,Coeffs_Periodic,additional_indizes_pp,0,true);
		cout << "Extended Coeffs: " << Coeffs_Periodic.size() << endl << Coeffs_Periodic << endl;
		cout << "Added Indizes: " << additional_indizes_pp.size() << additional_indizes_pp << endl << endl;


		cout << "===== PERIODIC x  INTERVAL : j = " << j << " ======= " << endl;
		getSparseGridIndexSet(periodicintervalbasis2d,Lambda_PeriodicInterval,j,0.2);
		FillWithZeros(Lambda_PeriodicInterval, Coeffs_PeriodicInterval);

		cout << "Sparse Grid 2D: " << Lambda_PeriodicInterval.size() << Lambda_Periodic << endl;
	    Index2D new_index1_pi(index1_x_p,index1_y_i);
		cout << "Adding Index " << new_index1_pi << " with support "
			 << periodicintervalbasis2d.first.generator(index1_x_p.xtype).support(index1_x_p.j, index1_x_p.k)
			 << " x " << periodicintervalbasis2d.second.generator(index1_y_i.xtype).support(index1_y_i.j, index1_y_i.k) << endl;
		completeMultiTree(periodicintervalbasis2d,new_index1_pi,Coeffs_PeriodicInterval,0,true);
		cout << "Extended Coeffs: " << Coeffs_PeriodicInterval.size() << endl << Coeffs_PeriodicInterval << endl;

		/*cout << "X1 Alignment " << endl;
	    XOneAlignedCoefficients x1aligned(6151,193);
	    x1aligned.align(Coeffs_PeriodicInterval,J+j0);
	    for (XOneAlignedCoefficients::const_map_prindex_it it=x1aligned.map.begin();
	                                                            it!=x1aligned.map.end(); ++it) {
	    	cout << (*it).first << (*it).second << endl;
	    }

		cout << "X2 Alignment " << endl;
	    XTwoAlignedCoefficients x2aligned(6151,193);
	    x2aligned.align(Coeffs_PeriodicInterval,J+j0);
	    for (XTwoAlignedCoefficients::const_map_prindex_it it=x2aligned.map.begin();
	                                                            it!=x2aligned.map.end(); ++it) {
	    	cout << (*it).first << (*it).second << endl;
	    }*/

		cout << "--------    2nd Completion -- return added indizes   ---------------" << endl << endl;

		IndexSet<Index2D> additional_indizes_pi;
	    Index2D new_index2_pi(index2_x_p,index2_y_i);
		cout << "Adding Index " << new_index2_pi << " with support "
			 << periodicintervalbasis2d.first.generator(index2_x_p.xtype).support(index2_x_p.j, index2_x_p.k)
			 << " x " << periodicintervalbasis2d.second.generator(index2_y_i.xtype).support(index2_y_i.j, index2_y_i.k) << endl;
		completeMultiTree(periodicintervalbasis2d,new_index2_pi,Coeffs_PeriodicInterval,additional_indizes_pi,0,true);
		cout << "Extended Coeffs: " << Coeffs_PeriodicInterval.size() << endl << Coeffs_PeriodicInterval << endl;
		cout << "Added Indizes: " << additional_indizes_pi.size() << additional_indizes_pi << endl << endl;


		cout << "===== INTERVAL x PERIODIC : j = " << j << " ======= " << endl;
		getSparseGridIndexSet(intervalperiodicbasis2d,Lambda_IntervalPeriodic,j,0.2);
		FillWithZeros(Lambda_IntervalPeriodic, Coeffs_IntervalPeriodic);

		cout << "Sparse Grid 2D: " << Lambda_IntervalPeriodic.size() << endl;// << Lambda_Periodic << endl;
	    Index2D new_index1_ip(index1_x_i,index1_y_p);
		cout << "Adding Index " << new_index1_ip << " with support "
			 << intervalperiodicbasis2d.first.generator(index1_x_i.xtype).support(index1_x_i.j, index1_x_i.k)
			 << " x " << intervalperiodicbasis2d.second.generator(index1_y_p.xtype).support(index1_y_p.j, index1_y_p.k) << endl;
		completeMultiTree(intervalperiodicbasis2d,new_index1_ip,Coeffs_IntervalPeriodic,0,true);
		cout << "Extended Coeffs: " << Coeffs_IntervalPeriodic.size() << endl << Coeffs_IntervalPeriodic << endl;

		cout << "--------    2nd Completion -- return added indizes   ---------------" << endl << endl;

		IndexSet<Index2D> additional_indizes_ip;
	    Index2D new_index2_ip(index2_x_i,index2_y_p);
		cout << "Adding Index " << new_index2_i << " with support "
			 << intervalperiodicbasis2d.first.generator(index2_x_i.xtype).support(index2_x_i.j, index2_x_i.k)
			 << " x " << intervalperiodicbasis2d.second.generator(index2_y_p.xtype).support(index2_y_p.j, index2_y_p.k) << endl;
		completeMultiTree(intervalperiodicbasis2d,new_index2_ip,Coeffs_IntervalPeriodic,additional_indizes_ip,0,true);
		cout << "Extended Coeffs: " << Coeffs_IntervalPeriodic.size() << endl << Coeffs_IntervalPeriodic << endl;
		cout << "Added Indizes: " << additional_indizes_ip.size() << additional_indizes_ip << endl << endl;

		cout << endl << endl;
	}


	return 0;
}


template <typename Basis2D>
void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, T gamma)
{
    int j0_1 = basis.first.j0;
    int j0_2 = basis.second.j0;
    for (long k1=basis.first.mra.rangeI(j0_1).firstIndex(); k1<=basis.first.mra.rangeI(j0_1).lastIndex(); ++k1) {
        for (long k2=basis.second.mra.rangeI(j0_2).firstIndex(); k2<=basis.second.mra.rangeI(j0_2).lastIndex(); ++k2) {
            Index1D row(j0_1,k1,XBSpline);
            Index1D col(j0_2,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0_2+i2-1;
            for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0_1,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (long k2=basis.second.mra.rangeI(j0_2).firstIndex(); k2<=basis.second.mra.rangeI(j0_2).lastIndex(); ++k2) {
        for (int i1=1; i1<=j; ++i1) {
            int j1=j0_1+i1-1;
            for (long k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.second.rangeJ(j1).lastIndex(); ++k1) {
                Index1D row(j1,k1,XWavelet);
                Index1D col(j0_2,k2,XBSpline);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0_1+i1-1;
        for (int i2=1; i1+i2<=j; ++i2) {
            if (T(i1+i2)-gamma*max(i1,i2)>(1-gamma)*j) continue;
            int j2=j0_2+i2-1;
            for (long k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    Lambda.insert(Index2D(row,col));
                }
            }
        }
    }
    return;
}
