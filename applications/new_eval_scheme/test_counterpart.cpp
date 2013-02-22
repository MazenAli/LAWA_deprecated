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

template <typename Basis2D_Origin, typename Basis2D_Target>
void
getCounterpart(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

template <typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

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
    PeriodicBasis periodicbasis(d, d, j0);
    IntervalBasis intervalbasis_bc(d, d, j0);
    intervalbasis_bc.enforceBoundaryCondition<DirichletBC>();
    IntervalBasis intervalbasis(d,d,j0);

    IntervalBasis2D 		intervalintervalbcbasis2d(intervalbasis,intervalbasis_bc);
    PeriodicIntervalBasis2D periodicintervalbcbasis2d(periodicbasis,intervalbasis_bc);

	
	for(int j = 0; j <= J-j0; ++j){

		// Indizes
		Index1D index1_x_p(j0+j+1,1,XWavelet);
	    Index1D index1_y_p(j0+j+1,4,XWavelet);
	    Index1D index1_x_i(j0+j+1,2,XWavelet);
	    Index1D index1_y_i(j0+j+1,3,XWavelet);

	    // Construct initial sets by completing the indizes

		IndexSet<Index2D>		indexset_pi;
		Coefficients<Lexicographical,T,Index2D> coeffs_pi, cp_coeffs_ii, ext_coeffs_ii ;

		cout << "===== PERIODIC x INTERVAL_BC: j = " << j << " ======= " << endl << endl;
		Index2D index1_pi(index1_x_p, index1_y_i);
		completeMultiTree(periodicintervalbcbasis2d,index1_pi,coeffs_pi,0,true);
		indexset_pi = supp(coeffs_pi);
		cout << "Initial Set v: " << indexset_pi.size() << indexset_pi << endl;

		getCounterpart(periodicintervalbcbasis2d, intervalintervalbcbasis2d, indexset_pi, cp_coeffs_ii);
		cout << "Counterpart in interval/intervalbc: " << cp_coeffs_ii.size() << endl;// << cp_coeffs_ii << endl;

		getStableExpansion(periodicintervalbcbasis2d, intervalintervalbcbasis2d, indexset_pi, ext_coeffs_ii);
		cout << "Extension in interval/intervalbc: " << ext_coeffs_ii.size() << endl; // << ext_coeffs_ii << endl;

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

// To a given indexset in basis_origin, find a corresponding indexset in basis_target.
// This is realized by taking the cone consisting of neighbours of all bfs in indexset_origin
// We have to make sure that this is a MT!!
template <typename Basis2D_Origin, typename Basis2D_Target>
void
getCounterpart(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{
    long k_first_x, k_last_x, k_first_y, k_last_y;
    int j_x, j_y;
    for(IndexSet<Index2D>::const_iterator it = indexset_origin.begin(); it != indexset_origin.end(); ++it){

		// Get neighbours for index_x
		if((*it).index1.xtype==XBSpline){
			basis_origin.first.getScalingNeighborsForScaling((*it).index1.j, (*it).index1.k,basis_target.first,
																j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j);
		}
		else{
			basis_origin.first.getWaveletNeighborsForWavelet((*it).index1.j, (*it).index1.k,basis_target.first,
																j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j);
		}

		// Get neighbours for index_y
		if((*it).index2.xtype==XBSpline){
			basis_origin.second.getScalingNeighborsForScaling((*it).index2.j, (*it).index2.k, basis_target.second,
																j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j);
		}
		else{
			basis_origin.second.getWaveletNeighborsForWavelet((*it).index2.j, (*it).index2.k, basis_target.second,
																j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j);
		}

		// Insert all combinations into res
			// Outer loop: all indizes in x-direction
		if(k_first_x < k_last_x){
			for(long k_new_x = k_first_x; k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					if((*it).index2.xtype==XBSpline){
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
				}
			}
		}
		else{
			if((*it).index1.xtype==XBSpline){
				for(long k_new_x = k_first_x; k_new_x <= basis_target.first.mra.rangeI(j_x).lastIndex(); ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
				for(long k_new_x = basis_target.first.mra.rangeI(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
			}
			else{
				for(long k_new_x = k_first_x; k_new_x <= basis_target.first.rangeJ(j_x).lastIndex(); ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
				for(long k_new_x = basis_target.first.rangeJ(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
			}
		}
    }
}

// To a given indexset in basis_origin, find a corresponding indexset in basis_target,
// so that the system is "stable" (or at least A^T A is approximated well enough).
// This is realized by taking the cone consisting of neighbours of all bfs in indexset_origin
// plus the HigherWavelet-Neighbours.
template <typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{

	// First we take the counterpart cone
	getCounterpart(basis_origin, basis_target, indexset_origin, coeffs_target);

	// Then we insert all HigherWaveletNeighbours
    long k_first_x, k_last_x, k_first_y, k_last_y;
    int j_x, j_y;
    for(IndexSet<Index2D>::const_iterator it = indexset_origin.begin(); it != indexset_origin.end(); ++it){

		// Get neighbours for index_x
		if((*it).index1.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.first.getWaveletNeighborsForScaling((*it).index1.j, (*it).index1.k,basis_origin.first,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index1.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_x = basis_target.first.rangeJ(j_s+1).lastIndex();
			k_last_x = basis_target.first.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.first.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				for(long k = basis_origin.first.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				assert(j_x == (*it).index1.j+1);
			}
		}
		else{
			basis_origin.first.getHigherWaveletNeighborsForWavelet((*it).index1.j, (*it).index1.k,basis_target.first,
																	j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j+1);

		}

		// Get neighbours for index_y
		if((*it).index2.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.second.getWaveletNeighborsForScaling((*it).index2.j, (*it).index2.k,basis_origin.second,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index2.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_y = basis_target.second.rangeJ(j_s+1).lastIndex();
			k_last_y = basis_target.second.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.second.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
				for(long k = basis_origin.second.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
				assert(j_y == (*it).index2.j+1);
			}
		}
		else{
			basis_origin.second.getHigherWaveletNeighborsForWavelet((*it).index2.j, (*it).index2.k, basis_target.second,
																	j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j+1);
		}

		// Insert all combinations into res
			// Outer loop: all indizes in x-direction
		if(k_first_x < k_last_x){
			for(long k_new_x = k_first_x; k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
					for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
			}
		}
		else{
			for(long k_new_x = k_first_x; k_new_x <= basis_target.first.rangeJ(j_x).lastIndex(); ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					if((*it).index2.xtype==XBSpline){
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
				}
			}
			for(long k_new_x = basis_target.first.rangeJ(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
					for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
			}
		}
    }
}

