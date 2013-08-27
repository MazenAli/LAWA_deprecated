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

typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XOne>             XOneAlignedCoefficients;
typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XTwo>             XTwoAlignedCoefficients;


void test_periodic_mt_completion(const PeriodicBasis& basis, int J);
void test_periodic_mt_completion_returningIndizes(const PeriodicBasis& basis, int J);


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

    PeriodicBasis periodicbasis(d, d, j0);
    PeriodicRefinementBasis &periodicrefinementbasis = periodicbasis.refinementbasis;

    //test_periodic_mt_completion(periodicbasis, j0+4);
    test_periodic_mt_completion_returningIndizes(periodicbasis, j0+4);

	return 0;
}

void test_periodic_mt_completion(const PeriodicBasis& basis, int J)
{
    cout << " ******** 1D MT Completion Periodic Basis  **********" << endl;
    for (int j_wavelet=basis.j0; j_wavelet<=J; ++j_wavelet) {
        for (long k_wavelet= basis.rangeJ(j_wavelet).firstIndex();
                  k_wavelet<=basis.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {

            Coefficients<Lexicographical,T,Index1D> coeffs;
    		completeMultiTree(basis,Index1D(j_wavelet, k_wavelet, XWavelet),coeffs,true);

            cout << "Wavelet (" << j_wavelet << "," << k_wavelet << "): " << coeffs << endl;


    		ofstream file("1d_mt_completion_periodic.txt");
    		if(file.is_open()){
    			cout << "Supports: " << endl;
    			for(auto& c : coeffs){
    	    		PeriodicSupport<T> supp = basis.generator(c.first.xtype).support(c.first.j, c.first.k);

    	    		cout << "   " << (c.first.xtype==XWavelet?"W(":"S(")<< c.first.j << "," << c.first.k << ") : " << supp << endl;

        			// set obj rect from x1,y1 to x2,y2
    	    		T y1 = c.first.xtype == XWavelet ? c.first.j : c.first.j - 1;
    	    		T y2 = y1 + 1;

    	    		if(supp.gaplength() > 0){
            			file << "set obj rect from " << supp.l1 <<"," << y1 << " to " << supp.li1 << "," << y2 << endl;
            			file << "set obj rect from " << supp.li2 <<"," << y1 << " to " << supp.l2 << "," << y2 << endl;
    	    		}
    	    		else{
            			file << "set obj rect from " << supp.l1 <<"," << y1 << " to " << supp.l2 << "," << y2 << endl;
    	    		}
    			}
    			file << endl;
    			file << "plot [0:1][0:"<< j_wavelet+2 << "] 1 " << endl;
    			file.close();
    		}
    		else{
    			cerr << "Couldn't open file for writing!" << endl;
    			cout << "Coeffs: " << coeffs << endl;
    		}

            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}

void test_periodic_mt_completion_returningIndizes(const PeriodicBasis& basis, int J)
{
    cout << " ******** 1D MT Completion Periodic Basis (returning indizes)  **********" << endl;
    for (int j_wavelet=basis.j0; j_wavelet<=J; ++j_wavelet) {
        for (long k_wavelet= basis.rangeJ(j_wavelet).firstIndex();
                  k_wavelet<=basis.rangeJ(j_wavelet).lastIndex(); ++k_wavelet) {

            Coefficients<Lexicographical,T,Index1D> coeffs;
            IndexSet<Index1D> diff_v;

    		completeMultiTree(basis,Index1D(j_wavelet, k_wavelet, XWavelet),coeffs,diff_v,true);

            cout << "Wavelet (" << j_wavelet << "," << k_wavelet << "): " << coeffs << endl;
            cout << "Indizes:" << diff_v << endl;


    		ofstream file("1d_mt_completion_periodic.txt");
    		if(file.is_open()){
    			cout << "Supports: " << endl;
    			for(auto& c : coeffs){
    	    		PeriodicSupport<T> supp = basis.generator(c.first.xtype).support(c.first.j, c.first.k);

    	    		cout << "   " << (c.first.xtype==XWavelet?"W(":"S(")<< c.first.j << "," << c.first.k << ") : " << supp << endl;

        			// set obj rect from x1,y1 to x2,y2
    	    		T y1 = c.first.xtype == XWavelet ? c.first.j : c.first.j - 1;
    	    		T y2 = y1 + 1;

    	    		if(supp.gaplength() > 0){
            			file << "set obj rect from " << supp.l1 <<"," << y1 << " to " << supp.li1 << "," << y2 << endl;
            			file << "set obj rect from " << supp.li2 <<"," << y1 << " to " << supp.l2 << "," << y2 << endl;
    	    		}
    	    		else{
            			file << "set obj rect from " << supp.l1 <<"," << y1 << " to " << supp.l2 << "," << y2 << endl;
    	    		}
    			}
    			file << endl;
    			file << "plot [0:1][0:"<< j_wavelet+2 << "] 1 " << endl;
    			file.close();
    		}
    		else{
    			cerr << "Couldn't open file for writing!" << endl;
    			cout << "Coeffs: " << coeffs << endl;
    		}

            cout << endl;
            getchar();
        }
    }
    cout << " ************************************************" << endl << endl;
}
