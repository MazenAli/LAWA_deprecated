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

string get_color(int level){
	switch(level){
	case 0:
		return "\"black\"";
	case 1:
		return "\"red\"";
	case 2:
		return "\"green\" ";
	case 3:
		return "\"blue\"";
	case 4:
		return "\"purple\"";
	case 5:
		return "\"magenta\"";
	case 6:
		return "\"dark-green\"";
	case 7:
		return "\"dark-blue\"";
	default:
		return "";
	}
};

void test_periodic_mt_completion(const PeriodicBasis2D& basis, int J);
void test_periodic_mt_completion_returningIndizes(const PeriodicBasis2D& basis, int J);


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
    PeriodicBasis2D periodicbasis2d(periodicbasis,periodicbasis);


    //test_periodic_mt_completion(periodicbasis2d, j0+2);
    test_periodic_mt_completion_returningIndizes(periodicbasis2d, j0+2);

	return 0;
}

void test_periodic_mt_completion(const PeriodicBasis2D& basis, int J)
{
    cout << " ******** 2D MT Completion Periodic Basis  **********" << endl;
    for (int j_wavelet_x=basis.first.j0; j_wavelet_x<=J; ++j_wavelet_x) {
    	for(int j_wavelet_y = basis.second.j0; j_wavelet_y <= j_wavelet_x; ++j_wavelet_y){
			for (long k_wavelet_x= basis.first.rangeJ(j_wavelet_x).firstIndex();
					  k_wavelet_x<=basis.first.rangeJ(j_wavelet_x).lastIndex(); ++k_wavelet_x) {

				Index1D ind_x = Index1D(j_wavelet_x, k_wavelet_x, XWavelet);

				for (long k_wavelet_y= basis.second.rangeJ(j_wavelet_y).firstIndex();
						  k_wavelet_y<=basis.second.rangeJ(j_wavelet_y).lastIndex(); ++k_wavelet_y) {

					Index1D ind_y = Index1D(j_wavelet_y, k_wavelet_y, XWavelet);

					Coefficients<Lexicographical,T,Index2D> coeffs;
					completeMultiTree(basis,Index2D(ind_x, ind_y),coeffs,0,true);

					cout << "Wavelet (" << j_wavelet_x << "," << k_wavelet_x << ") x Wavelet ("
										<< j_wavelet_y <<"," << k_wavelet_y << " ): " << coeffs << endl;


					ofstream file("2d_mt_completion_periodic.txt");
					if(file.is_open()){
						file << "set terminal x11" << endl << endl;

						cout << "Supports: " << endl;
						for(auto& c : coeffs){
							PeriodicSupport<T> supp_x = basis.first.generator(c.first.index1.xtype).support(c.first.index1.j, c.first.index1.k);
							PeriodicSupport<T> supp_y = basis.second.generator(c.first.index2.xtype).support(c.first.index2.j, c.first.index2.k);

							cout << "   " << (c.first.index1.xtype==XWavelet?"W(":"S(")<< c.first.index1.j << "," << c.first.index1.k << ") x "
								 << (c.first.index2.xtype==XWavelet?"W(":"S(")<< c.first.index2.j << "," << c.first.index2.k << ") : " << supp_x << " x "
								 << supp_y << endl;

							// set obj rect from x1,y1 to x2,y2
							int level = std::max(c.first.index1.j, c.first.index2.j);
							if(c.first.index1.xtype==XBSpline && c.first.index2.xtype==XBSpline) level -= 1;

							if(supp_x.gaplength() > 0){
								if(supp_y.gaplength() > 0){
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.l1 << " to " << supp_x.li1 <<"," << supp_y.li1
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.li2 <<"," << supp_y.l1 << " to " << supp_x.l2 <<"," << supp_y.li1
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.li2 << " to " << supp_x.li1 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.li2 <<"," << supp_y.li2 << " to " << supp_x.l2 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
								}
								else{
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.l1 << " to " << supp_x.li1 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.li2 <<"," << supp_y.l1 << " to " << supp_x.l2 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
								}
							}
							else{
								if(supp_y.gaplength() > 0){
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.l1 << " to " << supp_x.l2 <<"," << supp_y.li1
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.li2 << " to " << supp_x.l2 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
								}
								else{
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.l1 << " to " << supp_x.l2 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
								}
							}
						}
						file << endl;
						file << "plot [0:1][0:1] 1 " << endl;
						file << "unset obj" << endl;
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
    	}
    }
    cout << " ************************************************" << endl << endl;
}

void test_periodic_mt_completion_returningIndizes(const PeriodicBasis2D& basis, int J)
{
    cout << " ******** 2D MT Completion Periodic Basis  **********" << endl;
    for (int j_wavelet_x=basis.first.j0; j_wavelet_x<=J; ++j_wavelet_x) {
    	for(int j_wavelet_y = basis.second.j0; j_wavelet_y <= j_wavelet_x; ++j_wavelet_y){
			for (long k_wavelet_x= basis.first.rangeJ(j_wavelet_x).firstIndex();
					  k_wavelet_x<=basis.first.rangeJ(j_wavelet_x).lastIndex(); ++k_wavelet_x) {

				Index1D ind_x = Index1D(j_wavelet_x, k_wavelet_x, XWavelet);

				for (long k_wavelet_y= basis.second.rangeJ(j_wavelet_y).firstIndex();
						  k_wavelet_y<=basis.second.rangeJ(j_wavelet_y).lastIndex(); ++k_wavelet_y) {

					Index1D ind_y = Index1D(j_wavelet_y, k_wavelet_y, XWavelet);

					Coefficients<Lexicographical,T,Index2D> coeffs;
		            IndexSet<Index2D> diff_v;

					completeMultiTree(basis,Index2D(ind_x, ind_y),coeffs,diff_v, 0,true);

					cout << "Wavelet (" << j_wavelet_x << "," << k_wavelet_x << ") x Wavelet ("
										<< j_wavelet_y <<"," << k_wavelet_y << " ): " << coeffs << endl;
		            cout << "Indizes:" << diff_v << endl;


					ofstream file("2d_mt_completion_periodic.txt");
					if(file.is_open()){
						file << "set terminal x11" << endl << endl;

						cout << "Supports: " << endl;
						for(auto& c : coeffs){
							PeriodicSupport<T> supp_x = basis.first.generator(c.first.index1.xtype).support(c.first.index1.j, c.first.index1.k);
							PeriodicSupport<T> supp_y = basis.second.generator(c.first.index2.xtype).support(c.first.index2.j, c.first.index2.k);

							cout << "   " << (c.first.index1.xtype==XWavelet?"W(":"S(")<< c.first.index1.j << "," << c.first.index1.k << ") x "
								 << (c.first.index2.xtype==XWavelet?"W(":"S(")<< c.first.index2.j << "," << c.first.index2.k << ") : " << supp_x << " x "
								 << supp_y << endl;

							// set obj rect from x1,y1 to x2,y2
							int level = std::max(c.first.index1.j, c.first.index2.j);
							if(c.first.index1.xtype==XBSpline && c.first.index2.xtype==XBSpline) level -= 1;

							if(supp_x.gaplength() > 0){
								if(supp_y.gaplength() > 0){
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.l1 << " to " << supp_x.li1 <<"," << supp_y.li1
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.li2 <<"," << supp_y.l1 << " to " << supp_x.l2 <<"," << supp_y.li1
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.li2 << " to " << supp_x.li1 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.li2 <<"," << supp_y.li2 << " to " << supp_x.l2 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
								}
								else{
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.l1 << " to " << supp_x.li1 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.li2 <<"," << supp_y.l1 << " to " << supp_x.l2 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
								}
							}
							else{
								if(supp_y.gaplength() > 0){
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.l1 << " to " << supp_x.l2 <<"," << supp_y.li1
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.li2 << " to " << supp_x.l2 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
								}
								else{
									file << "set obj rect from " << supp_x.l1 <<"," << supp_y.l1 << " to " << supp_x.l2 <<"," << supp_y.l2
										  << " fs transparent pattern " << level << " fc rgbcolor " << get_color(level) << endl;
								}
							}
						}
						file << endl;
						file << "plot [0:1][0:1] 1 " << endl;
						file << "unset obj" << endl;
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
    	}
    }
    cout << " ************************************************" << endl << endl;
}

