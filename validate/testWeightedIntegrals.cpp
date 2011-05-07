/*
 *  testWeightedIntegrals.cpp
 *  Xcode-LAWA
 *
 *  Created by Kristina Steih on 27.04.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

// Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>		IntervalBasis;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> > 	DenseVectorT;

T
weight(T x)
{
	return (x < 0.5) ? 1 : 0;
}

T
weight2(T x)
{
    return 1.;
}

T
exponentialweight(T x)
{
    return exp(-4*fabs(x));
}

int main() {
    cout.precision(16);
	int d  = 2;
    int d_ = 2;
    int j0 = 2;
	IntervalBasis basis1(d, d_, j0);
	IntervalBasis basis2(d, d_, j0);
    /*
    DenseVectorT singpts(3);
    singpts = 0., 0.5, 1.;
	Function<T> weightFct(weight, singpts);    

	// 1D: works fine

    IntegralF<Gauss, IntervalBasis> integral1(weightFct, basis1);

    cout << integral1(j0, 1, XBSpline, 0) << endl;

    
    // 2D: compiler error
    IntegralF<Gauss, IntervalBasis> integral2(weightFct, basis1, basis2);
    
    cout << integral2(j0, 1, XBSpline, 0, j0, 1, XBSpline, 0) << endl;
	
    std::cout << "Check for correct scaling in IntegralF" << std::endl;
    DenseVectorT weight2_singpts;
    Function<T> weightFct2(weight2, weight2_singpts);
    IntegralF<Gauss, IntervalBasis> integral_ref(weightFct2, basis1, basis1);
    IntegralF<Gauss, IntervalBasis> integral_test(weightFct2, basis1, basis1, -1., 1.);

    cout << integral_ref(j0, 1, XBSpline, 0, j0, 1, XBSpline, 0) << endl;
    cout << integral_test(j0, 1, XBSpline, 0, j0, 1, XBSpline, 0) << endl;

    cout << integral_ref(j0+3, 1, XWavelet, 0, j0, 1, XBSpline, 0) << endl;
    cout << integral_test(j0+3, 1, XWavelet, 0, j0, 1, XBSpline, 0) << endl;

    cout << integral_ref(j0+3, 1, XWavelet, 0, j0+2, 1, XWavelet, 0) << endl;
    cout << integral_test(j0+3, 1, XWavelet, 0, j0+2, 1, XWavelet, 0) << endl;
    */
    std::cout << "Check for weighted integrals in IntegralF" << std::endl;
    DenseVectorT exponentialweight_singpts(1);
    exponentialweight_singpts = 0.;
    Function<T> exponentialweightFct(exponentialweight, exponentialweight_singpts);
    IntegralF<Gauss, IntervalBasis> integral_exp_num(exponentialweightFct, basis1, basis1,-1.,1.);
    integral_exp_num.quadrature.setOrder(4);
    IntegralF<ExpWeighted, IntervalBasis> integral_exp_exct(exponentialweightFct, basis1, basis1,-1.,1.);
    /*
    Timer time;
    int count=0;
    T tmp1 = 0.;
    time.start();
    for (int k=basis1.rangeJ(8).firstIndex(); k<=basis1.rangeJ(8).lastIndex();++k) {
        tmp1 += integral_exp_num(8, k, XWavelet, 1, 8, k, XWavelet, 1);
        ++count;
    }
    time.stop();
    cout << "Result: " << tmp1 << ", time elapsed Legendre: " << time.elapsed() << endl;
    cout << "count = " << count << endl;
    count = 0;
    T tmp2 = 0.;
    time.start();
    for (int k=basis1.rangeJ(8).firstIndex(); k<=basis1.rangeJ(8).lastIndex();++k) {
        tmp2 += integral_exp_exct(8, k, XWavelet, 1, 8, k, XWavelet, 1);
        ++count;
    }
    time.stop();
    cout << "Result: " << tmp2 << ", time elapsed exact formula: " << time.elapsed() << endl;
    cout << "count = " << count << endl;
    */
    cout << integral_exp_num(j0+5, 2, XWavelet, 1, j0+5, 3, XWavelet, 1) << endl;
    cout << integral_exp_exct(j0+5, 2, XWavelet, 1, j0+5, 3, XWavelet, 1) << endl;
/*
    cout << integral_exp_num(j0+10, 1, XWavelet, 1, j0+10, 1, XWavelet, 1) << endl;
    cout << integral_exp_exct(j0+10, 1, XWavelet, 1, j0+10, 1, XWavelet, 1) << endl;
*/


	return 0;
}
