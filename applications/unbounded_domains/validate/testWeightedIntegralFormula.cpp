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
typedef Basis<T,Primal,Interval,Dijkema>        IntervalBasis;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

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
    IntervalBasis basis(d, d_, j0);

    DenseVectorT exponentialweight_singpts(1);
    exponentialweight_singpts = 0.;
    Function<T> exponentialweightFct(exponentialweight, exponentialweight_singpts);
    IntegralF<Gauss, IntervalBasis> integral_exp_num(exponentialweightFct, basis, basis,-2.,1.);
    integral_exp_num.quadrature.setOrder(5);
    IntegralF<ExpWeighted, IntervalBasis> integral_exp_exct(exponentialweightFct, basis, basis,-2.,1.);


    Timer time;
    int count=0;
    T tmp1 = 0.;
    time.start();
    for (int k=basis.rangeJ(12).firstIndex(); k<=basis.rangeJ(12).lastIndex();++k) {
        tmp1 += integral_exp_num(12, 5, XWavelet, 1, 12, k, XWavelet, 1);
        ++count;
    }
    time.stop();
    cout << "Result: " << tmp1 << ", time elapsed Legendre: " << time.elapsed() << endl;
    cout << "count = " << count << endl;
    count = 0;
    T tmp2 = 0.;
    time.start();
    for (int k=basis.rangeJ(12).firstIndex(); k<=basis.rangeJ(12).lastIndex();++k) {
        tmp2 += integral_exp_exct(12, 5, XWavelet, 1, 12, k, XWavelet, 1);
        ++count;
    }
    time.stop();
    cout << "Result: " << tmp2 << ", time elapsed exact formula: " << time.elapsed() << endl;
    cout << "count = " << count << endl;


    return 0;
}
