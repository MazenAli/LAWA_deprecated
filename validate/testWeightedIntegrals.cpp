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

int main() {


    cout.precision(16);
	int d  = 2;
    int d_ = 2;
    int j0 = 2;
	IntervalBasis basis1(d, d_, j0);
	IntervalBasis basis2(d, d_, j0);


    DenseVectorT singpts(3);
    singpts = 0., 0.5, 1.;
	Function<T> weightFct(weight, singpts);    

	// 1D: works fine

    IntegralF<Gauss, IntervalBasis> integral1(weightFct, basis1);

    cout << integral1(j0, 1, XBSpline, 0) << endl;

    // 2D: compiler error
    IntegralF<Gauss, IntervalBasis> integral2(weightFct, basis1, basis2);
    
    cout << integral2(j0, 1, XBSpline, 0, j0, 1, XBSpline, 0) << endl;
	

	return 0;
}
