/*
 * abstractlocaloperator2d.h
 *
 *  Created on: 25.02.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_ABSTRACTLOCALOPERATOR2D_H_
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_ABSTRACTLOCALOPERATOR2D_H_

namespace lawa {

template <typename _T>
struct AbstractLocalOperator2D {

    typedef _T T;

    virtual void
    eval(const Coefficients<Lexicographical,T,Index2D> &input,
    		   Coefficients<Lexicographical,T,Index2D> &output) = 0;

};


} // namespace lawa

#endif /* LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_ABSTRACTLOCALOPERATOR2D_H_ */
