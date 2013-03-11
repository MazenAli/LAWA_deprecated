/*
 * affinelocaloperator.h
 *
 *  Created on: 07.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_OPERATORS_AFFINELOCALOPERATOR_H_
#define LAWA_METHODS_RB_OPERATORS_AFFINELOCALOPERATOR_H_

#include <lawa/methods/adaptive/operators/localoperators/flexiblecompoundlocaloperator.h>
#include <lawa/methods/rb/datastructures/thetastructure.h>

namespace lawa {

template <typename Index, typename LocalOperatorType, size_t PDim>
class AffineLocalOperator : public FlexibleCompoundLocalOperator<Index,LocalOperatorType>{

public:
	typedef typename LocalOperatorType::T T;

	AffineLocalOperator(const ThetaStructure<T,PDim>& _thetas, std::vector<LocalOperatorType*>& _localops);

private:

	const ThetaStructure<T,PDim>& thetas;

};

} // namespace lawa

#include <lawa/methods/rb/operators/affinelocaloperator.tcc>

#endif /* LAWA_METHODS_RB_OPERATORS_AFFINELOCALOPERATOR_H_ */
