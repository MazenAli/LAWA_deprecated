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

	AffineLocalOperator(ThetaStructure<T,PDim>& _thetas, std::vector<LocalOperatorType*>& _localops);

	void
	eval(const Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av);

	template <typename Preconditioner>
	void
	eval(Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P);

	template <typename RightPrec, typename LeftPrec>
	void
	eval(Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av, RightPrec& rightP, LeftPrec& leftP);

	void
	set_param(std::array<T,PDim>& mu);

private:

	ThetaStructure<T,PDim>& thetas;

};

} // namespace lawa

#include <lawa/methods/rb/operators/affinelocaloperator.tcc>

#endif /* LAWA_METHODS_RB_OPERATORS_AFFINELOCALOPERATOR_H_ */
