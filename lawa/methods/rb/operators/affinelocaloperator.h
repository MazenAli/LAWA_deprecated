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

template <typename Index, typename LocalOperatorType, typename ParamType>
class AffineLocalOperator : public FlexibleCompoundLocalOperator<Index,LocalOperatorType>{

public:
	typedef typename LocalOperatorType::T T;

	AffineLocalOperator(ThetaStructure<ParamType>& _thetas, std::vector<LocalOperatorType*>& _localops);

	void
	eval(const Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av);

	void
	eval(size_t i, const Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av, bool eval_mu = false);

	template <typename Preconditioner>
	void
	eval(Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P);

	template <typename RightPrec, typename LeftPrec>
	void
	eval(Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av, RightPrec& rightP, LeftPrec& leftP);

	void
	set_param(ParamType& mu);

private:

	ThetaStructure<ParamType>& thetas;

};

} // namespace lawa

#include <lawa/methods/rb/operators/affinelocaloperator.tcc>

#endif /* LAWA_METHODS_RB_OPERATORS_AFFINELOCALOPERATOR_H_ */
