/*
 * residualrhs.h
 *
 *  Created on: 30.08.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_RIGHTHANDSIDES_RESIDUALRHS_H_
#define LAWA_METHODS_RB_RIGHTHANDSIDES_RESIDUALRHS_H_

#include <vector>
#include <lawa/methods/adaptive/datastructures/datastructures.h>


namespace lawa {

/** Operator to evaluate the residual Res(v;mu) = f(v;mu) - a(u_N,v;mu)
 *  as right hand side of a linear system.
 */

template <typename Index, typename LHSType, typename RHSType, typename ParamType, typename DataType>
class ResidualRhs {

public:

	typedef typename LHSType::T T;

	ResidualRhs(LHSType& _lhs_bilform, RHSType& _rhs_fct);

	T
    operator()(const Index &index);

	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset);

	void
	set_param(ParamType& mu);

	void
	set_active_u(DataType const* u); // supposed to be the reconstructed approximate solution u_RB

private:

    LHSType&				lhs_bilform;
    RHSType&				rhs_fct;

    ResidualRhs(const ResidualRhs& rhs);
};

} // namespace lawa

#include <lawa/methods/rb/righthandsides/residualrhs.tcc>

#endif /* LAWA_METHODS_RB_RIGHTHANDSIDES_RESIDUALRHS_H_ */
