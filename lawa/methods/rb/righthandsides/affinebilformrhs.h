/*
 * affinebilformrhs.h
 *
 *  Created on: 11.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINEBILFORMRHS_H_
#define LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINEBILFORMRHS_H_

#include <vector>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

#include <lawa/methods/rb/righthandsides/flexiblebilformrhs.h>

namespace lawa {

template <typename Index, typename LocalOperatorType, typename ParamType>
class AffineBilformRhs : public FlexibleBilformRhs<Index,LocalOperatorType> {

public:

	typedef typename LocalOperatorType::T T;

	AffineBilformRhs(ThetaStructure<ParamType>& _thetas, std::vector<LocalOperatorType*>& _bilformvec);

	T
    operator()(const Index &index);

	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset);

	void
	set_param(ParamType& mu);

private:

	ThetaStructure<ParamType>& thetas;

	AffineBilformRhs(const AffineBilformRhs& rhs);
};

} // namespace lawa

#include <lawa/methods/rb/righthandsides/affinebilformrhs.tcc>

#endif /* LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINEBILFORMRHS_H_ */
