/*
 * affinerhs.h
 *
 *  Created on: 08.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINERHS_H_
#define LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINERHS_H_

#include <lawa/methods/adaptive/righthandsides/flexiblecompoundrhs.h>

namespace lawa {

template <typename T, typename Index, typename RHSType, typename ParamType>
class AffineRhs : public FlexibleCompoundRhs<T, Index,RHSType>{

public:
	AffineRhs(ThetaStructure<ParamType>& _thetas,
			  std::vector<RHSType*>& _rhsvec);

	T
	operator()(const Index& index);

	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset);

	void
	set_param(ParamType& mu);

private:

	ThetaStructure<ParamType>& thetas;

	AffineRhs(const AffineRhs& rhs);
};

} // namespace lawa

#include <lawa/methods/rb/righthandsides/affinerhs.tcc>

#endif /* LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINERHS_H_ */
