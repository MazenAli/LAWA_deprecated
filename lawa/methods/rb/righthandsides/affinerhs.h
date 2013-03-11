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

template <typename T, typename Index, typename RHSType, size_t PDim>
class AffineRhs : public FlexibleCompoundRhs<T, Index,RHSType>{

public:
	AffineRhs(const ThetaStructure<T,PDim>& _thetas,
			  std::vector<RHSType*>& _rhsvec);

private:

	const ThetaStructure<T,PDim>& thetas;

	AffineRhs(const AffineRhs& rhs);
};

} // namespace lawa

#include <lawa/methods/rb/righthandsides/affinerhs.tcc>

#endif /* LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINERHS_H_ */
