/*
 * flexiblebilformrhs.h
 *
 *  Created on: 11.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_RIGHTHANDSIDES_FLEXIBLEBILFORMRHS_H_
#define LAWA_METHODS_RB_RIGHTHANDSIDES_FLEXIBLEBILFORMRHS_H_

#include <vector>
#include <lawa/methods/adaptive/datastructures/datastructures.h>


namespace lawa {

template <typename Index, typename LocalOperatorType>
class FlexibleBilformRhs {

public:

	typedef typename LocalOperatorType::T T;

	FlexibleBilformRhs(std::vector<LocalOperatorType*>& _bilformvec);

	virtual
	T
    operator()(const Index &index);

	virtual
	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset);

	virtual
	void
	set_active_comp(int i);

	virtual
	void
	set_active_u(Coefficients<Lexicographical,T,Index> const* _u);

protected:

    std::vector<LocalOperatorType*>& 				bilformvec;
    std::vector<int> 								active_comp;

    Coefficients<Lexicographical,T,Index> const* 	active_u;

	FlexibleBilformRhs(const FlexibleBilformRhs& rhs);
};

} // namespace lawa

#include <lawa/methods/rb/righthandsides/flexiblebilformrhs.tcc>

#endif /* LAWA_METHODS_RB_RIGHTHANDSIDES_FLEXIBLEBILFORMRHS_H_ */
