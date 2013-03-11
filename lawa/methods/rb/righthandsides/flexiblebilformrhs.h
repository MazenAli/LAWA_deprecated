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

	T
    operator()(const Index &index);

	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset);

	void
	set_active_comp(int i);

	void
	set_active_u(Coefficients<Lexicographical,T,Index>* _u);

private:

    std::vector<LocalOperatorType*>& 		bilformvec;
    std::vector<int> 						active_comp;

    Coefficients<Lexicographical,T,Index>* 	active_u;

	FlexibleBilformRhs(const FlexibleBilformRhs& rhs);
};

} // namespace lawa

#include <lawa/methods/rb/righthandsides/flexiblebilformrhs.tcc>

#endif /* LAWA_METHODS_RB_RIGHTHANDSIDES_FLEXIBLEBILFORMRHS_H_ */