/*
 * thetastructure.h
 *
 *  Created on: 07.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_

#include <vector>
#include <array>

namespace lawa {

/* ThetaStructure:
 * 		Manage a vector of theta functions.
 * 		Helpful in order to store them in only one place.
 *
 * 		Addionally a pointer to the current parameter, so that
 * 		affine operators/rhs that only have an operator()(Index)
 * 		can still evaluate the functions
 */
template<typename ParamType>
class ThetaStructure{

public:
	typedef typename ParamType::value_type T;
    typedef T (*ThetaFct)(const ParamType& param); // Argumente -> eher auch RBThetaData-Objekt?

    ThetaStructure();

    ThetaStructure(const std::vector<ThetaFct>& _thetas);

    size_t
    size() const;

    void
    set_param(const ParamType& _param);

    ParamType&
    get_param();

    T
    eval(size_t i, const ParamType& mu) const;

    T
    eval(size_t i) const;

private:
	std::vector<ThetaFct> 	thetas;
	ParamType    			current_param;

	ThetaStructure(const ThetaStructure& thetastructure);
};

} // namespace lawa

#include <lawa/methods/rb/datastructures/thetastructure.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_ */
