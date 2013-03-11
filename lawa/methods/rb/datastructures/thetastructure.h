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
template<typename T,size_t PDim>
class ThetaStructure{

public:
    typedef T (*ThetaFct)(const std::array<T,PDim>& params); // Argumente -> eher auch RBThetaData-Objekt?

    ThetaStructure();

    ThetaStructure(const std::vector<ThetaFct>& _thetas);

    unsigned int
    size() const;

    void
    set_current_param(const std::array<T, PDim>& _param);

    std::array<T, PDim>&
    get_current_param();

    T
    eval(int i, std::array<T,PDim>& mu);

    T
    eval(int i);

private:
	std::vector<ThetaFct> 	thetas;
	std::array<T,PDim>*    	current_param;

	ThetaStructure(const ThetaStructure& thetastructure);
};

} // namespace lawa

#include <lawa/methods/rb/datastructures/thetastructure.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_ */
