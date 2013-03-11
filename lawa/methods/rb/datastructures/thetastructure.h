/*
 * thetastructure.h
 *
 *  Created on: 07.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_

#include <vector>

namespace lawa {

/* ThetaStructure:
 * 		Manage a vector of theta functions.
 * 		Helpful in order to store them in only one place.
 */
template<typename T,size_t PDim>
class ThetaStructure{

public:
    typedef T (*ThetaFct)(const std::array<T,PDim>& params); // Argumente -> eher auch RBThetaData-Objekt?

    ThetaStructure();

    ThetaStructure(const std::vector<ThetaFct>& _thetas);

    unsigned int
    size();

private:
	std::vector<ThetaFct> thetas;
};

} // namespace lawa

#include <lawa/methods/rb/datastructures/thetastructure.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_ */
