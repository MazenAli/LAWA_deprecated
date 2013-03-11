/*
 * rb_base.h
 *
 *  Created on: 07.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_RB_BASE_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_RB_BASE_H_

#include <vector>
#include <lawa/methods/rb/datastructures/rb_system.h>
#include <lawa/methods/rb/datastructures/thetastructure.h>

namespace lawa {

/* RB_Base:
 * 		General offline computations for the construction of
 * 		a reduced basis model.
 */
template <typename RB_Model, typename TruthModel, typename DataType>
class RB_Base{

public:
	RB_Base(RB_Model& _rb_system, TruthModel& _rb_truth);

private:
	RB_Model& 		rb_system;
	TruthModel& 	rb_truth;

	std::vector<DataType> 	rb_basisfunctions;

};

} // namespace lawa

#include <lawa/methods/rb/datastructures/rb_base.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_RB_BASE_H_ */
