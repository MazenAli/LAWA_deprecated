/*
 * rb_multitree.cpp
 *
 *  Created on: 07.03.2013
 *      Author: ksteih
 */

#include <iostream>
#include <array>
#include <lawa/methods/rb/datastructures/rb_system.h>
#include <lawa/methods/rb/datastructures/thetastructure.h>

using namespace std;
using namespace lawa;

typedef double T;

int main(){
// This code is deprecated
// TODO: update or remove!
/*
	const int p = 2;

	ThetaStructure<T> empty_theta;

	RB_System<T,p> rb_system(empty_theta,empty_theta);

	array<T,p> mu = {1.,2.5};

	rb_system.set_current_param(mu);

	array<T,p> curr_mu = rb_system.get_current_param();
	for(auto& m : curr_mu){
		cout << m << " ";
	}
	cout << endl;

*/
	return 0;
}

