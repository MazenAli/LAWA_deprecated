/*
 * rb_parameters.h
 *
 *  Created on: 12.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_RB_PARAMETERS_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_RB_PARAMETERS_H_

#include <array>

namespace lawa {

/* Helper struct to get size of an array
 */
template<typename ParamType>
struct ParamInfo {};

template<typename T,std::size_t N>
struct ParamInfo<std::array<T,N> >
{
	static const std::size_t dim = N;

	static void
	print(std::array<T,N> param);
};

/* Parameters for an RB Training
 *
 */
template<typename ParamType>
struct RB_Greedy_Parameters{

	typedef typename std::array<size_t,ParamInfo<ParamType>::dim> intArray;
	double 		tol;
	size_t 		Nmax;
	ParamType 	min_param;
	ParamType 	max_param;
	intArray	nb_training_params;

	bool 		print_info;
	std::string print_file;
	bool 		verbose;
	bool 		write_during_training;
	std::string trainingdata_folder;
	bool 		print_paramset;
	bool 		erase_snapshot_params;

	RB_Greedy_Parameters(double _tol = 1e-2,
						size_t _Nmax = 20,
						ParamType _min_param = ParamType(),
						ParamType _max_param = ParamType(),
						intArray  _training_params_per_dim = intArray(),
						bool _print_info = true,
						std::string _print_file = "greedy_info.txt",
						bool _verbose = true,
						bool _write_during_training = true,
						std::string _trainingdata_folder = "training_data",
						bool _print_paramset = false,
						bool _erase_snapshot_params = false);

	void print();
};

template<typename ParamType>
struct RB_Greedy_Information{
	std::vector<double> 	 			   greedy_errors;
	std::vector<ParamType>	 			   snapshot_params;
	std::vector<std::size_t>		 	   u_size;				// Dim 1 x n_bf
	std::vector<std::size_t>   			   repr_f_size;		    // Dim 1 x Q_f
	std::vector<std::vector<std::size_t> > repr_a_size;			// Dim n_bf x Q_a

	void print(const char* filename = "greedy_info.txt");
};

/* Parameters for a RB solution
 */
template<typename ParamType>
struct RB_Parameters {

	SolverCall  call;
	ParamType 	ref_param;

	bool verbose;

	RB_Parameters(SolverCall _call = call_cg,
				  ParamType _ref_param = ParamType(),
				  bool _verbose = true);

	void print();
};


} // namespace lawa

#include <lawa/methods/rb/datastructures/rb_parameters.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_RB_PARAMETERS_H_ */
