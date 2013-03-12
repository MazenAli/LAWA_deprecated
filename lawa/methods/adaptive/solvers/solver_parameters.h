/*
 * solver_parameters.h
 *
 *  Created on: 05.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_SOLVER_PARAMETERS_H_
#define LAWA_METHODS_ADAPTIVE_SOLVERS_SOLVER_PARAMETERS_H_

#include <vector>
#include <string>
#include <cstdlib>

namespace lawa {

/**
 * Parameters for the Adaptive Wavelet-Galerkin Method
 *  (= the outer solver), Petrov-Galerkin version
 */
struct AWGM_PG_Parameters{

	double 		tol;
	double 		alpha;
	size_t 		max_its;
	size_t 		max_basissize;
	bool 		reset_res;

	bool		print_info;
	bool 		verbose;
	bool    	plot_solution;
	bool		verbose_extra;
	size_t 		hashmapsize_trial;
	size_t 		hashmapsize_test;
	std::string info_filename;
	std::string plot_filename;

	AWGM_PG_Parameters(double _tol = 5e-03,
					double _alpha = 0.7,
					size_t _max_its = 100,
					size_t _max_basissize = 400000,
					bool _reset_res = false,
					bool _print_info = true,
					bool _verbose = true,
					bool _plot_solution = false,
					bool _verbose_extra = false,
					size_t _hashmapsize_trial = 10,
					size_t _hashmapsize_test = 10,
					std::string _info_filename = "awgm_cgls_conv_info.txt",
					std::string _plot_filename = "awgm_cgls_u_plot");

	void print();
};

/**
 * Parameters for the Adaptive Wavelet-Galerkin Method
 *  (= the outer solver), Galerkin version
 */
struct AWGM_Parameters{

	double 	tol;
	double 	alpha;
	size_t 	max_its;
	size_t 	max_basissize;

	bool	print_info;
	bool 	verbose;
	bool    plot_solution;
	bool	verbose_extra;
	size_t 	hashmapsize;
	std::string info_filename;
	std::string plot_filename;

	AWGM_Parameters(double _tol = 5e-03,
					double _alpha = 0.7,
					size_t _max_its = 100,
					size_t _max_basissize = 400000,
					bool _print_info = true,
					bool _verbose = true,
					bool _plot_solution = false,
					bool _verbose_extra = false,
					size_t _hashmapsize = 10,
					std::string _info_filename = "awgm_cg_conv_info.txt",
					std::string _plot_filename = "awgm_cg_u_plot");

	void print();
};

/**
 * Parameters for the inner solver (cg/cgls)
 */
struct IS_Parameters{
	bool adaptive_tol;
	size_t max_its;
	double init_tol;
	double res_reduction;
	double absolute_tol;

	bool verbose;

	IS_Parameters(bool _adaptive_tol = true,
				  size_t _max_its = 100,
				  double _init_tol = 0.001,
				  double _res_reduction = 0.01,
				  double _absolute_tol = 1e-8,
				  bool _verbose = true);

	void print();
};

/**
 * Gathers information that is interesting during a cgls-solve
 * and which can later be printed out
 */
struct AWGM_PG_Information{
	std::vector<double> awgm_res, awgm_resNE,
						sizeLambdaTrial, sizeLambdaTest,
						sizeLambdaResNE, sizeLambdaRes,
						cgls_its;

	void print(const char* filename = "awgm_cgls_conv_info.txt");
};

/**
 * Gathers information that is interesting during a cg-solve
 * and which can later be printed out
 */
struct AWGM_Information{
	std::vector<double> awgm_res, sizeLambda,
						sizeLambdaRes, cg_its;

	void print(const char* filename = "awgm_cg_conv_info.txt");
};

} // namespace lawa

#include <lawa/methods/adaptive/solvers/solver_parameters.tcc>

#endif /* LAWA_METHODS_ADAPTIVE_SOLVERS_SOLVER_PARAMETERS_H_ */
