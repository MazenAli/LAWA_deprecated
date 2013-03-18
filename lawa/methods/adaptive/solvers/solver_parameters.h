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
	std::size_t 		max_its;
	std::size_t 		max_basissize;
	bool 		reset_res;

	bool		print_info;
	bool 		verbose;
	bool    	plot_solution;
	bool		verbose_extra;
	std::size_t 		hashmapsize_trial;
	std::size_t 		hashmapsize_test;
	std::string info_filename;
	std::string plot_filename;

	AWGM_PG_Parameters(double _tol = 5e-03,
					double _alpha = 0.7,
					std::size_t _max_its = 100,
					std::size_t _max_basissize = 400000,
					bool _reset_res = false,
					bool _print_info = true,
					bool _verbose = true,
					bool _plot_solution = false,
					bool _verbose_extra = false,
					std::size_t _hashmapsize_trial = 10,
					std::size_t _hashmapsize_test = 10,
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
	std::size_t 	max_its;
	std::size_t 	max_basissize;

	bool	print_info;
	bool 	verbose;
	bool    plot_solution;
	bool	verbose_extra;
	std::size_t 	hashmapsize;
	std::string info_filename;
	std::string plot_filename;

	AWGM_Parameters(double _tol = 5e-03,
					double _alpha = 0.7,
					std::size_t _max_its = 100,
					std::size_t _max_basissize = 400000,
					bool _print_info = true,
					bool _verbose = true,
					bool _plot_solution = false,
					bool _verbose_extra = false,
					std::size_t _hashmapsize = 10,
					std::string _info_filename = "awgm_cg_conv_info.txt",
					std::string _plot_filename = "awgm_cg_u_plot");

	void print();
};

/**
 * Parameters for the inner solver (cg/cgls)
 */
struct IS_Parameters{
	bool 			adaptive_tol;
	std::size_t 	max_its;
	double 			init_tol;
	double 			res_reduction;
	double 			absolute_tol;

	bool 			verbose;

	IS_Parameters(bool _adaptive_tol = true,
				  std::size_t _max_its = 100,
				  double _init_tol = 0.001,
				  double _res_reduction = 0.01,
				  double _absolute_tol = 1e-8,
				  bool _verbose = true);

	void print();
};

/**
 * Parameters for a indexset-based Wavelet-Galerkin Method
 *  (= the outer solver)
 */
struct ISWGM_Parameters{

	bool	print_info;
	bool 	verbose;
	bool    plot_solution;
	std::string info_filename;
	std::string plot_filename;

	ISWGM_Parameters(bool _print_info = true,
					bool _verbose = true,
					bool _plot_solution = false,
					std::string _info_filename = "iswgm_conv_info.txt",
					std::string _plot_filename = "iswgm_u_plot");

	void print();
};

/**
 * Gathers information that is interesting during a awgm-solve
 * and which can later be printed out (Petrov-Galerkin version)
 */
struct AWGM_PG_Information{
	std::vector<double> 		awgm_res, awgm_resNE;
	std::vector<std::size_t>	sizeLambdaTrial, sizeLambdaTest,
								sizeLambdaResNE, sizeLambdaRes,
								cgls_its;

	void print(const char* filename = "awgm_cgls_conv_info.txt");

	void reset();
};

/**
 * Gathers information that is interesting during a awgm-solve
 * and which can later be printed out (alerkin version)
 */
struct AWGM_Information{
	std::vector<double> 		awgm_res;
	std::vector<std::size_t>	sizeLambda, sizeLambdaRes,
								cg_its;

	void print(const char* filename = "awgm_cg_conv_info.txt");

	void reset();
};

/**
 * Gathers information that is interesting during a iswgm-solve
 * and which can later be printed out (Galerkin version)
 */
struct ISWGM_Information{
	std::vector<double> 		is_res;

	void print(const char* filename = "iswgm_is_conv_info.txt");

	void reset();
};

/**
 * Gathers information that is interesting during a iswgm-solve
 * and which can later be printed out (Petrov-Galerkin version)
 */
struct ISWGM_PG_Information{
	std::vector<double> 		is_res, is_resNE;

	void print(const char* filename = "iswgm_pg_is_conv_info.txt");

	void reset();
};

} // namespace lawa

#include <lawa/methods/adaptive/solvers/solver_parameters.tcc>

#endif /* LAWA_METHODS_ADAPTIVE_SOLVERS_SOLVER_PARAMETERS_H_ */
