#include <iostream>
#include <fstream>
#include <iomanip>

namespace lawa {

AWGM_PG_Parameters::AWGM_PG_Parameters(double _tol, double _alpha, size_t _max_its,
		size_t _max_basissize, bool _reset_res, bool _print_info,
		bool _verbose, bool _plot_solution, bool _verbose_extra,
		size_t _hashmapsize_trial, size_t _hashmapsize_test,
		std::string _info_filename,	std::string _plot_filename)
: tol(_tol), alpha(_alpha), max_its(_max_its), max_basissize(_max_basissize),
  reset_res(_reset_res), print_info(_print_info),
  verbose(_verbose), plot_solution(_plot_solution), verbose_extra(_verbose_extra),
  hashmapsize_trial(_hashmapsize_trial), hashmapsize_test(_hashmapsize_test),
  info_filename(_info_filename), plot_filename(_plot_filename)
{}

void
AWGM_PG_Parameters::print()
{
	std::cout << "###### AWGM Parameters #################" << std::endl;
	std::cout << std::left << std::setw(24) << "# tol:" << std::setw(20) <<  tol << std::endl;
	std::cout << std::left << std::setw(24) << "# alpha:" << std::setw(20) <<  alpha << std::endl;
	std::cout << std::left << std::setw(24) << "# max_its:" << std::setw(20) <<  max_its << std::endl;
	std::cout << std::left << std::setw(24) << "# max_basissize:" << std::setw(20) <<  max_basissize << std::endl;
	std::cout << std::left << std::setw(24) << "# reset_res:" << std::setw(20) <<  (reset_res?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# print_info:" << std::setw(20) <<  (print_info?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << (verbose_extra?" (extra)":"") << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_solution:" << std::setw(20) <<  (plot_solution?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# hashmapsize trial:" << std::setw(20) <<  hashmapsize_trial << std::endl;
	std::cout << std::left << std::setw(24) << "# hashmapsize test:" << std::setw(20) <<  hashmapsize_test << std::endl;
	std::cout << std::left << std::setw(24) << "# info_filename:" << std::setw(20) <<  info_filename << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_filename:" << std::setw(20) <<  plot_filename << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}


IS_Parameters::IS_Parameters(bool _adaptive_tol, size_t _max_its, double _init_tol,
			  double _res_reduction, double _absolute_tol, bool _verbose)
: adaptive_tol(_adaptive_tol), max_its(_max_its), init_tol(_init_tol),
  res_reduction(_res_reduction), absolute_tol(_absolute_tol), verbose(_verbose)
{}

void
IS_Parameters::print()
{
	std::cout << "###### Inner Solver Parameters ##########" << std::endl;
	std::cout << std::left << std::setw(18) << "# adaptive_tol:" << std::setw(16) <<  (adaptive_tol?"true":"false") << std::endl;
	std::cout << std::left << std::setw(18) << "# max_its:" << std::setw(16) <<  max_its << std::endl;
	std::cout << std::left << std::setw(18) << "# init_tol:" << std::setw(16) <<  init_tol << std::endl;
	std::cout << std::left << std::setw(18) << "# res_reduction:" << std::setw(16) <<  res_reduction << std::endl;
	std::cout << std::left << std::setw(18) << "# absolute_tol:" << std::setw(16) <<  absolute_tol << std::endl;
	std::cout << std::left << std::setw(18) << "# verbose:" << std::setw(16) <<  (verbose?"true":"false") << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}

AWGM_Parameters::AWGM_Parameters(double _tol, double _alpha, size_t _max_its,
		size_t _max_basissize, bool _print_info,
		bool _verbose, bool _plot_solution, bool _verbose_extra,
		size_t _hashmapsize, std::string _info_filename, std::string _plot_filename)
: tol(_tol), alpha(_alpha), max_its(_max_its), max_basissize(_max_basissize),
  print_info(_print_info), verbose(_verbose), plot_solution(_plot_solution),
  verbose_extra(_verbose_extra), hashmapsize(_hashmapsize),
  info_filename(_info_filename), plot_filename(_plot_filename)
{}

void
AWGM_Parameters::print()
{
	std::cout << "###### AWGM Parameters #################" << std::endl;
	std::cout << std::left << std::setw(24) << "# tol:" << std::setw(20) <<  tol << std::endl;
	std::cout << std::left << std::setw(24) << "# alpha:" << std::setw(20) <<  alpha << std::endl;
	std::cout << std::left << std::setw(24) << "# max_its:" << std::setw(20) <<  max_its << std::endl;
	std::cout << std::left << std::setw(24) << "# max_basissize:" << std::setw(20) <<  max_basissize << std::endl;
	std::cout << std::left << std::setw(24) << "# print_info:" << std::setw(20) <<  (print_info?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << (verbose_extra?" (extra)":"") << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_solution:" << std::setw(20) <<  (plot_solution?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# hashmapsize:" << std::setw(20) <<  hashmapsize << std::endl;
	std::cout << std::left << std::setw(24) << "# info_filename:" << std::setw(20) <<  info_filename << std::endl;
	std::cout << std::left << std::setw(24) << "# plot_filename:" << std::setw(20) <<  plot_filename << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}


void
AWGM_PG_Information::print(const char* filename)
{
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# It Res Res_NE SizeTrial SizeTest SizeTestResNE SizeTestRes CGLS_Its" << std::endl;
    	for(size_t i=0; i < awgm_res.size(); ++i){
    		infofile << i << " " << awgm_res[i] << " " << awgm_resNE[i] << " "
    				<< sizeLambdaTrial[i] << " " << sizeLambdaTest[i] << " "
    				<< sizeLambdaResNE[i] << " " << sizeLambdaRes[i] << " "
    				<< cgls_its[i] << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}

void
AWGM_Information::print(const char* filename)
{
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# It Res SizeLambda SizeTestRes CG_Its" << std::endl;
    	for(size_t i=0; i < awgm_res.size(); ++i){
    		infofile << i << " " << awgm_res[i] << " " << sizeLambda[i] << " "
    		        << sizeLambdaRes[i] << " " << cg_its[i] << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}


} // namespace lawa
