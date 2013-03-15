#include <iostream>
#include <fstream>
#include <array>

namespace lawa {

template<typename T,std::size_t N>
void
ParamInfo<std::array<T,N> >::print(std::array<T,N> mu)
{
	std::cout << "[";
	for(auto& el : mu){
		std::cout << " " << std::fixed << std::right << std::setw(18) << el;
	}
	std::cout << "] ";
}

template<typename ParamType>
RB_Greedy_Parameters<ParamType>::RB_Greedy_Parameters(double _tol, size_t _Nmax,
		ParamType _min_param, ParamType _max_param,
		intArray _training_params_per_dim, bool _print_info,
		std::string _print_file, bool _verbose,
		bool _write_during_training, std::string _trainingdata_folder,
		bool _print_paramset, bool _erase_snapshot_params)
 : tol(_tol), Nmax(_Nmax), min_param(_min_param), max_param(_max_param),
   nb_training_params(_training_params_per_dim), print_info(_print_info),
   print_file(_print_file), verbose(_verbose),
   write_during_training(_write_during_training), trainingdata_folder(_trainingdata_folder),
   print_paramset(_print_paramset),
   erase_snapshot_params(_erase_snapshot_params)
{}

template<typename ParamType>
void
RB_Greedy_Parameters<ParamType>::print()
{
	std::cout << "###### RB Training Parameters #################" << std::endl;
	std::cout << std::left << std::setw(24) << "# tol:" << std::setw(20) <<  tol << std::endl;
	std::cout << std::left << std::setw(24) << "# Nmax:" << std::setw(20) <<  Nmax << std::endl;
	std::cout << std::left << std::setw(24) << "# Min_Param:" << std::setw(2) << "[ ";
		for(auto& el : min_param){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(24) << "# Max_Param:" << std::setw(2) << "[ ";
		for(auto& el : max_param){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(24) << "# Nb of train. params:" << std::setw(2) << "[ ";
		for(auto& el : nb_training_params){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(24) << "# print_info:" << std::setw(20) <<  (print_info?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# print in folder:" << std::setw(20) <<  trainingdata_folder << std::endl;
	std::cout << std::left << std::setw(24) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# write during training:" << std::setw(20) <<  (write_during_training?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# print_paramset:" << std::setw(20) <<  (print_paramset?"true":"false") << std::endl;
	std::cout << std::left << std::setw(24) << "# erase snapshot params:" << std::setw(20) <<  (erase_snapshot_params?"true":"false") << std::endl;
	std::cout << "###############################################" << std::endl << std::endl;

}

template<typename ParamType>
void
RB_Greedy_Information<ParamType>::RB_Greedy_Information::print(const char* filename)
{
	size_t pdim = ParamInfo<ParamType>::dim;
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# N Error Mu Size(u) Size(ReprF_qf) ... Size(Repr_A_bf_qa)" << std::endl;
    	for(size_t i=0; i < greedy_errors.size(); ++i){
    		// Write Error
    		infofile << i << " " << greedy_errors[i] << " ";
    		// Write Parameter
    		for(size_t d = 0; d < pdim; ++d){
    			infofile << snapshot_params[i][d] << " ";
    		}
    		// Write Size of Solution
    		infofile << u_size[i] << " ";
    		// Write Size of Riesz Reprs of F
    		for(auto& el : repr_f_size){
    			infofile << el << " ";
    		}
    		// Write Size of Riesz Reprs of A
			for(auto& el : repr_a_size[i]){
				infofile << el << " ";
			}
    		infofile << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}

template<typename ParamType>
RB_Parameters<ParamType>::RB_Parameters(SolverCall _call, ParamType _ref_param, bool _verbose)
 : call(_call), ref_param(_ref_param), verbose(_verbose)
{}

template<typename ParamType>
void
RB_Parameters<ParamType>::print()
{
	std::cout << "###### RB Parameters #################" << std::endl;
	std::cout << std::left << std::setw(24) << "# SolverType:" << std::setw(20);
		switch(call){
			case call_cg: 	std::cout << "cg" << std::endl;
							break;
			case call_gmres:std::cout << "gmres" << std::endl;
							break;
			default:		std::cout << "<unrecognized>" << std::endl;
							break;
		}
	std::cout << std::left << std::setw(24) << "# Ref_Param:" << std::setw(2) << "[ ";
		for(auto& el : ref_param){
			std::cout << std::right << std::setw(4) << el << " ";
		}
		std::cout << std::left << " ]" << std::endl;
	std::cout << std::left << std::setw(24) << "# verbose:" << std::setw(20) <<  (verbose?"true":"false") << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}

} // namespace lawa
