#include <algorithm>
#include <iomanip>
#include <sys/stat.h>

namespace lawa {

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
RB_Base<RB_Model,TruthModel, DataType, ParamType>::RB_Base(RB_Model& _rb_system, TruthModel& _rb_truth)
 : rb_system(_rb_system), rb_truth(_rb_truth)
{}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
std::size_t
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
n_bf()
{
	return rb_basisfunctions.size();
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
DataType&
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
get_bf(std::size_t i)
{
	assert(i < rb_basisfunctions.size());
	return rb_basisfunctions[i];
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
train_Greedy(std::size_t N)
{
	if(greedy_params.verbose){
		std::cout.precision(12);
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=========      OFFLINE TRAINING  " << (greedy_params.training_type==weak?"(WEAK)":"(STRONG)") << "            ================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl << std::endl;

	    std::cout << "Starting with N = " << N << std::endl << std::endl;

	}

	std::vector<ParamType> Xi_train = generate_uniform_paramset(greedy_params.min_param,
																greedy_params.max_param,
																greedy_params.nb_training_params);
	if(greedy_params.print_paramset){
		std::cout << "Training Parameters: " << std::endl;
		print_paramset(Xi_train);
	}


	// Preparations for strong Greedy (Output Folder / Calculation of Truth Sols)
	std::map<ParamType, DataType> truth_sols;
	if(greedy_params.training_type == strong){
		std::string truthdatafolder;
		std::string paraminfo_filename;
		std::ofstream paraminfo_file;
		if(greedy_params.write_during_training){
			truthdatafolder = greedy_params.trainingdata_folder + "/truthsolutions";

			// Make a directory to store all the basisfunction files
			if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
			{
				if(rb_system.rb_params.verbose){
					  std::cerr << "         [In RB_Base::train_strong_Greedy: Directory "
							    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
				}
			}
			if(mkdir(truthdatafolder.c_str(), 0777) == -1)
			{
				if(rb_system.rb_params.verbose){
					  std::cerr << "         [In RB_Base::train_strong_Greedy: Directory "
							    << truthdatafolder << " already exists, overwriting contents.]" << std::endl;
				}
			}

			// Open file for parameter information
			paraminfo_filename = truthdatafolder + "/paraminfo.txt";
		}

		// Calculate all truth solutions
		std::size_t count = 0;
	    for (auto& mu : Xi_train) {
			DataType u = rb_truth.get_truth_solution(mu);
	    	truth_sols.insert(std::make_pair(mu, u));
	    	count++;

	    	if(greedy_params.write_during_training){
	    		std::stringstream filename;
	    		filename << truthdatafolder << "/truthsol_" << count << ".txt";
			    saveCoeffVector2D(u, rb_truth.get_trialbasis(), filename.str().c_str());

			    paraminfo_file.open(paraminfo_filename.c_str(), std::fstream::app);
			    paraminfo_file << count;
			    for(auto& d : mu){
			    	paraminfo_file << " " << d;
			    }
			    paraminfo_file << std::endl;
			    paraminfo_file.close();
	    	}
	    }
	}

	if(N == 0){
		// In order to be able to calculate empty error bounds,
		// we have to calculate the Riesz Representors for F
		calculate_Riesz_RHS_information();
	}

	ParamType current_param;
	T max_error = 0;
	do {

		if(greedy_params.verbose){
		    std::cout << "||---------------------------------------------------------------------||" << std::endl;
		    std::cout << "||-------------- Greedy Search for new parameter  ---------------------||" << std::endl;
		    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}


		// Find parameter with maximum error estimator
		T error_est;
		max_error = 0;
		for(auto& mu : Xi_train){
            typename RB_Model::DenseVectorT u_N = rb_system.get_rb_solution(N, mu);

            if(greedy_params.training_type == weak){
                error_est = rb_system.get_errorbound(u_N,mu);
            }
            else{
                DataType diff;
                if(N > 0){
                	diff = reconstruct_u_N(u_N,N);
                }
                diff -= (*truth_sols.find(mu)).second;
                error_est = std::sqrt(rb_truth.innprod_Y_u_u(diff,diff));
            }

            if(greedy_params.verbose){
            	std::cout << "    u_N = " << u_N;
            	std::cout << "    Mu = ";
            	ParamInfo<ParamType>::print(mu);
            	std::cout << " : Error = " << std::scientific << std::setw(15) << error_est << std::endl;
            }

            if(error_est > max_error){
            	max_error = error_est;
            	current_param = mu;
            }
		}

		// Check if we have to tighten the snapshot tolerance
		bool update_snapshot = false;
		if(greedy_params.tighten_tol){
			// find (first occurence of) current parameter
			for(auto& mu : greedy_info.snapshot_params){
				bool is_current_param = true;
				for(std::size_t i = 0; i < mu.size(); ++i){
					if(mu[i]!=current_param[i]){
						is_current_param = false;
						break;
					}
				}
				if(is_current_param){
					rb_truth.access_solver().access_params().tol *= greedy_params.tighten_tol_reduction;
					update_snapshot = true;
					if(greedy_params.verbose){
						std::cout << std::endl<< "====> -------------------------------------------------------------- <=====" << std::endl;
						std::cout << "====>     Reduced snapshot tolerance to " << std::scientific << rb_truth.access_solver().access_params().tol << std::endl;
						std::cout << "====> -------------------------------------------------------------- <=====" << std::endl;

					}

					if(greedy_params.tighten_tol_rieszA){
						rb_truth.access_RieszSolver_A().access_params().tol *= greedy_params.tighten_tol_reduction;
						if(greedy_params.verbose){
							std::cout << std::endl<< "====> -------------------------------------------------------------- <=====" << std::endl;
							std::cout << "====>     Reduced Riesz Repr A tolerance to " << std::scientific << rb_truth.access_RieszSolver_A().access_params().tol << std::endl;
							std::cout << "====> -------------------------------------------------------------- <=====" << std::endl;
						}
					}
					if(greedy_params.tighten_tol_rieszF){
						rb_truth.access_RieszSolver_F().access_params().tol *= greedy_params.tighten_tol_reduction;
						if(greedy_params.verbose){
							std::cout << std::endl<< "====> -------------------------------------------------------------- <=====" << std::endl;
							std::cout << "====>     Reduced Riesz Repr F tolerance to " << std::scientific << rb_truth.access_RieszSolver_F().access_params().tol << std::endl;
							std::cout << "====> -------------------------------------------------------------- <=====" << std::endl;
						}
						// Save old Riesz Representor
						if(greedy_params.write_during_training){
							if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
							{
								if(greedy_params.verbose){
									  std::cerr << "         [In RB_Base::train_Greedy: Directory "
											    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
								}
							}
							std::string repr_folder = greedy_params.trainingdata_folder + "/representors";
							if(mkdir(repr_folder.c_str(), 0777) == -1)
							{
								if(greedy_params.verbose){
									  std::cerr << "         [In RB_Base::train_Greedy: Directory "
											    << repr_folder << " already exists, overwriting contents.]" << std::endl;
								}
							}
							std::stringstream old_repr_folder;
							old_repr_folder << greedy_params.trainingdata_folder << "/representors/F_repr_before_it_" << N+1;
							write_rieszrepresentors(old_repr_folder.str(), 0);
						}
						// Recalculate F Representors
						bool update = false;
						if(greedy_params.update_rieszF){
							update = true;
						}
						calculate_Riesz_RHS_information(update);
						recalculate_A_F_norms();
					}
					break;
				}
			}
		}

		greedy_info.greedy_errors.push_back(max_error);
		greedy_info.snapshot_params.push_back(current_param);

		// Remove parameter from indexset
		if(greedy_params.erase_snapshot_params){
			// find current parameter
			for(auto& mu : Xi_train){
				bool is_current_param = true;
				for(std::size_t i = 0; i < mu.size(); ++i){
					if(mu[i]!=current_param[i]){
						is_current_param = false;
						break;
					}
				}
				if(is_current_param){
					auto mu_it = std::find(Xi_train.begin(), Xi_train.end(), mu);
					Xi_train.erase(mu_it);
					break;
				}
			}
		}



		if(greedy_params.verbose){
			std::cout << std::endl << "Greedy Error = " << std::scientific << max_error << std::endl << std::endl;
		}

		if(greedy_params.verbose){
		    std::cout << "||=====================================================================||" << std::endl;
		    std::cout << "||       SNAPSHOT  at ";  ParamInfo<ParamType>::print(current_param); std::cout << std::endl;
		    std::cout << "||=====================================================================||" << std::endl << std::endl;
		}

		if(rb_truth.access_solver().access_params().print_info){
			std::stringstream filename;
			filename << greedy_params.print_file << "_bf" << N+1 << ".txt";
			rb_truth.access_solver().access_params().info_filename = filename.str();
		}

		// Compute next snapshot
		DataType u;
		if(update_snapshot && greedy_params.update_snapshot){
			// Find last corresponding snapshot:
			// Find (last occurence of) current parameter
			int last_snapshot_index;
			for(auto& mu : greedy_info.snapshot_params){
				bool is_current_param = true;
				for(std::size_t i = 0; i < mu.size(); ++i){
					if(mu[i]!=current_param[i]){
						is_current_param = false;
						break;
					}
				}
				if(is_current_param){
					if(&mu - &greedy_info.snapshot_params[0] < greedy_info.snapshot_params.size()-1){
						last_snapshot_index = &mu - &greedy_info.snapshot_params[0];
					}
				}
			}

			if(greedy_params.verbose){
				std::cout << " Starting truth computation with snapshot nb " << last_snapshot_index+1 << std::endl << std::endl;
			}

			// Start solver with copy of old snapshot
			u = rb_basisfunctions[last_snapshot_index];
			rb_truth.get_truth_solution(current_param, u);

			if(greedy_params.training_type == strong){
				truth_sols.erase(current_param);
		    	truth_sols.insert(std::make_pair(current_param, u));
			}
		}
		else{
			if(greedy_params.training_type == strong){
				for(auto& entry : truth_sols){
					bool is_current_mu = true;
					for(std::size_t i = 0; i < current_param.size(); ++i){
						if(entry.first[i]!=current_param[i]){
							is_current_mu = false;
							break;
						}
					}
					if(is_current_mu){
						u = entry.second;
					}
				}
			}
			else{
				u = rb_truth.get_truth_solution(current_param);
			}
		}
		add_to_basis(u);
		N++;

		greedy_info.u_size.push_back(u.size());
		std::vector<std::size_t> f_sizes;
		for(auto& el : F_representors){
			f_sizes.push_back(el.size());
		}
		greedy_info.repr_f_size.push_back(f_sizes);
		std::vector<std::size_t> a_sizes;
		for(auto& el : A_representors[n_bf()-1]){
			a_sizes.push_back(el.size());
		}
		greedy_info.repr_a_size.push_back(a_sizes);

		if(greedy_params.write_during_training){

			if(greedy_params.verbose){
				std::cout << "=====>  Writing RB Greedy Training information to file " << std::endl << std::endl;
			}
			if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
			{
				if(greedy_params.verbose){
					  std::cerr << "         [In RB_Base::train_Greedy: Directory "
							    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
				}
			}

			// Write Basis Functions
			std::string bf_folder = greedy_params.trainingdata_folder + "/bf";
			write_basisfunctions(bf_folder, (int)N);

			// Write Riesz Representors
			std::string repr_folder = greedy_params.trainingdata_folder + "/representors";
			write_rieszrepresentors(repr_folder);

			// Write RB Data
			rb_system.write_rb_data(greedy_params.trainingdata_folder);

			// Write Training Information
			greedy_info.print(greedy_params.print_file.c_str());
		}

	} while( (N < greedy_params.Nmax) && (max_error > greedy_params.tol) );

	if(greedy_params.print_info){
		greedy_info.print();
	}
}

/*template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
train_strong_Greedy()
{
	if(greedy_params.verbose){
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=========      OFFLINE TRAINING                      ================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl;
	    std::cout << "||=====================================================================||" << std::endl << std::endl;

	}

	std::vector<ParamType> Xi_train = generate_uniform_paramset(greedy_params.min_param,
																greedy_params.max_param,
																greedy_params.nb_training_params);
	if(greedy_params.print_paramset){
		std::cout << "Training Parameters: " << std::endl;
		print_paramset(Xi_train);
	}

	std::string truthdatafolder;
	std::string paraminfo_filename;
	std::ofstream paraminfo_file;
	if(greedy_params.write_during_training){
		truthdatafolder = greedy_params.trainingdata_folder + "/truthsolutions";

		// Make a directory to store all the basisfunction files
		if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
		{
			if(rb_system.rb_params.verbose){
				  std::cerr << "         [In RB_Base::train_strong_Greedy: Directory "
						    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
			}
		}
		if(mkdir(truthdatafolder.c_str(), 0777) == -1)
		{
			if(rb_system.rb_params.verbose){
				  std::cerr << "         [In RB_Base::train_strong_Greedy: Directory "
						    << truthdatafolder << " already exists, overwriting contents.]" << std::endl;
			}
		}

		// Open file for parameter information
		paraminfo_filename = truthdatafolder + "/paraminfo.txt";
	}

	// Calculate all truth solutions
	std::map<ParamType, DataType> truth_sols;
	std::size_t count = 0;
    for (auto& mu : Xi_train) {
		DataType u = rb_truth.get_truth_solution(mu);
    	truth_sols.insert(std::make_pair(mu, u));
    	count++;

    	if(greedy_params.write_during_training){
    		std::stringstream filename;
    		filename << truthdatafolder << "/truthsol_" << count << ".txt";
		    saveCoeffVector2D(u, rb_truth.get_trialbasis(), filename.str().c_str());

		    paraminfo_file.open(paraminfo_filename.c_str(), std::fstream::app);
		    paraminfo_file << count;
		    for(auto& d : mu){
		    	paraminfo_file << " " << d;
		    }
		    paraminfo_file << std::endl;
		    paraminfo_file.close();
    	}
    }

	// In order to be able to calculate empty error bounds,
	// we have to calculate the Riesz Representors for F
	calculate_Riesz_RHS_information();

	for(auto& el : F_representors){
		greedy_info.repr_f_size.push_back(el.size());
	}

	ParamType current_param;
	std::size_t N = 0;
	T max_error = 0;
	do {

		if(greedy_params.verbose){
		    std::cout << "||---------------------------------------------------------------------||" << std::endl;
		    std::cout << "||-------------- Greedy Search for new parameter  ---------------------||" << std::endl;
		    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}

		T error_est;
		max_error = 0;
		for(auto& mu : Xi_train){
            typename RB_Model::DenseVectorT u_N = rb_system.get_rb_solution(N, mu);

            //error_est = rb_system.get_errorbound(u_N,mu);
            DataType diff;
            if(N > 0){
            	diff = reconstruct_u_N(u_N,N);
            }
            diff -= (*truth_sols.find(mu)).second;
            error_est = std::sqrt(rb_truth.innprod_Y_u_u(diff,diff));

            if(greedy_params.verbose){
            	std::cout << "    u_N = " << u_N;
            	std::cout << "    Mu = ";
            	ParamInfo<ParamType>::print(mu);
            	std::cout << " : Error = " << std::scientific << std::setw(15) << error_est << std::endl;
            }

            if(error_est > max_error){
            	max_error = error_est;
            	current_param = mu;
            }
		}

		greedy_info.greedy_errors.push_back(max_error);
		greedy_info.snapshot_params.push_back(current_param);

		// Remove parameter from indexset
		if(greedy_params.erase_snapshot_params){
			for(auto& mu : Xi_train){
				bool is_current_param = true;
				for(std::size_t i = 0; i < mu.size(); ++i){
					if(mu[i]!=current_param[i]){
						is_current_param = false;
						break;
					}
				}
				if(is_current_param){
					auto mu_it = std::find(Xi_train.begin(), Xi_train.end(), mu);
					Xi_train.erase(mu_it);
					break;
				}
			}
		}


		if(greedy_params.verbose){
			std::cout << std::endl << "Greedy Error = " << std::scientific << max_error << std::endl << std::endl;
		}

		if(greedy_params.verbose){
		    std::cout << "||=====================================================================||" << std::endl;
		    std::cout << "||       SNAPSHOT  at ";  ParamInfo<ParamType>::print(current_param); std::cout << std::endl;
		    std::cout << "||=====================================================================||" << std::endl << std::endl;
		}

		if(rb_truth.access_solver().access_params().print_info){
			std::stringstream filename;
			filename << greedy_params.print_file << "_bf" << N+1 << ".txt";
			rb_truth.access_solver().access_params().info_filename = filename.str();
		}

		//DataType u = rb_truth.get_truth_solution(current_param);
		DataType u;
		for(auto& entry : truth_sols){
			bool is_current_mu = true;
			for(std::size_t i = 0; i < current_param.size(); ++i){
				if(entry.first[i]!=current_param[i]){
					is_current_mu = false;
					break;
				}
			}
			if(is_current_mu){
				u = entry.second;
			}
		}
		add_to_basis(u);
		N++;

		greedy_info.u_size.push_back(u.size());
		std::vector<std::size_t> a_sizes;
		for(auto& el : A_representors[n_bf()-1]){
			a_sizes.push_back(el.size());
		}
		greedy_info.repr_a_size.push_back(a_sizes);

		if(greedy_params.write_during_training){

			if(greedy_params.verbose){
				std::cout << "=====>  Writing RB Greedy Training information to file " << std::endl << std::endl;
			}
			if(mkdir(greedy_params.trainingdata_folder.c_str(), 0777) == -1)
			{
				if(greedy_params.verbose){
					  std::cerr << "         [In RB_Base::train_Greedy: Directory "
							    << greedy_params.trainingdata_folder << " already exists, overwriting contents.]" << std::endl;
				}
			}

			// Write Basis Functions
			std::string bf_folder = greedy_params.trainingdata_folder + "/bf";
			write_basisfunctions(bf_folder, (int)N);

			// Write Riesz Representors
			std::string repr_folder = greedy_params.trainingdata_folder + "/representors";
			write_rieszrepresentors(repr_folder);
			if(N == 1){
				write_rieszrepresentors(repr_folder);
			}

			// Write RB Data
			rb_system.write_rb_data(greedy_params.trainingdata_folder);

			// Write Training Information
			greedy_info.print(greedy_params.print_file.c_str());
		}

	} while( (N < greedy_params.Nmax) && (max_error > greedy_params.tol) );

	if(greedy_params.print_info){
		greedy_info.print();
	}
}*/

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
DataType
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
reconstruct_u_N(typename RB_Model::DenseVectorT u, std::size_t N)
{
	assert(N <= n_bf());
	assert(u.length() > 0);
	assert(u.length() >= (int)N);

	DataType u_full;
	for (unsigned int i = 1; i <= N; ++i) {
		u_full +=  u(i) * rb_basisfunctions[i-1];
	}

	return u_full;
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
std::vector<ParamType>
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
generate_uniform_paramset(ParamType min_param, ParamType max_param, intArrayType param_nb)
{
	std::size_t pdim = ParamInfo<ParamType>::dim;

	// Calculate step lengths in each parameter dimension
    std::vector<T> h(pdim);
    for(std::size_t d = 0; d < pdim; ++d) {
    	if(param_nb[d] > 1){
            h[d] = (max_param[d] - min_param[d]) / ((T)param_nb[d]-1.);
    	}
    	else{
    		h[d] = 0;
    	}
    }

    std::vector<ParamType> Xi_train;

    std::vector<std::size_t>	   index(pdim,0);

    while(true){
        ParamType new_mu;
        for(std::size_t i = 0; i < pdim; ++i){
        	new_mu[i] = std::min(min_param[i] + index[i]*h[i], max_param[i]);
        }
        Xi_train.push_back(new_mu);

        for(int i = pdim-1; ; --i){
        	if(i<0){
        		return Xi_train;
        	}
        	index[i]++;
        	if(index[i]==param_nb[i]){
        		index[i] = 0;
        	}
        	else{
        		break;
        	}
        }
    }


    /*
    // Generate Parameter Grid
    if(pdim == 1) {
        for (std::size_t i = 0; i < greedy_params.nb_training_params[0]; ++i){
            ParamType new_mu;
            new_mu[0]  = std::min(greedy_params.min_param[0] + i*h[0], greedy_params.max_param[0]);
            Xi_train.push_back(new_mu);
        }
    }
    else{
        if (pdim == 2) {
            for (std::size_t i = 0; i < greedy_params.nb_training_params[0]; ++i){
                for (std::size_t j = 0; j < greedy_params.nb_training_params[1]; ++j){
                	ParamType new_mu;
                    new_mu[0]  = std::min(greedy_params.min_param[0] + i*h[0], greedy_params.max_param[0]);
                    new_mu[1]  = std::min(greedy_params.min_param[1] + j*h[1], greedy_params.max_param[1]);
                    Xi_train.push_back(new_mu);
                }
            }
        }
        else {
            std::cerr << "Generate Trainingsset for dim = " << pdim << " : Not implemented yet " << std::endl;
        }
    }
    */

    return Xi_train;
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
print_paramset(std::vector<ParamType> PSet)
{
	for(auto& el : PSet){
		ParamInfo<ParamType>::print(el);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
add_to_basis(const DataType& u)
{
	if(greedy_params.verbose){
	    std::cout << "||---------------------------------------------------------------------||" << std::endl;
	    std::cout << "||-------------- Adding Snapshot to Basis   ---------------------------||" << std::endl;
	    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
	}

	// =========== Orthogonalization ======================= //

	DataType new_bf = u;

	if(greedy_params.orthonormalize_bfs){
		for(auto& bf : rb_basisfunctions){
			new_bf = new_bf - bf * rb_truth.innprod_Y_u_u(bf, u);
		}
	}

	T new_bf_norm_sq = rb_truth.innprod_Y_u_u(new_bf, new_bf);
	new_bf.scale(1./std::sqrt(new_bf_norm_sq));


	rb_basisfunctions.push_back(new_bf);

	if(greedy_params.verbose){
	    std::cout << "||------- GRAM-SCHMIDT Norm: " << std::setw(10) << std::sqrt(new_bf_norm_sq) << " --------------------------||" << std::endl;

	    std::cout << std::endl << "||------- UPDATE RB Structures ----------------------------------------||" << std::endl;
	}


	add_to_RB_structures(new_bf);

	update_Riesz_LHS_information(new_bf);

}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
add_to_RB_structures(const DataType& bf)
{
	std::size_t Qa = rb_system.Q_a();
	std::size_t Qf = rb_system.Q_f();

	// ===== Update RB_A ====== //

	// We need some matrix to start with...
	bool first_bf = false;
    if (rb_system.RB_A_matrices.size() < Qa) {
        for (std::size_t q_a = 0; q_a < Qa; ++q_a) {
        	typename RB_Model::FullColMatrixT A(1,1);
            rb_system.RB_A_matrices.push_back(A);
        }
        first_bf = true;
    }

    for (std::size_t q_a = 0; q_a < Qa; ++q_a) {
    	// "Pad" RB_A_matrices with zeros to new size
    	std::size_t n;
		if (first_bf == true) {
			rb_system.RB_A_matrices[q_a].engine().resize(1, 1);
			n = 1;
		}
		else {
			typename RB_Model::FullColMatrixT tmp(rb_system.RB_A_matrices[q_a]);
			rb_system.RB_A_matrices[q_a].engine().resize(tmp.numRows()+1, tmp.numCols()+1);
			rb_system.RB_A_matrices[q_a](tmp.rows(), tmp.cols()) = tmp;

			n = rb_system.RB_A_matrices[q_a].numRows();
			for(unsigned int i = 1; i <= n; ++i) {
				rb_system.RB_A_matrices[q_a](i, n) = 0.;
				rb_system.RB_A_matrices[q_a](n, i) = 0.;
			}

		}

		// Compute new entries (one column, one row)
        for (unsigned int i = 1; i <= n; ++i) {
        	rb_system.RB_A_matrices[q_a](n,i) = rb_truth.lhs_u_u(q_a, bf, rb_basisfunctions[i-1]);
        }
        for (unsigned int i = 1; i < n; ++i) {
        	rb_system.RB_A_matrices[q_a](i,n) = rb_truth.lhs_u_u(q_a, rb_basisfunctions[i-1], bf);
        }

        if(greedy_params.verbose){
    		std::cout << std::endl << "||------- RB_A (" << q_a << ")  -----------------------------------------||" << std::endl;
    		std::cout << rb_system.RB_A_matrices[q_a] << std::endl;
        }
    }

	// ===== Update RB_F ====== //

	// We need some matrix to start with...
    if (rb_system.RB_F_vectors.size() < Qf) {
    	rb_system.RB_F_vectors.resize(Qf);
    	first_bf = true;
    }
    for (std::size_t q_f = 0; q_f < Qf; ++q_f) {
    	// Add one zero entry to RB_F_vectors
    	std::size_t n;
    	if (first_bf == true) {
    		rb_system.RB_F_vectors[q_f].engine().resize(1);
    		n = 1;
    	}
    	else {
    		typename RB_Model::DenseVectorT tmp(rb_system.RB_F_vectors[q_f]);
    		rb_system.RB_F_vectors[q_f].engine().resize(tmp.length()+1);
    		n = rb_system.RB_F_vectors[q_f].length();
    		rb_system.RB_F_vectors[q_f](tmp.range()) = tmp;
    		rb_system.RB_F_vectors[q_f](n) = 0.;

    	}

    	rb_system.RB_F_vectors[q_f](n) = rb_truth.rhs_u(q_f, bf);

    	if(greedy_params.verbose){
    		std::cout << std::endl << "||------- RB_F (" << q_f << ")  -----------------------------------------||" << std::endl;
    		std::cout << rb_system.RB_F_vectors[q_f] << std::endl;
    	}
    }

	// ===== Update RB_InnerProduct ====== //

    std::size_t n;
    if (rb_system.RB_inner_product.numRows()==0) {
        rb_system.RB_inner_product.engine().resize(1, 1);
        n=1;
    }
    else {
    	typename RB_Model::FullColMatrixT tmp(rb_system.RB_inner_product);
    	rb_system.RB_inner_product.engine().resize(tmp.numRows()+1, tmp.numCols()+1);
    	rb_system. RB_inner_product(tmp.rows(), tmp.cols()) = tmp;
    	n = rb_system. RB_inner_product.numRows();
        for(unsigned int i = 1; i <= n; ++i) {
        	rb_system.RB_inner_product(i, n) = 0.;
        	rb_system.RB_inner_product(n, i) = 0.;
        }
    }

    for (unsigned int i = 1; i <= n; ++i) {
    	rb_system.RB_inner_product(n, i) = rb_truth.innprod_Y_u_u(bf, rb_basisfunctions[i-1]);
		rb_system.RB_inner_product(i, n) = rb_system.RB_inner_product(n,i);
    }

	if(greedy_params.verbose){
		std::cout << std::endl << "||------- RB_InnerProduct  ----------------------------------||" << std::endl;
		std::cout << rb_system.RB_inner_product << std::endl;
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
calculate_Riesz_RHS_information(bool update)
{

	if(!update){
		F_representors.clear();

		if(greedy_params.verbose){
			std::cout << "||---------------------------------------------------------------------||" << std::endl;
			std::cout << "||-------------- Calculate Riesz Representors for F -------------------||" << std::endl;
			std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}
	}
	else{

		if(greedy_params.verbose){
			std::cout << "||---------------------------------------------------------------------||" << std::endl;
			std::cout << "||-------------- Update Riesz Representors for F --------   -----------||" << std::endl;
			std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}
	}
	// Calculate the Riesz Representors for F
	std::size_t Qf = rb_system.Q_f();
	for (unsigned int i = 0; i < Qf; ++i) {

		if(greedy_params.verbose){
			std::cout << "||------- Component Nb " << std::setw(3) << i << "  -------------------------------------------||" << std::endl << std::endl;
		}

		if(update){
			rb_truth.get_riesz_representor_f(i, F_representors[i]);
		}
		else{
			DataType c = rb_truth.get_riesz_representor_f(i);
			F_representors.push_back(c);
		}
	}

	// Update the Riesz Representor Norms
	rb_system.F_F_representor_norms.engine().resize((int)rb_system.Q_f(), (int)rb_system.Q_f());
	DataType vec1, vec2;
	for(std::size_t qf1 = 1; qf1 <= Qf; ++qf1) {
		for (std::size_t qf2 = qf1; qf2 <= Qf; ++qf2) {
			vec1 = F_representors[qf1-1];
			vec2 = F_representors[qf2-1];

			rb_system.F_F_representor_norms(qf1,qf2) = rb_truth.innprod_Y_v_v(vec1, vec2);

			if(qf1 != qf2) {
				rb_system.F_F_representor_norms(qf2,qf1) = rb_system.F_F_representor_norms(qf1,qf2);
			}
		}
	}

	if(greedy_params.verbose){
		std::cout << std::endl << "||------- Representor Norms F x F  ------------------------------------||" << std::endl;
		std::cout << rb_system.F_F_representor_norms << std::endl;
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
update_Riesz_LHS_information(const DataType& bf)
{

	// Check if we only do an update and search for "old" bf index
	bool update = false;
	int last_snapshot_index = -1;
	if(greedy_params.update_rieszA){
		for(auto& mu : greedy_info.snapshot_params){
			bool is_last_param = true;
			for(std::size_t i = 0; i < mu.size(); ++i){
				if(mu[i]!=greedy_info.snapshot_params[greedy_info.snapshot_params.size()-1][i]){
					is_last_param = false;
					break;
				}
			}
			if(is_last_param){
				if(&mu - &greedy_info.snapshot_params[0] < greedy_info.snapshot_params.size()-1){
					last_snapshot_index = &mu - &greedy_info.snapshot_params[0];
				}
			}
		}
		if(last_snapshot_index >= 0){
			update=true;

			if(greedy_params.verbose){
				std::cout << "||---------------------------------------------------------------------||" << std::endl;
				std::cout << "||-------------- Update Riesz Representors for A -------------------||" << std::endl;
				std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;

				std::cout << " Starting representor computation with snapshot nb " << last_snapshot_index+1 << std::endl << std::endl;
			}
		}
		else{
			if(greedy_params.verbose){
				std::cout << "||---------------------------------------------------------------------||" << std::endl;
				std::cout << "||-------------- Calculate Riesz Representors for A -------------------||" << std::endl;
				std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
			}
		}
	}
	else{
		if(greedy_params.verbose){
			std::cout << "||---------------------------------------------------------------------||" << std::endl;
			std::cout << "||-------------- Calculate Riesz Representors for A -------------------||" << std::endl;
			std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}
	}


	// Calculate the Riesz Representors for A
	std::size_t Qa = rb_system.Q_a();
	std::vector<DataType> new_A_reprs(Qa);
	for (unsigned int i = 0; i < Qa; ++i) {

		if(greedy_params.verbose){
			std::cout << "||------- Component Nb " << std::setw(3) << i << "  -------------------------------------------||" << std::endl << std::endl;
		}

		DataType c;
		if(update && greedy_params.update_rieszA){
			// Start solver with copy of old Riesz representor
			c = A_representors[last_snapshot_index][i];
			std::cout << "Old Representor: " << c.size() << " indizes" << std::endl;
			rb_truth.get_riesz_representor_a(i, bf, c, greedy_params.coarsen_rieszA_for_update);
		}
		else{
			c = rb_truth.get_riesz_representor_a(i, bf);
		}
		new_A_reprs[i] = c;
	}
	A_representors.push_back(new_A_reprs);

	// Update the Riesz Representor Norms A x A
	int N = rb_basisfunctions.size();
	DataType vec1, vec2;
    for(int n1 = 0; n1 < N; ++n1) {
    	typename RB_Model::FullColMatrixT A_n1_N(Qa, Qa);
        for(std::size_t qa1 = 1; qa1 <= Qa; ++qa1) {
        	for(std::size_t qa2 = qa1; qa2 <= Qa; ++qa2) {
        		vec1 = A_representors[n1][qa1-1];
        		vec2 = A_representors[N-1][qa2-1];
        		A_n1_N(qa1, qa2) = rb_truth.innprod_Y_v_v(vec1, vec2);

        		if(qa1 != qa2){
        			if(n1 == N-1){
        				A_n1_N(qa2, qa1) = A_n1_N(qa1, qa2);
        			}
        			else{
        				vec1 = A_representors[n1][qa2-1];
        				vec2 = A_representors[N-1][qa1-1];
        				A_n1_N(qa2, qa1) = rb_truth.innprod_Y_v_v(vec1, vec2);
        			}
        		}
        	}
        }
        if(n1 == N-1){
        	std::vector<typename RB_Model::FullColMatrixT> newvec;
        	newvec.push_back(A_n1_N);
        	rb_system.A_A_representor_norms.push_back(newvec);
        }
        else{
        	rb_system.A_A_representor_norms[n1].push_back(A_n1_N);
        }

    	if(greedy_params.verbose){
    		std::cout << std::endl << "||------- Representor Norms A["<< n1 <<"] x A["<< N-1 << "]  -----------------------------||" << std::endl;
    		std::cout << rb_system.A_A_representor_norms[n1][N-1-n1] << std::endl;
    	}
    }

	// Update the Riesz Representor Norms A x F
    std::size_t Qf = rb_system.Q_f();
    typename RB_Model::FullColMatrixT A_F(Qa, Qf);
    for(std::size_t qa = 1; qa <= Qa; ++qa) {
        for(std::size_t qf = 1; qf <= Qf; ++qf) {
            vec1 = A_representors[N-1][qa-1];
            vec2 = F_representors[qf-1];
            A_F(qa, qf) = rb_truth.innprod_Y_v_v(vec1, vec2);
        }
    }
    rb_system.A_F_representor_norms.push_back(A_F);
	if(greedy_params.verbose){
		std::cout << std::endl << "||------- Representor Norms A["<< N-1 <<"] x F  --------------------------------||" << std::endl;
		std::cout << rb_system.A_F_representor_norms[N-1] << std::endl;
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
recalculate_A_F_norms(){

	rb_system.A_F_representor_norms.clear();

	for(std::size_t N = 1; N <= n_bf(); ++N){
		// Update the Riesz Representor Norms A x F
	    std::size_t Qa = rb_system.Q_a();
	    std::size_t Qf = rb_system.Q_f();
	    typename RB_Model::FullColMatrixT A_F(Qa, Qf);
	    for(std::size_t qa = 1; qa <= Qa; ++qa) {
	        for(std::size_t qf = 1; qf <= Qf; ++qf) {
	            DataType vec1 = A_representors[N-1][qa-1];
	            DataType vec2 = F_representors[qf-1];
	            A_F(qa, qf) = rb_truth.innprod_Y_v_v(vec1, vec2);
	        }
	    }
	    rb_system.A_F_representor_norms.push_back(A_F);
		if(greedy_params.verbose){
			std::cout << std::endl << "||------- Representor Norms A["<< N-1 <<"] x F  --------------------------------||" << std::endl;
			std::cout << rb_system.A_F_representor_norms[N-1] << std::endl;
		}
	}

}


template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
write_basisfunctions(const std::string& directory_name, int nb){

	if(rb_system.rb_params.verbose){
		std::cout << "=====>  Writing RB BasisFunctions to file " << std::endl << std::endl;
	}

	// Make a directory to store all the data files
	if(mkdir(directory_name.c_str(), 0777) == -1)
	{
		if(rb_system.rb_params.verbose){
			  std::cerr << "         [In RB_Base::write_basisfunctions: Directory "
					    << directory_name << " already exists, overwriting contents.]" << std::endl;
		}
	}

	std::string n_bf_filename = directory_name + "/n_bf.txt";
	std::ofstream n_bf_file(n_bf_filename.c_str());
	n_bf_file <<  rb_basisfunctions.size() << std::endl;
	n_bf_file.close();

	if(nb < 0){
		for(std::size_t i = 0; i < rb_basisfunctions.size(); ++i){
		    std::stringstream filename;
		    filename << directory_name << "/bf_" << i+1 << ".txt";
		    saveCoeffVector2D(rb_basisfunctions[i], rb_truth.get_trialbasis(), filename.str().c_str());
		}
	}
	else{
		assert((std::size_t)nb-1 < rb_basisfunctions.size());
	    std::stringstream filename;
	    filename << directory_name << "/bf_" << nb << ".txt";
	    saveCoeffVector2D(rb_basisfunctions[nb-1], rb_truth.get_trialbasis(), filename.str().c_str());
	}


}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
read_basisfunctions(const std::string& directory_name, int nb){

	unsigned int n_bf;
	if(nb < 0){
		std::string n_bf_filename = directory_name + "/n_bf.txt";

		std::ifstream n_bf_file(n_bf_filename.c_str());
		if(n_bf_file.is_open()){
			n_bf_file >> n_bf;
			n_bf_file.close();
		}
		else{
			std::cerr << "Unable to read number of basis functions: " << strerror(errno) << std::endl;
			exit(1);
		}
	}
	else{
		n_bf = nb;
	}


	rb_basisfunctions.clear();
	for(unsigned int i = 1; i <= n_bf; ++i){
		std::stringstream filename;
		filename << directory_name << "/bf_" << i << ".txt";

		DataType bf_coeffs;
		readCoeffVector2D(bf_coeffs, filename.str().c_str(),false);
		rb_basisfunctions.push_back(bf_coeffs);

		if(rb_system.rb_params.verbose){
			std::cout << " Read " << filename.str() << std::endl;
		}
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
write_rieszrepresentors(const std::string& directory_name, int nb)
{
	// Make a directory to store all the data files
	if(mkdir(directory_name.c_str(), 0777) == -1)
	{
		if(rb_system.rb_params.verbose){
			  std::cerr << "         [In RB_Base::write_rieszrepresentors: Directory "
					    << directory_name << " already exists, overwriting contents.]" << std::endl;
		}
	}

	if(nb <= 0){
		for(std::size_t i = 0; i < rb_system.Q_f(); ++i){
			std::stringstream filename;
			filename << directory_name << "/F_representor_" << i+1 << ".txt";
			saveCoeffVector2D(F_representors[i], rb_truth.get_testbasis(), filename.str().c_str());
		}
		if(nb < 0){
			for(std::size_t n = 0; n < rb_basisfunctions.size(); ++n){
				for(std::size_t i = 0; i < rb_system.Q_a(); ++i){
					std::stringstream filename;
					filename << directory_name << "/A_representor_" << i+1 << "_" << n+1 << ".txt";
					saveCoeffVector2D(A_representors[n][i], rb_truth.get_testbasis(), filename.str().c_str());
				}
			}
		}
	}
	else{
		assert((std::size_t)nb < rb_basisfunctions.size());

		for(std::size_t i = 0; i < rb_system.Q_a(); ++i){
			std::stringstream filename;
			filename << directory_name << "/A_representor_" << i+1 << "_" << nb+1 << ".txt";
			saveCoeffVector2D(A_representors[nb][i], rb_truth.get_testbasis(), filename.str().c_str());
		}
	}

}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
read_rieszrepresentors(const std::string& directory_name, int nb)
{

	F_representors.clear();
	for(std::size_t i = 0; i < rb_system.Q_f(); ++i){
		std::stringstream filename;
		filename << directory_name << "/F_representor_" << i+1 << ".txt";

		DataType rieszf_coeffs;
		readCoeffVector2D(rieszf_coeffs, filename.str().c_str(),false);
		F_representors.push_back(rieszf_coeffs);

		if(rb_system.rb_params.verbose){
			std::cout << " Read " << filename.str() << std::endl;
		}
	}

	A_representors.clear();
	std::size_t n_bf = (nb < 0) ? rb_basisfunctions.size() : nb;

	for(std::size_t n = 0; n < n_bf; ++n){
		std::vector<DataType> A_reprs_n;
		for(std::size_t i = 0; i < rb_system.Q_a(); ++i){
			std::stringstream filename;
			filename << directory_name << "/A_representor_" << i+1 << "_" << n+1 << ".txt";

			DataType riesza_coeffs;
			readCoeffVector2D(riesza_coeffs, filename.str().c_str(),false);
			A_reprs_n.push_back(riesza_coeffs);

			if(rb_system.rb_params.verbose){
				std::cout << " Read " << filename.str() << std::endl;
			}
		}
		A_representors.push_back(A_reprs_n);
	}
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
read_greedy_info(const char* filename, int nb){
	greedy_info.read(filename, rb_system.Q_f(), rb_system.Q_a(), nb);
}

} // namespace lawa
