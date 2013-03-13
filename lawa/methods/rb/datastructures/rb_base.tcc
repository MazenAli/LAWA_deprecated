namespace lawa {

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
RB_Base<RB_Model,TruthModel, DataType, ParamType>::RB_Base(RB_Model& _rb_system, TruthModel& _rb_truth)
 : rb_system(_rb_system), rb_truth(_rb_truth)
{}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
void
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
train_Greedy()
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

	// In order to be able to calculate empty error bounds,
	// we have to calculate the Riesz Representors for F

	if(greedy_params.verbose){
	    std::cout << "||---------------------------------------------------------------------||" << std::endl;
	    std::cout << "||-------------- Calculate Riesz Representors for F -------------------||" << std::endl;
	    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
	}

	// Calculate the Riesz Representors for F
	size_t Qf = rb_system.Q_f();
    for (unsigned int i = 0; i < Qf; ++i) {

    	if(greedy_params.verbose){
    	    std::cout << "||------- Component Nb " << std::setw(3) << i << "  -------------------------------------------||" << std::endl << std::endl;
    	}

        DataType c = rb_truth.get_riesz_representor_f(i);
        F_representors.push_back(c);
    }

	// Update the Riesz Representor Norms
    rb_system.F_F_representor_norms.engine().resize((int)rb_system.Q_f(), (int)rb_system.Q_f());
    DataType vec1, vec2;
    for(size_t qf1 = 1; qf1 <= Qf; ++qf1) {
        for (size_t qf2 = qf1; qf2 <= Qf; ++qf2) {
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


	ParamType current_param;
	T max_error = 0;
	size_t N = 0;
	do {

		if(greedy_params.verbose){
		    std::cout << "||---------------------------------------------------------------------||" << std::endl;
		    std::cout << "||-------------- Greedy Search for new parameter  ---------------------||" << std::endl;
		    std::cout << "||---------------------------------------------------------------------||" << std::endl << std::endl;
		}

		T error_est;
		for(auto& mu : Xi_train){
            typename RB_Model::DenseVectorT u_N = rb_system.get_rb_solution(N, mu);
            error_est = rb_system.get_errorbound(u_N,mu);

            if(greedy_params.verbose){
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

		if(greedy_params.verbose){
			std::cout << std::endl << "Greedy Error = " << std::scientific << max_error << std::endl << std::endl;
		}

		if(greedy_params.verbose){
		    std::cout << "||=====================================================================||" << std::endl;
		    std::cout << "||       SNAPSHOT  at ";  ParamInfo<ParamType>::print(current_param); std::cout << std::endl;
		    std::cout << "||=====================================================================||" << std::endl << std::endl;
		}

		DataType u = rb_truth.get_truth_solution(current_param);
		add_to_basis(u);
		N++;

	} while( (N < greedy_params.Nmax) && (max_error > greedy_params.tol) );
}

template <typename RB_Model, typename TruthModel, typename DataType, typename ParamType>
std::vector<ParamType>
RB_Base<RB_Model,TruthModel, DataType, ParamType>::
generate_uniform_paramset(ParamType min_param, ParamType max_param, intArrayType param_nb)
{
	size_t pdim = ParamInfo<ParamType>::dim;

	// Calculate step lengths in each parameter dimension
    std::vector<T> h(pdim);
    for(size_t d = 0; d < pdim; ++d) {
        h[d] = (max_param[d] - min_param[d]) / ((T)param_nb[d]-1.);
    }

    std::vector<ParamType> Xi_train;

    /*
    std::vector<size_t>	   index(pdim,0);

    while(true){
        ParamType new_mu;
        for(size_t i = 0; i < pdim; ++i){
        	new_mu[i] = std::min(min_param[index[i]] + i*h[index[i]], max_param[index[i]]);
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
    */


    // Generate Parameter Grid
    if(pdim == 1) {
        for (size_t i = 0; i < greedy_params.nb_training_params[0]; ++i){
            ParamType new_mu;
            new_mu[0]  = std::min(greedy_params.min_param[0] + i*h[0], greedy_params.max_param[0]);
            Xi_train.push_back(new_mu);
        }
    }
    else{
        if (pdim == 2) {
            for (size_t i = 0; i < greedy_params.nb_training_params[0]; ++i){
                for (size_t j = 0; j < greedy_params.nb_training_params[1]; ++j){
                	ParamType new_mu;
                    new_mu[0]  = std::min(greedy_params.min_param[0] + i*h[0], greedy_params.max_param[0]);
                    new_mu[1]  = std::min(greedy_params.min_param[1] + j*h[0], greedy_params.max_param[1]);
                    Xi_train.push_back(new_mu);
                }
            }
        }
        else {
            std::cerr << "Generate Trainingsset for dim = " << pdim << " : Not implemented yet " << std::endl;
        }
    }

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
	for(auto& bf : rb_basisfunctions){
		new_bf = new_bf - bf * rb_truth.innprod_Y_u_u(bf, u);
	}

	T new_bf_norm_sq = rb_truth.innprod_Y_u_u(new_bf, new_bf);
	new_bf.scale(1./std::sqrt(new_bf_norm_sq));
	rb_basisfunctions.push_back(new_bf);

	if(greedy_params.verbose){
	    std::cout << "||------- GRAM-SCHMIDT Norm: " << std::setw(10) << std::sqrt(new_bf_norm_sq) << "------------||" << std::endl;
	}

	exit(0);
}

} // namespace lawa
