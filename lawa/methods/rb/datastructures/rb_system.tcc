#include <cfloat>
#include <cassert>
#include <iomanip>
#include <sys/stat.h>

namespace lawa {

template<typename T, typename ParamType>
RB_System<T,ParamType>::RB_System(ThetaStructure<ParamType>& _thetas_a,
							 ThetaStructure<ParamType>& _thetas_f)
 : thetas_a(_thetas_a), thetas_f(_thetas_f)
{}

template<typename T, typename ParamType>
std::size_t
RB_System<T,ParamType>::Q_f()
{
	return thetas_f.size();
}

template<typename T, typename ParamType>
std::size_t
RB_System<T,ParamType>::Q_a()
{
	return thetas_a.size();
}

template<typename T, typename ParamType>
flens::DenseVector<flens::Array<T> >
RB_System<T,ParamType>::
get_rb_solution(std::size_t N, ParamType& mu)
{
	if(N==0){
		return DenseVectorT(0);
	}

	assert(RB_F_vectors.size() > 0);
	assert(RB_F_vectors[0].length() >= (int)N);
	assert(RB_A_matrices[0].numRows() >= (int)N);

	FullColMatrixT A(N, N);

	for (std::size_t i = 0; i < thetas_a.size(); ++i) {
		//A += thetas_a.eval(i,mu) * RB_A_matrices[i](_(1,N), _(1,N));
		flens::axpy(cxxblas::NoTrans, thetas_a.eval(i,mu), RB_A_matrices[i](_(1,N), _(1,N)), A);
	}

	DenseVectorT F(N);
	for (unsigned int i = 0; i < thetas_f.size(); ++i) {
		//F += thetas_f.eval(i,mu) * RB_F_vectors[i](_(1,N));
		flens::axpy(thetas_f.eval(i,mu), RB_F_vectors[i](_(1,N)), F);
	}

	DenseVectorT u(N);

	int its;
	switch (rb_params.call) {
		case call_cg:
			its = cg(A, u, F);
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " cg iterations" << std::endl;
			}
			break;
		case call_gmres:
			std::cout << "A = " << A << std::endl;
			std::cout << "F = " << F << std::endl;
			its = gmres(A, u, F);
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " gmres iterations" << std::endl;
			}
			break;
		default:
			if(rb_params.verbose){
				std::cerr << "RB solution: unrecognized solver call " << std::endl;
			}
			exit(1);
			break;
	}

	return u;
}

template<typename T, typename ParamType>
flens::DenseVector<flens::Array<T> >
RB_System<T,ParamType>::
get_rb_solution(std::vector<std::size_t> indices, ParamType& mu)
{
	std::size_t N = indices.size();
	if(N==0){
		return DenseVectorT(0);
	}

	// First assemble all up to largest index...
	std::size_t Nmax = *(std::max_element(indices.begin(), indices.end()))+1;

	assert(RB_F_vectors.size() > 0);
	assert(RB_F_vectors[0].length() >= (int)Nmax);
	assert(RB_A_matrices[0].numRows() >= (int)Nmax);

	FullColMatrixT Afull(Nmax, Nmax);

	for (std::size_t i = 0; i < thetas_a.size(); ++i) {
		//A += thetas_a.eval(i,mu) * RB_A_matrices[i](_(1,N), _(1,N));
		flens::axpy(cxxblas::NoTrans, thetas_a.eval(i,mu), RB_A_matrices[i](_(1,Nmax), _(1,Nmax)), Afull);
	}

	DenseVectorT Ffull(Nmax);
	for (std::size_t i = 0; i < thetas_f.size(); ++i) {
		//F += thetas_f.eval(i,mu) * RB_F_vectors[i](_(1,N));
		flens::axpy(thetas_f.eval(i,mu), RB_F_vectors[i](_(1,Nmax)), Ffull);
	}

	// Then throw away unnecessary info
	FullColMatrixT  A(N,N);
	DenseVectorT	F(N);

	std::size_t count_i=1;
	for(auto& i : indices){
		F(count_i) = Ffull(i+1);
		std::size_t count_j=1;
		for(auto& j : indices){
			A(count_i, count_j) = Afull(i+1,j+1);
			count_j++;
		}
		count_i++;
	}

	// Then use the rest
	DenseVectorT u(N);

	int its;
	switch (rb_params.call) {
		case call_cg:
			its = cg(A, u, F);
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " cg iterations" << std::endl;
			}
			break;
		case call_gmres:
//			std::cout << "A = " << A << std::endl;
//			std::cout << "F = " << F << std::endl;
			its = gmres(A, u, F);
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " gmres iterations" << std::endl;
			}
			break;
		default:
			if(rb_params.verbose){
				std::cerr << "RB solution: unrecognized solver call " << std::endl;
			}
			exit(1);
			break;
	}

	return u;
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
get_errorbound(const DenseVectorT& u_N, ParamType& mu)
{
    return  residual_dual_norm(u_N, mu) / alpha_LB(mu);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
get_errorbound_accuracy(const DenseVectorT& u_N, ParamType& mu, T eps_f, T eps_a){

	T sum_theta_f = 0;
	T sum_theta_a = 0;
	T sum_u_coeffs = 0;
    for (std::size_t i = 1; i <= thetas_f.size(); ++i) {
    	sum_theta_f += thetas_f.eval(i-1,mu);
    }
    for (std::size_t i = 1; i <= thetas_a.size(); ++i) {
    	sum_theta_a += thetas_a.eval(i-1,mu);
    }
    for(int i = 1; i <= u_N.length(); ++i){
    	sum_u_coeffs += u_N(i);
    }

    std::cout << "Sum_Thetas_F = " << sum_theta_f << ", Sum_Thetas_A = " << sum_theta_a << ", sum_u_coeffs = " << sum_u_coeffs << std::endl;

    return std::sqrt(eps_f*eps_f*sum_theta_f*sum_theta_f
    		+ 2*eps_f*eps_a*sum_theta_f*sum_u_coeffs*sum_theta_a
    		+ eps_a*eps_a*sum_u_coeffs*sum_u_coeffs*sum_theta_a*sum_theta_a);
}


template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
residual_dual_norm(const DenseVectorT& u_N, ParamType& mu)
{

    T res_dual_norm = 0;

    int N = u_N.length();

    std::size_t Qf = thetas_f.size();
    std::size_t Qa = thetas_a.size();

	assert(F_F_representor_norms.numRows() >= (int)Qf);
	assert(A_F_representor_norms.size() >= (std::size_t)N);
	assert(A_A_representor_norms.size() >= (std::size_t)N);


    DenseVectorT ThetaF(Qf);
    DenseVectorT ThetaA(Qa);
    for (std::size_t i = 1; i <= Qf; ++i) {
        ThetaF(i) = thetas_f.eval(i-1,mu);
    }
    for (std::size_t i = 1; i <= Qa; ++i) {
        ThetaA(i) = thetas_a.eval(i-1,mu);
    }

    DenseVectorT FF_T = F_F_representor_norms * ThetaF;

    res_dual_norm = ThetaF * FF_T;

	if(N==0){
		return std::sqrt(res_dual_norm);
	}

    DenseVectorT T_AF_T(N);
    FullColMatrixT T_AA_T(N,N);
    for (int n1 = 1; n1 <= N; ++n1) {
        DenseVectorT AF_T = A_F_representor_norms[n1-1] * ThetaF;
        T_AF_T(n1) = ThetaA * AF_T;
        for(int n2 = n1; n2 <= N; ++n2) {
            DenseVectorT AA_T = A_A_representor_norms[n1-1][n2-n1] * ThetaA;
            DenseVectorT T_AA = transpose(A_A_representor_norms[n1-1][n2-n1]) * ThetaA;
            T_AA_T(n1, n2) = ThetaA * AA_T;
            T_AA_T(n2, n1) = ThetaA * T_AA;
        }
    }

    //std::cout << " Residual Dual Norm: size(u) = " << u_RB.length() << ", size(T_AF_T) = " << T_AF_T.length() << std::endl;
    //res_dual_norm += 2 * u_RB * T_AF_T;

    DenseVectorT T_AA_T_u = T_AA_T * u_N;
    res_dual_norm += u_N * T_AA_T_u;
    res_dual_norm += 2. * u_N * T_AF_T;

    if(res_dual_norm < 0){
      std::cout << "Warning: Residual dual norm negative: " << std::setprecision(20) << res_dual_norm << std::endl;
      res_dual_norm = std::fabs(res_dual_norm);
    }

	std::cout << " res_dual_norm = " << std::sqrt(res_dual_norm) << std::endl;

    return std::sqrt(res_dual_norm);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
residual_dual_norm(std::vector<std::size_t> indices, const DenseVectorT& u_N, ParamType& mu)
{

    T res_dual_norm = 0;

	std::size_t N = indices.size();
	std::size_t Nmax = *(std::max_element(indices.begin(), indices.end()))+1;

    std::size_t Qf = thetas_f.size();
    std::size_t Qa = thetas_a.size();

	assert(F_F_representor_norms.numRows() >= (int)Qf);
	assert(A_F_representor_norms.size() >= (std::size_t)Nmax);
	assert(A_A_representor_norms.size() >= (std::size_t)Nmax);


    DenseVectorT ThetaF(Qf);
    DenseVectorT ThetaA(Qa);
    for (std::size_t i = 1; i <= Qf; ++i) {
        ThetaF(i) = thetas_f.eval(i-1,mu);
    }
    for (std::size_t i = 1; i <= Qa; ++i) {
        ThetaA(i) = thetas_a.eval(i-1,mu);
    }

    DenseVectorT FF_T = F_F_representor_norms * ThetaF;

    res_dual_norm = ThetaF * FF_T;

	if(N==0){
		return std::sqrt(res_dual_norm);
	}

    DenseVectorT T_AF_T(N);
    FullColMatrixT T_AA_T(N,N);
    std::vector<std::size_t> indices_1based;
    for(auto& i : indices){
    	indices_1based.push_back(i+1);
    }

    std::size_t count_n1=1;
    for (auto& n1 : indices_1based) {
        DenseVectorT AF_T = A_F_representor_norms[n1-1] * ThetaF;
        T_AF_T(count_n1) = ThetaA * AF_T;
        //for(int n2 = n1; n2 <= N; ++n2) {
        std::size_t count_n2 = 1;
        for(auto& n2 : indices_1based){
        	if(n2 < n1){
        		count_n2++;
        		continue;
        	}

            DenseVectorT AA_T = A_A_representor_norms[n1-1][n2-n1] * ThetaA;
            DenseVectorT T_AA = transpose(A_A_representor_norms[n1-1][n2-n1]) * ThetaA;
            T_AA_T(count_n1, count_n2) = ThetaA * AA_T;
            T_AA_T(count_n2, count_n1) = ThetaA * T_AA;
            count_n2++;
        }
        count_n1++;
    }

    //std::cout << " Residual Dual Norm: size(u) = " << u_RB.length() << ", size(T_AF_T) = " << T_AF_T.length() << std::endl;
    //res_dual_norm += 2 * u_RB * T_AF_T;

    DenseVectorT T_AA_T_u = T_AA_T * u_N;
    res_dual_norm += u_N * T_AA_T_u;
    res_dual_norm += 2. * u_N * T_AF_T;

    if(res_dual_norm < 0){
      std::cout << "Warning: Residual dual norm negative: " << std::setprecision(20) << res_dual_norm << std::endl;
      res_dual_norm = std::fabs(res_dual_norm);
    }

	//std::cout << " res_dual_norm = " << std::sqrt(res_dual_norm) << std::endl;

    return std::sqrt(res_dual_norm);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
alpha_LB(ParamType& mu)
{
    T alpha_lb = DBL_MAX;
    for(std::size_t qa = 0; qa < thetas_a.size(); ++qa){
        T reftheta = thetas_a.eval(qa, rb_params.ref_param);
        T mutheta = thetas_a.eval(qa, mu);
        if ((mutheta / reftheta) < alpha_lb) {
            alpha_lb = mutheta / reftheta;
        }
    }

	std::cout << " alpha_LB = " << alpha_lb << std::endl;
    return alpha_lb;
}

template<typename T, typename ParamType>
void
RB_System<T,ParamType>::
write_rb_data(const std::string& directory_name){

	const unsigned int highprecision = 40;

	// Make a directory to store all the data files
	if(mkdir(directory_name.c_str(), 0777) == -1)
	{
		if(rb_params.verbose){
			  std::cerr << "         [In RB_System::write_rb_data: Directory "
					    << directory_name << " already exists, overwriting contents.]" << std::endl;
		}
	}

	std::size_t n_bf = RB_inner_product.numRows();
	std::string n_bf_filename = directory_name + "/n_bf.txt";
	std::ofstream n_bf_file(n_bf_filename.c_str());
	n_bf_file << n_bf << std::endl;
	n_bf_file.close();

	// Write RB_A_matrices
	for(std::size_t i = 0; i < Q_a(); ++i){
		std::stringstream filename;
		filename << directory_name << "/RB_A_" << i+1 << ".txt";
		std::ofstream file(filename.str().c_str());
		file.precision(highprecision);
		file << std::scientific << RB_A_matrices[i] << std::endl;
		file.close();
	}

	// Write RB_F_vectors
	for(std::size_t i = 0; i < Q_f(); ++i){
		std::stringstream filename;
		filename << directory_name << "/RB_F_" << i+1 << ".txt";
		std::ofstream file(filename.str().c_str());
		file.precision(highprecision);
		file << std::scientific << RB_F_vectors[i] << std::endl;
		file.close();
	}

	// Write RB_inner_product
	std::stringstream filename;
	filename << directory_name << "/RB_inner_product.txt";
	std::ofstream file(filename.str().c_str());
	file.precision(highprecision);
	file << std::scientific << RB_inner_product << std::endl;
	file.close();

	// Write F_F_representor_norms
	std::stringstream repr_F_filename;
	repr_F_filename << directory_name << "/F_F_representor_norms.txt";
	std::ofstream repr_F_file(repr_F_filename.str().c_str());
	repr_F_file.precision(highprecision);
	repr_F_file << std::scientific << F_F_representor_norms << std::endl;
	repr_F_file.close();

	// Write A_F_representor_norms
	std::stringstream repr_A_F_filename;
	repr_A_F_filename << directory_name << "/A_F_representor_norms.txt";
	std::ofstream repr_A_F_file(repr_A_F_filename.str().c_str());
	repr_A_F_file.precision(highprecision);
	for(std::size_t i = 0; i < n_bf; ++i){
		repr_A_F_file << std::scientific << A_F_representor_norms[i] << std::endl;
	}
	repr_A_F_file.close();

 	// Write A_A_representor_norms
	std::stringstream repr_A_A_filename;
	repr_A_A_filename << directory_name << "/A_A_representor_norms.txt";
	std::ofstream repr_A_A_file(repr_A_A_filename.str().c_str());
	repr_A_A_file.precision(highprecision);
	for(std::size_t i = 0; i < n_bf; ++i){
		for(std::size_t j = i; j < n_bf; ++j){
			repr_A_A_file << std::scientific << A_A_representor_norms[i][j-i] << std::endl;
		}
	}
	repr_A_A_file.close();

}

template<typename T, typename ParamType>
void
RB_System<T,ParamType>::
read_rb_data(const std::string& directory_name, int nb)
{
	// Read Nr of BasisFunctions
	unsigned int n_bf;
	std::string n_bf_filename = directory_name + "/n_bf.txt";

	std::ifstream n_bf_file(n_bf_filename.c_str());
	if(n_bf_file.is_open()){
		n_bf_file >> n_bf;
	    n_bf_file.close();

	    if(rb_params.verbose){
		    std::cout << "Number of Basis Functions: " << n_bf << std::endl;
	    }
	}
	else{
	    std::cerr << "Unable to read number of basis functions: " << strerror(errno) << std::endl;
	    exit(1);
	}

	// If N < 0 or > n_bf, set it to n_bf
	nb = (nb < 0) ? n_bf : nb;
	nb = (nb > (int)n_bf) ? n_bf : nb;

	// Read RB_A_matrices
	RB_A_matrices.clear();
	for(std::size_t i = 0; i < Q_a(); ++i){
	    std::stringstream filename;
	    filename << directory_name << "/RB_A_" << i+1 << ".txt";
	    std::ifstream file(filename.str().c_str());
	    if(file.is_open()){
	    	FullColMatrixT RB_A(n_bf, n_bf);
	    	for(unsigned int n1 = 1; n1 <= n_bf; n1++){
	    		for(unsigned int n2 = 1; n2 <= n_bf; n2++){
	    			file >> RB_A(n1, n2);
	    		}
	    	}
	    	file.close();
	    	RB_A_matrices.push_back(RB_A(_(1,nb), _(1,nb)));
		    if(rb_params.verbose){
		    	std::cout << " Read " << filename.str() << std::endl;
		    }
	    }
	    else{
	    	std::cerr << "Unable to open file " << filename.str() << " for reading: "<< strerror(errno) << std::endl;
	    	exit(1);
	    }
	}

	// Read RB_F_vectors
	RB_F_vectors.clear();
	for(unsigned int i = 0; i < Q_f(); ++i){
		std::stringstream filename;
	    filename << directory_name << "/RB_F_" << i+1 << ".txt";
	    std::ifstream file(filename.str().c_str());
	    if(file.is_open()){
	    	DenseVectorT RB_F(n_bf);
	    	for(unsigned int n = 1; n <= n_bf; n++){
	    		file >> RB_F(n);
	    	}
	    	file.close();
	    	RB_F_vectors.push_back(RB_F(_(1,nb)));
		    if(rb_params.verbose){
		    	std::cout << " Read " << filename.str() << std::endl;
		    }
	    }
	    else{
	    	std::cerr << "Unable to open file " << filename.str() << " for reading: " << strerror(errno) << std::endl;
	    	exit(1);
	    }
	 }

	 // Read RB_inner_product
	 RB_inner_product.engine().resize((int)n_bf, (int)n_bf);
	 std::stringstream inner_product_filename;
	 inner_product_filename << directory_name << "/RB_inner_product.txt";
	 std::ifstream inner_product_file(inner_product_filename.str().c_str());
	 if(inner_product_file.is_open()){
		 for(unsigned int n1 = 1; n1 <= n_bf; n1++){
			 for(unsigned int n2 = 1; n2 <= n_bf; n2++){
				 inner_product_file >> RB_inner_product(n1, n2);
			 }
		 }
		 if(nb < (int)n_bf){
			 FullColMatrixT tmp(RB_inner_product);
			 RB_inner_product.engine().resize(nb, nb);
			 RB_inner_product = tmp(_(1,nb), _(1,nb));
		 }
		 inner_product_file.close();
		 if(rb_params.verbose){
			 std::cout << " Read " << inner_product_filename.str() << std::endl;
		 }
	 }
	 else{
		 std::cerr << "Unable to open file " << inner_product_filename.str() << " for reading: " << strerror(errno) << std::endl;
		 exit(1);
	 }

	 // Read F_F_representor norms
	 F_F_representor_norms.engine().resize((int)Q_f(), (int)Q_f());
	 std::stringstream F_F_filename;
	 F_F_filename << directory_name << "/F_F_representor_norms.txt";
	 std::ifstream F_F_file(F_F_filename.str().c_str());
	 if(F_F_file.is_open()){
		 for(unsigned int i = 1; i <= Q_f(); ++i){
			 for(unsigned int j = 1; j <= Q_f(); ++j){
				 F_F_file >> F_F_representor_norms(i,j);
			 }
		 }
		 F_F_file.close();
		 if(rb_params.verbose){
			 std::cout << " Read " << F_F_filename.str() << std::endl;
		 }
	 }
	 else{
		 std::cerr << "Unable to open file " << F_F_filename.str() << " for reading: " << strerror(errno) << std::endl;
		 exit(1);
	 }

	  /*// Read output_output_representor norms
	  output_output_representor_norms.engine().resize((int)Q_output(), (int)Q_output());
	  std::stringstream output_output_filename;
	  output_output_filename << directory_name << "/output_output_representor_norms.dat";
	  std::ifstream output_output_file(output_output_filename.str().c_str());
	  if(output_output_file.is_open()){
	    for(unsigned int i = 1; i <= Q_output(); ++i){
	      for(unsigned int j = 1; j <= Q_output(); ++j){
	        output_output_file >> output_output_representor_norms(i,j);
	      }
	    }
	    output_output_file.close();
	    std::cout << " Read " << output_output_filename.str() << std::endl;
	  }
	  else{
	    std::cerr << "Unable to open file " << output_output_filename.str() << " for reading!" << std::endl;
	    exit(1);
	  }*/

	 // Read A_F_representor norms
	 A_F_representor_norms.clear();
	 std::stringstream A_F_filename;
	 A_F_filename << directory_name << "/A_F_representor_norms.txt";
	 std::ifstream A_F_file(A_F_filename.str().c_str());
	 if(A_F_file.is_open()){
		 for(int n = 1; n <= nb; ++n){
			 FullColMatrixT A_F_matrix(Q_a(), Q_f());
			 for(unsigned int i = 1; i <= Q_a(); ++i){
				 for(unsigned int j = 1; j <= Q_f(); ++j){
					 A_F_file >> A_F_matrix(i,j);
				 }
			 }
			 A_F_representor_norms.push_back(A_F_matrix);
		 }
		 A_F_file.close();
		 if(rb_params.verbose){
			 std::cout << " Read " << A_F_filename.str() << std::endl;
		 }
	 }
	 else{
		 std::cerr << "Unable to open file " << A_F_filename.str() << " for reading: " << strerror(errno) << std::endl;
		 exit(1);
	 }

	 // Read A_A_representor norms
	 A_A_representor_norms.clear();
	 std::stringstream A_A_filename;
	 A_A_filename << directory_name << "/A_A_representor_norms.txt";
	 std::ifstream A_A_file(A_A_filename.str().c_str());
	 if(A_A_file.is_open()){
		 for(int n1 = 1; n1 <= nb; ++n1){
			 std::vector<FullColMatrixT> A_A_n_vec;
			 for(unsigned int n2 = n1; n2 <= n_bf; ++n2){	// We have to read every information...
				 FullColMatrixT A_A_matrix(Q_a(), Q_a());
				 for(unsigned int i = 1; i <= Q_a(); ++i){
					 for(unsigned int j = 1; j <= Q_a(); ++j){
						 A_A_file >> A_A_matrix(i,j);
					 }
				 }
				 if((int)n2 <= nb){				// .. but safe only for n <= nb
					 A_A_n_vec.push_back(A_A_matrix);
				 }
			 }
			 A_A_representor_norms.push_back(A_A_n_vec);
		 }
		 A_A_file.close();
		 if(rb_params.verbose){
			 std::cout << " Read " << A_A_filename.str() << std::endl;
		 }
	 }
	 else{
		 std::cerr << "Unable to open file " << A_A_filename.str() << " for reading: " << strerror(errno) << std::endl;
		 exit(1);
	 }
}

} // namespace lawa
