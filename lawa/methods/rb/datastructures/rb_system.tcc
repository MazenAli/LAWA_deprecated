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
get_rb_solution(size_t N, ParamType& mu)
{
	if(N==0){
		return DenseVectorT(0);
	}

	assert(RB_F_vectors.size() > 0);
	assert(RB_F_vectors[0].length() >= (int)N);
	assert(RB_A_matrices[0].numRows() >= (int)N);

	FullColMatrixT A(N, N);

	for (size_t i = 0; i < thetas_a.size(); ++i) {
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
residual_dual_norm(const DenseVectorT& u_N, ParamType& mu)
{

    T res_dual_norm = 0;

    int N = u_N.length();

    size_t Qf = thetas_f.size();
    size_t Qa = thetas_a.size();

	assert(F_F_representor_norms.numRows() >= (int)Qf);
	assert(A_F_representor_norms.size() >= (size_t)N);
	assert(A_A_representor_norms.size() >= (size_t)N);


    DenseVectorT ThetaF(Qf);
    DenseVectorT ThetaA(Qa);
    for (size_t i = 1; i <= Qf; ++i) {
        ThetaF(i) = thetas_f.eval(i-1,mu);
    }
    for (size_t i = 1; i <= Qa; ++i) {
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

    return std::sqrt(res_dual_norm);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
alpha_LB(ParamType& mu)
{
    T alpha_lb = DBL_MAX;
    for(size_t qa = 0; qa < thetas_a.size(); ++qa){
        T reftheta = thetas_a.eval(qa, rb_params.ref_param);
        T mutheta = thetas_a.eval(qa, mu);
        if ((mutheta / reftheta) < alpha_lb) {
            alpha_lb = mutheta / reftheta;
        }
    }

    return alpha_lb;
}

} // namespace lawa
