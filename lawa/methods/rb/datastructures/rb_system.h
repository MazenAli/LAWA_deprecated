/*
 * rb_system.h
 *
 *  Created on: 07.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_

#include <flens/flens.h>
#include <array>
#include <vector>
#include <lawa/methods/rb/datastructures/thetastructure.h>
#include <lawa/methods/rb/datastructures/rb_parameters.h>


namespace lawa {

/* RB_System:
 * 		This class contains all N-dependent data and functions needed
 * 		for the computation of a RB approximation.
 */
template <typename T, typename ParamType>
class RB_System {

public:
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
    typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

	RB_System(ThetaStructure<ParamType>& _thetas_a,
			  ThetaStructure<ParamType>& _thetas_f);

	DenseVectorT
	get_rb_solution(std::size_t N, ParamType& mu);

	DenseVectorT
	get_rb_solution(std::vector<std::size_t> indices, ParamType& mu);

	virtual T
	get_errorbound(const DenseVectorT& u_N, ParamType& mu);

	T
	residual_dual_norm(const DenseVectorT& u_N, ParamType& mu);

	T
	residual_dual_norm(std::vector<std::size_t> indices, const DenseVectorT& u_N, ParamType& mu);

	virtual T
	alpha_LB(ParamType& mu);

	std::size_t
	Q_f();

	std::size_t
	Q_a();

	void
	write_rb_data(const std::string& directory_name = "offline_data");

	void
	read_rb_data(const std::string& directory_name = "offline_data", int nb = -1);

	RB_Parameters<ParamType> 					rb_params;

    ThetaStructure<ParamType>& 					thetas_a;
    ThetaStructure<ParamType>& 					thetas_f;

    std::vector<FullColMatrixT>     			RB_A_matrices;
    std::vector<DenseVectorT>       			RB_F_vectors;
    FullColMatrixT  							RB_inner_product;

    FullColMatrixT                              F_F_representor_norms; // Speicherbedarf kann verringert werden..
    std::vector<FullColMatrixT>                 A_F_representor_norms;
    std::vector<std::vector<FullColMatrixT> >   A_A_representor_norms; //.. Ausnutzen der Symmetrie (Matrix als Vektor)

};


} // namespace lawa

#include <lawa/methods/rb/datastructures/rb_system.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_ */
