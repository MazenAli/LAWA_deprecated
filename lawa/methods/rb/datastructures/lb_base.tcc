#include <iostream>

namespace lawa {

template <typename ParamType, typename Truth>
LB_Base<ParamType,Truth>::LB_Base(Truth& _truth, ThetaStructure<ParamType>& _thetas_lhs)
 : truth(_truth), thetas(_thetas_lhs)
{}

template <typename ParamType, typename Truth>
template <typename IndexSetType>
void
LB_Base<ParamType,Truth>::
assemble_matrices_for_alpha_computation(IndexSetType& indexset)
{
	A_u_u_matrices.clear();
	A_u_u_matrices.resize(thetas.size());
	for(auto& matrix : A_u_u_matrices){
		matrix.resize((int)indexset.size(), (int)indexset.size());
	}
	I_Y_u_u_matrix.resize((int)indexset.size(), (int)indexset.size());

    int col_count = 1;
    for (auto& ind_col : indexset) {
    	int row_count = 1;
        for (auto& ind_row : indexset) {

        	// Get entry in A_u_u_matrices
        	for(std::size_t i = 0; i < thetas.size(); ++i){
        		T tmp = truth.lhs_u_u(i,ind_row, ind_col);
        		if(tmp > 1e-12){
        			A_u_u_matrices[i](row_count, col_count) = tmp;
        		}
        	}
        	// Get entry for I_u_u_matrix
        	T tmp = truth.innprod_Y_u_u(ind_row, ind_col);
        	if(tmp > 1e-12){
        		I_Y_u_u_matrix(row_count,col_count) = tmp;
        	}

            row_count++;
        }
        col_count++;
    }

	for(std::size_t i = 0; i < A_u_u_matrices.size(); ++i){
		A_u_u_matrices[i].finalize();
	}
    I_Y_u_u_matrix.finalize();

}

template <typename ParamType, typename Truth>
typename Truth::T
LB_Base<ParamType,Truth>::
calculate_alpha(ParamType& mu)
{
	assert(I_Y_u_u_matrix.numRows() > 0);
	assert(A_u_u_matrices.size() > 0);
	assert(A_u_u_matrices[0].numRows() == I_Y_u_u_matrix.numRows());

	std::size_t N = I_Y_u_u_matrix.numRows();

	SparseMatrixT A(N,N);

	std::cout << "    Calculating alpha: " << std::endl;
	std::cout << "      Assembling Sparse Matrix .... " << std::endl;

	for(std::size_t i = 0; i < thetas.size(); ++i){
		for(auto it = A_u_u_matrices[i].begin(); it != A_u_u_matrices[i].end(); ++it){
			// Find Transposed entry
			auto it_transp = A_u_u_matrices[i].begin();
			for(; it_transp != A_u_u_matrices[i].end(); ++it_transp){
				if( (*it_transp).first.first == (*it).first.second && (*it_transp).first.second == (*it).first.first){
					break;
				}
			}

			if(it_transp != A_u_u_matrices[i].end()){
				A((*it).first.first, (*it).first.second) += thetas.eval(i, mu) * ((*it).second + (*it_transp).second);

			}
			else{
				A((*it).first.first, (*it).first.second) += thetas.eval(i, mu) * (*it).second;
			}

		}
	}
	A.finalize();
	A /= 2.;

	FullColMatrixT A_dense, I_dense;
	std::cout << "      Densify A .... " << std::endl;
	densify(cxxblas::NoTrans, A, A_dense);
	std::cout << "      Densify I .... " << std::endl;
	densify(cxxblas::NoTrans, I_Y_u_u_matrix, I_dense);

	std::cout << "      Solving EVP .... " << std::endl;
    DenseVectorT    wr(N), wi(N), beta(N);
    FullColMatrixT  vl, vr;
    gv(false, false, A_dense, I_dense, wr, wi, beta, vl, vr);
    T eps = 10e-8;
    std::vector<T> evals;
    for(std::size_t n = 1; n <= N; n++){
        if(wi(n) < eps){
            evals.push_back(wr(n)/beta(n));
        }
    }
    std::sort(evals.begin(), evals.end());

	return evals[0];
}


} // namespace lawa
