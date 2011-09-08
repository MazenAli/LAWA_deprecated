#include <lawa/methods/rb/postprocessing/plotting.h>

namespace  lawa {

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::
AdaptiveRBTruth2D_PG(TrialBasis& _trialbasis, TestBasis& _testbasis, TrialPrec& _trialprec, TestPrec& _testprec,
                  bool _use_inner_product, bool _use_A_matrix)
    : trialbasis(_trialbasis), basis(_trialbasis), testbasis(_testbasis), lhs_op(this), rhs_op(this), repr_lhs_op(this), repr_rhs_A_op(this), repr_rhs_F_op(this),
      use_inner_product_matrix(_use_inner_product), use_A_operator_matrices(_use_A_matrix), 
      assembled_inner_product_matrix(false), assembled_prec_vec(false), assembled_A_operator_matrices(false),
      trial_prec(_trialprec), test_prec(_testprec), trial_prec_data(), test_prec_data()
{}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::
attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q)
{
    rb->theta_a.push_back(theta_a_q);
    A_operators.push_back(&A_q);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::attach_A_q(Operator2D<T>& A_q)
{
    A_operators.push_back(&A_q);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::
attach_F_q(theta_fctptr theta_f_q, AdaptiveRhs<T, Index2D>& F_q)
{
    rb->theta_f.push_back(theta_f_q);
    F_operators.push_back(&F_q);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::attach_F_q(AdaptiveRhs<T, Index2D>& F_q)
{
    F_operators.push_back(&F_q);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::
attach_inner_product_op(Operator2D<T>& _trial_inner_product_op, Operator2D<T>& _test_inner_product_op)
{
    trial_inner_product_op = &_trial_inner_product_op;
    test_inner_product_op = &_test_inner_product_op;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void 
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::set_truthsolver(TruthSolver& _truthsolver)
{
    solver = &_truthsolver;
    solver->set_model(*this);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::
set_rb_model(RBModel2D<T, AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression> >& _rb)
{
    rb = &_rb;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
RBModel2D<T, AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression> >&
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::get_rb_model()
{
    return *rb;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::truth_solve()
{
    CoeffVector u = solver->truth_solve();
    
    undo_trial_prec(u);
	    
    return u;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::undo_trial_prec(CoeffVector& u_trial)
{    
    // Undo preconditioning
    typename CoeffVector::iterator it;
    for(it = u_trial.begin(); it != u_trial.end(); ++it){
        (*it).second *= get_trial_prec((*it).first);
    }     
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::undo_test_prec(CoeffVector& u_test)
{
    // Undo preconditioning
    typename CoeffVector::iterator it;
    for(it = u_test.begin(); it != u_test.end(); ++it){
        (*it).second *= get_test_prec((*it).first);
    } 

}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::get_trial_prec(const Index2D& index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    
    if (!flens::IsSame<NoPreconditioner<T,Index2D>, TrialPrec>::value) {
        // Right precondioning:
        const_coeff_it it_index   = trial_prec_data.find(index);
        //  Entry has already been computed:
        if (it_index != trial_prec_data.end()) {
            prec *= (*it_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = trial_prec(index);
            trial_prec_data[index] = tmp;
            prec *= tmp;
        }
    }
    return prec;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::get_test_prec(const Index2D& index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    
    if (!flens::IsSame<NoPreconditioner<T,Index2D>, TestPrec>::value) {
        // Left precondioning:
        const_coeff_it it_index   = test_prec_data.find(index);
        //  Entry has already been computed:
        if (it_index != test_prec_data.end()) {
            prec *= (*it_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = test_prec(index);
            test_prec_data[index] = tmp;
            prec *= tmp;
        }
    }

    return prec;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::add_new_basis_function(const CoeffVector& sol)
{
    typename std::vector<CoeffVector>::iterator it;
    CoeffVector new_bf = sol;

    CoeffVector bf;
    for (it = rb->rb_basis_functions.begin(); it != rb->rb_basis_functions.end(); ++it) {
        bf = (*it);
        new_bf = new_bf - bf * trial_inner_product(bf, sol); 
    }
    
    std::cout << "Gram - Schmidt: norm^2 = " << trial_inner_product(new_bf, new_bf) << std::endl;
    new_bf.scale(1./std::sqrt(trial_inner_product(new_bf, new_bf)));
    rb->rb_basis_functions.push_back(new_bf);    
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::update_representors()
{
    int N = rb->n_bf();
    if(N == 1) {
        F_representors.resize(rb->Q_f());
        for (unsigned int i = 0; i < rb->Q_f(); ++i) {
            repr_rhs_F_op.set_current_op(*F_operators[i]);
            std::cout<< std::endl << "  --- Solving for Riesz Representor of F_" << i+1 << "---" << std::endl<< std::endl;
            CoeffVector c = solver->repr_solve_F();
            
            undo_test_prec(c);
            
            F_representors[i] = c;
        } 
    }
    
    A_representors.resize(N);
    repr_rhs_A_op.set_current_bf(rb->rb_basis_functions[N-1]);
    A_representors[N-1].resize(rb->Q_a());
    for (unsigned int i = 0; i < rb->Q_a(); ++i) {
        repr_rhs_A_op.set_current_op(*A_operators[i]);
        std::cout << std::endl<< "  --- Solving for Riesz Representor of A_" << i+1 
                  << "(" << N << ")" << " --- " << std::endl<< std::endl;
         
        CoeffVector c = solver->repr_solve_A();   
        
        undo_test_prec(c);
            
        A_representors[N-1][i] = c;
    }
 }

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::update_representor_norms()
{
  Timer timer;
    int N = rb->n_bf();
    int Qf = rb->Q_f();
    
    CoeffVector vec1, vec2;

    if(N == 1){
      std::cout << "  ... F x F .... " << std::endl;
       // F_F representor norms 
       rb->F_F_representor_norms.engine().resize(Qf, Qf);

       timer.start();
       for(int qf1 = 1; qf1 <= Qf; ++qf1) {
           for (int qf2 = qf1; qf2 <= Qf; ++qf2) {
               
               vec1 = F_representors[qf1-1];
               vec2 = F_representors[qf2-1];

               rb->F_F_representor_norms(qf1,qf2) = test_inner_product(vec1, vec2);

               if(qf1 != qf2) {
                   rb->F_F_representor_norms(qf2,qf1) = rb->F_F_representor_norms(qf1,qf2);
               }
           }
       }   
       timer.stop();
       std::cout << "  " << timer.elapsed() << " seconds" << std::endl << std::endl;
       std::cout << rb->F_F_representor_norms << std::endl;
    }
     
     
    std::cout << "  ...  A x A .... " << std::endl;    
     
     
    // A_A representor norms 
    int Qa = rb->Q_a();
    timer.start();
    for(int n1 = 0; n1 < N; ++n1) {
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_n1_N(Qa, Qa);
        for(int qa1 = 1; qa1 <= Qa; ++qa1) {
          for(int qa2 = qa1; qa2 <= Qa; ++qa2) {  
          
            vec1 = A_representors[n1][qa1-1];
            vec2 = A_representors[N-1][qa2-1];
                              
            A_n1_N(qa1, qa2) = test_inner_product(vec1, vec2);
            
            if(qa1 != qa2){
              if(n1 == N-1){
                A_n1_N(qa2, qa1) = A_n1_N(qa1, qa2);
              }
              else{
              
                vec1 = A_representors[n1][qa2-1];
                vec2 = A_representors[N-1][qa1-1];
                  
                A_n1_N(qa2, qa1) = test_inner_product(vec1, vec2);                              
              }
            }
          }
        }
        if(n1 == N-1){
          std::vector<flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > > newvec;
          newvec.push_back(A_n1_N);
          rb->A_A_representor_norms.push_back(newvec);
        }
        else{
          rb->A_A_representor_norms[n1].push_back(A_n1_N);          
        }
        std::cout << "n1 = " << n1 << ", n2 = " << N-1 << " " <<  rb->A_A_representor_norms[n1][N-1-n1] << std::endl;
    }
    timer.stop();
    std::cout << "  " << timer.elapsed() << " seconds" << std::endl << std::endl;

    // A_F representor norms 
    std::cout << "  ... A x F .... " << std::endl;
    timer.start();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_F(Qa, Qf);
    for(int qa = 1; qa <= Qa; ++qa) {
        for(int qf = 1; qf <= Qf; ++qf) {
           // std::cout << "....A_" << qa <<"(" << n << "), F_" << qf << std::endl;
           
            vec1 = A_representors[N-1][qa-1];
            vec2 = F_representors[qf-1];
            
            A_F(qa, qf) = test_inner_product(vec1, vec2);    
        }
    }
    rb->A_F_representor_norms.push_back(A_F);
    std::cout << "n = " << N-1 << " "<< rb->A_F_representor_norms[N-1] << std::endl ;

    timer.stop();
    std::cout << "  " << timer.elapsed() << " seconds" << std::endl << std::endl;
    
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::write_riesz_representors(const std::string& directory_name)
{
  // Make a directory to store all the data files
  if( mkdir(directory_name.c_str(), 0777) == -1)
  {
    std::cerr << "In RBModel::write_RB_data, directory "
                 << directory_name << " already exists, overwriting contents." << std::endl;
  }
  
  for(unsigned int i = 0; i < rb->Q_f(); ++i){
    std::stringstream filename;
    filename << directory_name << "/F_representor_" << i+1 << ".dat";
    saveCoeffVector2D(F_representors[i], testbasis, filename.str().c_str());
  }
  
  for(unsigned int n = 0; n < rb->n_bf(); ++n){
    for(unsigned int i = 0; i < rb->Q_a(); ++i){
      std::stringstream filename;
      filename << directory_name << "/A_representor_" << i+1 << "_" << n+1 << ".dat";
      saveCoeffVector2D(A_representors[n][i], testbasis, filename.str().c_str()); 
    }
  }
  
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::trial_inner_product(const CoeffVector& v1, const CoeffVector& v2)
{
    T val = 0;
    
    if(use_inner_product_matrix && assembled_inner_product_matrix){
        // Assumption here: both vectors and the matrix have the same indexset
        
        assert(v1.size() == v2.size());
        typename CoeffVector::const_iterator it1, it2;
        
        // Build dense vectors
        DenseVectorT v1_dense(v1.size()), v2_dense(v2.size());
        int index_count = 1;
        for (it1 = v1.begin(), it2 = v2.begin(); it1 != v1.end(); ++it1, ++it2, ++index_count) {
            v1_dense(index_count) = (*it1).second;
            v2_dense(index_count) = (*it2).second;
            
        }
        
        DenseVectorT I_v2 = trial_inner_product_matrix * v2_dense;
        val = v1_dense * I_v2;
    }
    else{
        typename CoeffVector::const_iterator it1, it2;
        for (it1 = v1.begin(); it1 != v1.end() ; ++it1) {
            for (it2 = v2.begin(); it2 != v2.end(); ++it2) {
                val += (*it1).second * (*trial_inner_product_op)((*it1).first, (*it2).first) * (*it2).second;
            }
        }
    }
    return val;   
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::test_inner_product(const CoeffVector& v1, const CoeffVector& v2)
{
    T val = 0;
    
    if(use_inner_product_matrix && assembled_inner_product_matrix){
        // Assumption here: both vectors and the matrix have the same indexset
        
        assert(v1.size() == v2.size());
        typename CoeffVector::const_iterator it1, it2;
        
        // Build dense vectors
        DenseVectorT v1_dense(v1.size()), v2_dense(v2.size());
        int index_count = 1;
        for (it1 = v1.begin(), it2 = v2.begin(); it1 != v1.end(); ++it1, ++it2, ++index_count) {
            v1_dense(index_count) = (*it1).second;
            v2_dense(index_count) = (*it2).second;
            
        }
        
        DenseVectorT I_v2 = test_inner_product_matrix * v2_dense;
        val = v1_dense * I_v2;
    }
    else{
        typename CoeffVector::const_iterator it1, it2;
        for (it1 = v1.begin(); it1 != v1.end() ; ++it1) {
            for (it2 = v2.begin(); it2 != v2.end(); ++it2) {
                val += (*it1).second * (*test_inner_product_op)((*it1).first, (*it2).first) * (*it2).second;
            }
        }
    }
    return val;   
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::inner_product(const CoeffVector& v1, const CoeffVector& v2)
{
    return trial_inner_product(v1, v2);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::assemble_inner_product_matrix(IndexSet<Index2D>& trial_indexset, IndexSet<Index2D>& test_indexset)
{    
    Timer timer;
    std::cout << "Assemble Inner Product Matrices ...." << std::endl;
    
    trial_inner_product_matrix.resize((int)trial_indexset.size(), (int)trial_indexset.size());
    test_inner_product_matrix.resize((int)test_indexset.size(), (int)test_indexset.size());
    timer.start();
    
    // ToFlensSparseMatrix Trial Indexset
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    std::map<Index2D,int,lt<Lexicographical,Index2D> > row_indices;
    int row_count = 1, col_count = 1;
    for (const_set_it row=trial_indexset.begin(); row!=trial_indexset.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    for (const_set_it col=trial_indexset.begin(); col!=trial_indexset.end(); ++col, ++col_count) {
        for (const_set_it row=trial_indexset.begin(); row!=trial_indexset.end(); ++row, ++row_count) {
            T tmp = (*trial_inner_product_op)(*row,*col);
            if (fabs(tmp)>0)                trial_inner_product_matrix(row_indices[*row],col_count) = tmp;
        }
    }
    trial_inner_product_matrix.finalize();
    
                //toFlensSparseMatrix(repr_lhs_op, test_indexset, test_indexset, test_inner_product_matrix);
    
    // ToFlensSparseMatrix Test Indexset
    row_indices.clear();    
    row_count = 1;
    col_count = 1;
    
    if(!assembled_prec_vec){
        test_prec_vec.engine().resize((int)test_indexset.size());
        for (const_set_it row=test_indexset.begin(); row!=test_indexset.end(); ++row, ++row_count) {
            row_indices[(*row)] = row_count;
            test_prec_vec(row_count) = get_test_prec((*row));
        }
        assembled_prec_vec = true;
    }
    else{
        for (const_set_it row=test_indexset.begin(); row!=test_indexset.end(); ++row, ++row_count) {
            row_indices[(*row)] = row_count;
        }    
    }
    
    repr_lhs_op.compression.setParameters(test_indexset);
    for (const_set_it col=test_indexset.begin(); col!=test_indexset.end(); ++col, ++col_count) {
        IndexSet<Index2D> LambdaRowSparse = repr_lhs_op.compression.SparsityPattern(*col, test_indexset);
        for (const_set_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row) {
            T tmp = (*test_inner_product_op)(*row,*col);
            if (fabs(tmp)>0)                test_inner_product_matrix(row_indices[*row],col_count) = tmp;
        }
    }
    test_inner_product_matrix.finalize();
    
    timer.stop();
    std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;
    
    assembled_inner_product_matrix = true;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::assemble_A_operator_matrices(IndexSet<Index2D>& trial_indexset, IndexSet<Index2D>& test_indexset)
{
  Timer timer;
  std::cout << "Assemble A Matrices ...." << std::endl;
  unsigned int Q_a = get_rb_model().Q_a();
  int N1 = (int)trial_indexset.size();
  int N2 = (int)test_indexset.size();
  
  timer.start();
  for(unsigned int qa = 0; qa < Q_a; ++qa){
    SparseMatrixT A_matrix(N1, N2);
    
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    std::map<Index2D,int,lt<Lexicographical,Index2D> > row_indices;
    int row_count = 1, col_count = 1;
    
        if((!assembled_prec_vec) && qa == 0){
          test_prec_vec.engine().resize((int)test_indexset.size());
          for (const_set_it row=test_indexset.begin(); row!=test_indexset.end(); ++row, ++row_count) {
              row_indices[(*row)] = row_count;
              test_prec_vec(row_count) = get_test_prec((*row));
          }
          assembled_prec_vec = true;
      }
      else{
          for (const_set_it row=test_indexset.begin(); row!=test_indexset.end(); ++row, ++row_count) {
              row_indices[(*row)] = row_count;
          }    
      }
      
    for (const_set_it col=trial_indexset.begin(); col!=trial_indexset.end(); ++col, ++col_count) {
        for (const_set_it row=trial_indexset.begin(); row!=trial_indexset.end(); ++row, ++row_count) {
            T tmp = (*A_operators[qa])(*row,*col);
            if (fabs(tmp)>0)                A_matrix(row_indices[*row],col_count) = tmp;
        }
    }
    A_matrix.finalize();
      
    A_operator_matrices.push_back(A_matrix);
  }
  timer.stop();
  std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;
  
  assembled_A_operator_matrices = true;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::uncached_residual_dual_norm(const DenseVectorT& u_RB, const std::vector<T>& mu)
{
  // Calculate residual
  Operator_Residual_Representor op_res_repr(this, mu, u_RB);

  // Solve for Riesz representor
  Coefficients<Lexicographical,T,Index2D> res_repr = solver->repr_solve_totalRes(op_res_repr);
  
  // Calculate Norm
  undo_test_prec(res_repr);
  
  return std::sqrt(inner_product(res_repr, res_repr));
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::
uncached_residual_dual_norm(const DenseVectorT& u_RB, const std::vector<T>& mu, Coefficients<Lexicographical,T,Index2D>& res_repr)
{
  // Calculate residual
  Operator_Residual_Representor op_res_repr(this, mu, u_RB);

  // Solve for Riesz representor
  res_repr = solver->repr_solve_totalRes(op_res_repr);
  
  // Calculate Norm
  undo_test_prec(res_repr);
  return std::sqrt(test_inner_product(res_repr, res_repr));
}

// ================================================================================================================ //
// ================================================================================================================ //

/*  Operator LHS */

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_LHS::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    if(qa < 0){
      T val = 0;
      for (unsigned int i = 0; i < thisTruth->A_operators.size(); ++i) {
          val += (*thisTruth->rb->theta_a[i])(thisTruth->rb->get_current_param()) 
               * (*thisTruth->A_operators[i])(row_index, col_index);
      }
      return thisTruth->get_trial_prec(col_index) * thisTruth->get_test_prec(row_index) * val;
    }
    else{
      return thisTruth->get_trial_prec(col_index) * thisTruth->get_test_prec(row_index) * (*thisTruth->A_operators[qa])(row_index, col_index);
    }
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_LHS::
mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &x)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS::mv not implemented."
              << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_LHS::
toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow,const IndexSet<Index2D> &LambdaCol,
                    SparseMatrixT &A, T tol)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS::toFlensSparseMatrix "
              << "not implemented."
              << std::endl;
    assert(0);
    exit(1);
}

/*  Operator RHS */

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS::operator()(const Index2D &lambda)
{
    T val = 0;
    for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
        val += (*thisTruth->rb->theta_f[i])(thisTruth->rb->get_current_param()) 
             * (*thisTruth->F_operators[i])(lambda);
    }
    
    return thisTruth->get_test_prec(lambda) * val;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    typename IndexSet<Index2D>::const_iterator it_Lambda;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
        coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
    }
    
    return coeffs;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS::operator()(T tol)
{
    Coefficients<Lexicographical,T,Index2D> c;
       
    std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
    exit(1);
    
    
    return c;
}

/*  Operator LHS_Representor */

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_LHS_Representor::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    return thisTruth->get_test_prec(row_index) * thisTruth->get_test_prec(col_index) * (*thisTruth->test_inner_product_op)(row_index, col_index);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_LHS_Representor::
mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &x)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS_Representor::mv "
              << "not implemented." << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_LHS_Representor::
toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol,
                    SparseMatrixT &A, T tol)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS_Representor::"
              <<"toFlensSparseMatrix not implemented." << std::endl;
    assert(0);
    exit(1);
}

/* Operator RHS_BilFormRepresentor */

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(const Index2D &lambda)
{
    T val = 0;
    typename Coefficients<Lexicographical,T,Index2D>::const_iterator it;
    for (it = current_bf->begin(); it != current_bf->end(); ++it) {
        //val += (*it).second * (*current_op)((*it).first, lambda) ;
        val += (*it).second * (*current_op)(lambda, (*it).first) ;
    }
    
    return - thisTruth->get_test_prec(lambda) * val;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    typename IndexSet<Index2D>::const_iterator it_Lambda;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
        coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
    }
    
    return coeffs;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(T tol)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
    exit(1);
    
    return coeffs;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::
set_current_op(Operator2D<T>& op)
{
    current_op = &op;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::
set_current_bf(Coefficients<Lexicographical,T,Index2D>& bf)
{
    current_bf = &bf;
}

/* Operator RHS_FunctionalRepresentor */

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(const Index2D &lambda)
{
    return thisTruth->get_test_prec(lambda) * (*current_op)(lambda);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(const IndexSet<Index2D> &Lambda)
{
    CoeffVector coeffs;
    
    typename IndexSet<Index2D>::const_iterator it_Lambda;
    typedef typename CoeffVector::value_type val_type;
    for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
        coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
    }
    
    return coeffs;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(T tol)
{
    std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
    exit(1);
    
    return (*current_op)(tol);
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::
set_current_op(AdaptiveRhs<T, Index2D>& op)
{
    current_op = &op;
}

/* Operator Residual_Representor */

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_Residual_Representor::operator()(const Index2D &lambda)
{
  T val = 0;

  typename Coefficients<Lexicographical,T,Index2D>::const_iterator it;

  for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
      val += (*thisTruth->rb->theta_f[i])(eval_mu) * (*thisTruth->F_operators[i])(lambda);
  }
   

  for (unsigned int i = 0; i < thisTruth->A_operators.size(); ++i) {
    T val_Au = 0;
    
    for (int j = 0; j < u_N.length();++j){
      T val_A_bf = 0;
      for (it = thisTruth->rb->rb_basis_functions[j].begin(); it != thisTruth->rb->rb_basis_functions[j].end(); ++it) {
          //val_A_bf += (*it).second *(*thisTruth->A_operators[i])((*it).first,lambda);
          val_A_bf += (*it).second *(*thisTruth->A_operators[i])(lambda, (*it).first);
      }
      val_Au += u_N(j+1) * val_A_bf;
    }
    
    val -= (*thisTruth->rb->theta_a[i])(eval_mu) * val_Au;
  }
  
  return thisTruth->get_test_prec(lambda) * val;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_Residual_Representor::operator()(const IndexSet<Index2D> &Lambda)
{
  Coefficients<Lexicographical, T, Index2D> coeffs;
  
  typename IndexSet<Index2D>::const_iterator it_Lambda;
  typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
  for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
      coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
  }
  
  return coeffs;
}

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec,  typename TestPrec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>::Operator_Residual_Representor::operator()(T tol)
{
  Coefficients<Lexicographical, T, Index2D> coeffs;
  
  std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
  exit(1);
  
  return coeffs;
}



} // namespace lawa