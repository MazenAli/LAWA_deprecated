#include <cassert>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>

namespace  lawa {

template <typename T, typename TruthModel>
RBModel2D<T, TruthModel>::RBModel2D()
 : assembled_A_operator_matrices(false)
{
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::set_current_param(const std::vector<T>& _param)
{
  current_param = _param;
}

template <typename T, typename TruthModel>
std::vector<T>&
RBModel2D<T, TruthModel>::get_current_param()
{
  return current_param;
}

template <typename T, typename TruthModel>
unsigned int
RBModel2D<T, TruthModel>::Q_a()
{
  return theta_a.size();  
}
        
template <typename T, typename TruthModel>
unsigned int
RBModel2D<T, TruthModel>::Q_f()
{
  return theta_f.size();  
}

template <typename T, typename TruthModel>
unsigned int
RBModel2D<T, TruthModel>::n_bf()
{
  return rb_basis_functions.size();  
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::attach_theta_a_q(theta_fctptr theta_a_q)
{
  theta_a.push_back(theta_a_q);
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::attach_theta_f_q(theta_fctptr theta_f_q)
{
  theta_f.push_back(theta_f_q);
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::set_truthmodel(TruthModel& _truthmodel)
{
  truth = &_truthmodel;
    truth->set_rb_model(*this);
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::add_to_basis(const CoeffVector& sol)
{
    std::cout << "Adding basis function .... " << std::endl;
    
    Timer timer;
    
    std::cout << "  Gram-Schmidt...  " << std::endl;
    timer.start();
    truth->add_new_basis_function(sol);
    timer.stop();
    std::cout << "  ... done: " << timer.elapsed() << " seconds" << std::endl << std::endl;
    
    timer.start();
    std::cout << "  Updating RB Matrices ..." << std::endl;
    update_RB_A_matrices();
    update_RB_F_vectors();
    update_RB_inner_product();
    timer.stop();
    std::cout << "  Done updating RB Matrices: " << timer.elapsed() << " seconds" << std::endl << std::endl;
    
    std::cout << "  Updating Representors ..." << std::endl;
    timer.start();
    truth->update_representors();
    timer.stop();
    std::cout << "  Done updating representors: " << timer.elapsed() << " seconds" << std::endl << std::endl;
    
    std::cout << "  Calculating Representor norms ..." << std::endl;
    timer.start();
    truth->update_representor_norms(); // Funktion kann eigetnlich hierher verschoben werden
    timer.stop();
    std::cout << "  Done calculating representor norms: " << timer.elapsed() << " seconds" << std::endl << std::endl;
    
    std::cout << "  .............................." << std::endl << std::endl;
}

template <typename T, typename TruthModel>
flens::DenseVector<flens::Array<T> >
RBModel2D<T, TruthModel>::RB_solve(unsigned int N, std::vector<T> param, SolverCall call)
{
  assert(RB_A_matrices.size() > 0);
    
  FullColMatrixT A(N, N);
    for (unsigned int i = 0; i < Q_a(); ++i) {
        //A += (*theta_a[i])(get_current_param()) * RB_A_matrices[i](_(1,N), _(1,N));
        flens::axpy(cxxblas::NoTrans, (*theta_a[i])(param), RB_A_matrices[i](_(1,N), _(1,N)), A);
    }

    DenseVectorT F(N);
    for (unsigned int i = 0; i < Q_f(); ++i) {
        //F += (*theta_f[i])(get_current_param()) * RB_F_vectors[i](_(1,N));
        flens::axpy((*theta_f[i])(param), RB_F_vectors[i](_(1,N)), F);
    }
    
    DenseVectorT u(N);
    switch (call) {
        case call_cg:
            std::cout << "RB solve: " << cg(A, u, F) << " cg iterations" << std::endl;
            break;
        case call_gmres:
            std::cout << "RB solve: " << gmres(A, u, F) << " gmres iterations" << std::endl;
            break;
        case call_cgls:
            std::cout << "RB solve: " << cgls(A, u, F) << " cgls iterations" << std::endl;
            break;
        default:
          std::cerr << "RB solve: Method not implemented! " << std::endl;
            break;
    }
    
    return u;
}

template <typename T, typename TruthModel>
flens::DenseVector<flens::Array<T> >
RBModel2D<T, TruthModel>::RB_solve(unsigned int N, SolverCall call)
{
  assert(RB_A_matrices.size() > 0);
    
  FullColMatrixT A(N, N);
    for (unsigned int i = 0; i < Q_a(); ++i) {
        //A += (*theta_a[i])(get_current_param()) * RB_A_matrices[i](_(1,N), _(1,N));
        flens::axpy(cxxblas::NoTrans, (*theta_a[i])(get_current_param()), RB_A_matrices[i](_(1,N), _(1,N)), A);
    }

    DenseVectorT F(N);
    for (unsigned int i = 0; i < Q_f(); ++i) {
        //F += (*theta_f[i])(get_current_param()) * RB_F_vectors[i](_(1,N));
        flens::axpy((*theta_f[i])(get_current_param()), RB_F_vectors[i](_(1,N)), F);
    }
    
    DenseVectorT u(N);
    switch (call) {
        case call_cg:
            std::cout << "RB solve: " << cg(A, u, F) << " cg iterations" << std::endl;
            break;
        case call_gmres:
            std::cout << "RB solve: " << gmres(A, u, F) << " gmres iterations" << std::endl;
            break;
        default:
          std::cerr << "RB solve: Method not implemented! " << std::endl;
            break;
    }
    
    return u;
}

template <typename T, typename TruthModel>
Coefficients<Lexicographical,T,Index2D>
RBModel2D<T, TruthModel>::reconstruct_u_N(DenseVectorT u, unsigned int N)
{
  assert(N <= n_bf());
  assert(u.length() > 0);
  assert(u.length() >= (int)N);
    
  CoeffVector u_full;
    for (unsigned int i = 1; i <= N; ++i) {
        u_full = u_full + u(i) * rb_basis_functions[i-1];
    }
    
    return u_full;
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::set_min_param(const std::vector<T>& _param)
{
    min_param = _param;
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::set_max_param(const std::vector<T>& _param)
{
    max_param = _param;
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::set_ref_param(const std::vector<T>& _param)
{
    ref_param = _param;
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::generate_uniform_trainingset(std::vector<int>& n_train)
{
    std::vector<T> h(current_param.size());
    for(unsigned int d = 0; d < current_param.size(); ++d) {
        h[d] = (max_param[d] - min_param[d]) / (n_train[d]-1);
    }
    
    if(n_train.size() == 1) {
        for (int i = 0; i < n_train[0]; ++i){
            std::vector<T> new_mu;
            new_mu.push_back(std::min(min_param[0] + i*h[0], max_param[0]));
            Xi_train.push_back(new_mu);   
        }
    }
    else {
        std::cerr << "Generate Trainingsset for dim = " << n_train.size() << " : Not implemented yet " << std::endl;
    }
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::train_Greedy(const std::vector<T>& init_param, T tol, int Nmax, const char* filename, SolverCall call)
{
    // Initial Snapshot
    set_current_param(init_param);
    
    std::ofstream error_file(filename); 
    error_file << "# N   mu   Error " << std::endl;    
    
    // Delete init parameter from training sample
    for(unsigned int i = 0; i < Xi_train.size(); ++i) {
        bool is_param = false;
        for(unsigned int l = 0; l < init_param.size(); ++l) {
            if(Xi_train[i][l] == init_param[l]){
                is_param = true;
            }
            else {
                is_param = false;
                break;
            }
        }
        if( is_param == true) {
            Xi_train.erase(Xi_train.begin() + i);
        }
    }
    
    T maxerr;
    int N = 0;
    do {
        std::cout << " ================================================= " << std::endl << std::endl;
        std::cout << "Adding Snapshot at mu = " << get_current_param()[0] << std::endl << std::endl;
        error_file << N+1 << " " << get_current_param()[0];
        
        CoeffVector u = truth->solver->truth_solve();
        add_to_basis(u);
        N++;
        
        /*
        // Scatterplot
        std::cout << "Plotting Scatterplot Truth .... " << std::endl;
        std::stringstream s;
        s << "truth_" << N;
        saveCoeffVector2D(u, truth->basis, s.str().c_str());
        std::cout << ".... finished " << std::endl;
        
        // Scatterplot
        std::cout << "Plotting Scatterplot BasisFunction .... " << std::endl;
        std::stringstream s_bf;
        s_bf << "bf_" << N;
        saveCoeffVector2D(rb_basis_functions[N-1], truth->basis, s_bf.str().c_str());
        std::cout << ".... finished " << std::endl;
        */
        
        std::cout << std::endl;
        maxerr = 0;
        int next_Mu = 0;
        std::cout << "Start Greedy search for new parameter " << std::endl;
        for (unsigned int n = 0; n < Xi_train.size(); ++n) {
            
            set_current_param(Xi_train[n]);
            std::cout << "    Parameter set to " << Xi_train[n][0] << std::endl;
            DenseVectorT u_N = RB_solve(N, Xi_train[n], call);
            
            T resnorm = residual_dual_norm(u_N, Xi_train[n]);
            T alpha =  alpha_LB(Xi_train[n]);
            T error_est = resnorm / alpha;
            
            std::cout << "Training parameter " << Xi_train[n][0]  << ": Error = " << std::scientific << error_est 
                      << " = " << resnorm << " / " << alpha << std::endl;
            if ( error_est > maxerr) {
                maxerr = error_est;
                next_Mu = n;
            }
        }
        
        error_file << " " << maxerr << std::endl;
        
        std::cout << std::endl << "Greedy Error = " << std::scientific << maxerr << std::endl << std::endl;
        
        set_current_param(Xi_train[next_Mu]);
        Xi_train.erase(Xi_train.begin() + next_Mu);
        
    } while ((N < Nmax) && (maxerr > tol));
    
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::write_basis_functions(const std::string& directory_name){
  // Make a directory to store all the data files
  if( mkdir(directory_name.c_str(), 0777) == -1)
  {
    std::cerr << "In RBModel::write_basis_functions, directory "
                 << directory_name << " already exists, overwriting contents." << std::endl;
  }
  
  std::stringstream n_bf_filename;
  n_bf_filename << directory_name << "/n_bf.dat";
  std::ofstream n_bf_file(n_bf_filename.str().c_str());
  n_bf_file <<  n_bf() << std::endl;
  n_bf_file.close();
  
  for(unsigned int i = 0; i < n_bf(); ++i){
    std::stringstream filename;
    filename << directory_name << "/bf_" << i+1 << ".dat";
    saveCoeffVector2D(rb_basis_functions[i], truth->trialbasis, filename.str().c_str());
  }
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::read_basis_functions(const std::string& directory_name){
  
  unsigned int n_bf;
  std::stringstream n_bf_filename;
  n_bf_filename << directory_name << "/n_bf.dat";
  
  std::ifstream n_bf_file(n_bf_filename.str().c_str());
  if(n_bf_file.is_open()){
    n_bf_file >> n_bf;
    n_bf_file.close();
  }
  else{
    std::cerr << "Unable to read number of basis functions! " << std::endl;
    exit(1);
  }
  
  rb_basis_functions.clear();
  for(unsigned int i = 1; i <= n_bf; ++i){
    std::stringstream filename;
    filename << directory_name << "/bf_" << i << ".dat";
    CoeffVector bf_coeffs;
    readCoeffVector2D(bf_coeffs, filename.str().c_str(),false); 
    rb_basis_functions.push_back(bf_coeffs);
    std::cout << " Read " << filename.str() << std::endl;
  }
}



template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::write_RB_data(const std::string& directory_name){
  
  const unsigned int precision = 16;
  
  // Make a directory to store all the data files
  if( mkdir(directory_name.c_str(), 0777) == -1)
  {
    std::cerr << "In RBModel::write_RB_data, directory "
                 << directory_name << " already exists, overwriting contents." << std::endl;
  }
  
  std::stringstream n_bf_filename;
  n_bf_filename << directory_name << "/n_bf.dat";
  std::ofstream n_bf_file(n_bf_filename.str().c_str());
  n_bf_file << n_bf() << std::endl;
  n_bf_file.close();
  
  // Write RB_A_matrices
  for(unsigned int i = 0; i < Q_a(); ++i){
    std::stringstream filename;
    filename << directory_name << "/RB_A_" << i+1 << ".dat";
    std::ofstream file(filename.str().c_str());
    file.precision(precision);
    file << std::scientific << RB_A_matrices[i] << std::endl;
    file.close();
  }
  
  // Write RB_F_vectors
  for(unsigned int i = 0; i < Q_f(); ++i){
    std::stringstream filename;
    filename << directory_name << "/RB_F_" << i+1 << ".dat";
    std::ofstream file(filename.str().c_str());
    file.precision(precision);
    file << std::scientific << RB_F_vectors[i] << std::endl;
    file.close();
  }
  
  // Write RB_inner_product
  std::stringstream filename;
  filename << directory_name << "/RB_inner_product.dat";
  std::ofstream file(filename.str().c_str());
  file.precision(precision);
  file << std::scientific << RB_inner_product << std::endl;
  file.close();
  
  // Write F_F_representor_norms
  std::stringstream repr_F_filename;
  repr_F_filename << directory_name << "/F_F_representor_norms.dat";
  std::ofstream repr_F_file(repr_F_filename.str().c_str());
  repr_F_file.precision(precision);
  repr_F_file << std::scientific << F_F_representor_norms << std::endl;
  repr_F_file.close();
  
  // Write A_F_representor_norms
  std::stringstream repr_A_F_filename;
  repr_A_F_filename << directory_name << "/A_F_representor_norms.dat";
  std::ofstream repr_A_F_file(repr_A_F_filename.str().c_str());
  repr_A_F_file.precision(precision);
  for(unsigned int i = 0; i < n_bf(); ++i){    
    repr_A_F_file << std::scientific << A_F_representor_norms[i] << std::endl;
  }
  repr_A_F_file.close();
  
  // Write A_A_representor_norms
  std::stringstream repr_A_A_filename;
  repr_A_A_filename << directory_name << "/A_A_representor_norms.dat";
  std::ofstream repr_A_A_file(repr_A_A_filename.str().c_str());
  repr_A_A_file.precision(precision);
  for(unsigned int i = 0; i < n_bf(); ++i){    
    for(unsigned int j = i; j < n_bf(); ++j){    
        repr_A_A_file << std::scientific << A_A_representor_norms[i][j-i] << std::endl;
    }
  }
  repr_A_A_file.close();
  
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::read_RB_data(const std::string& directory_name){
  
  // Read Nr of BasisFunctions
  unsigned int n_bf;
  std::stringstream n_bf_filename;
  n_bf_filename << directory_name << "/n_bf.dat";
  
  std::ifstream n_bf_file(n_bf_filename.str().c_str());
  if(n_bf_file.is_open()){
    n_bf_file >> n_bf;
    n_bf_file.close();
    
    std::cout << "Number of Basis Functions: " << n_bf << std::endl;
  }
  else{
    std::cerr << "Unable to read number of basis functions! " << std::endl;
    exit(1);
  }
  
  // Read RB_A_matrices
  RB_A_matrices.clear();
  for(unsigned int i = 0; i < Q_a(); ++i){
    std::stringstream filename;
    filename << directory_name << "/RB_A_" << i+1 << ".dat";
    std::ifstream file(filename.str().c_str());
    if(file.is_open()){
      FullColMatrixT RB_A(n_bf, n_bf);
      for(unsigned int n1 = 1; n1 <= n_bf; n1++){
        for(unsigned int n2 = 1; n2 <= n_bf; n2++){
          file >> RB_A(n1, n2);
        }
      }
      file.close();
      RB_A_matrices.push_back(RB_A);
      std::cout << " Read " << filename.str() << std::endl;
    }
    else{
      std::cerr << "Unable to open file " << filename.str() << " for reading!" << std::endl;
      exit(1);
    }
  }
  
  // Read RB_F_vectors
  RB_F_vectors.clear();
  for(unsigned int i = 0; i < Q_f(); ++i){
    std::stringstream filename;
    filename << directory_name << "/RB_F_" << i+1 << ".dat";
    std::ifstream file(filename.str().c_str());
    if(file.is_open()){
      DenseVectorT RB_F(n_bf);
      for(unsigned int n = 1; n <= n_bf; n++){
          file >> RB_F(n);
      }
      file.close();
      RB_F_vectors.push_back(RB_F);
      std::cout << " Read " << filename.str() << std::endl;
    }
    else{
      std::cerr << "Unable to open file " << filename.str() << " for reading!" << std::endl;
      exit(1);
    }
  }
  
  // Read RB_inner_product
  RB_inner_product.engine().resize((int)n_bf, (int)n_bf);
  std::stringstream inner_product_filename;
  inner_product_filename << directory_name << "/RB_inner_product.dat";
  std::ifstream inner_product_file(inner_product_filename.str().c_str());
  if(inner_product_file.is_open()){
    for(unsigned int n1 = 1; n1 <= n_bf; n1++){
      for(unsigned int n2 = 1; n2 <= n_bf; n2++){
        inner_product_file >> RB_inner_product(n1, n2);
      }
    }
    inner_product_file.close();
    std::cout << " Read " << inner_product_filename.str() << std::endl;

  }
  else{
    std::cerr << "Unable to open file " << inner_product_filename.str() << " for reading!" << std::endl;
    exit(1);
  }
  
  // Read F_F_representor norms
  F_F_representor_norms.engine().resize((int)Q_f(), (int)Q_f());
  std::stringstream F_F_filename;
  F_F_filename << directory_name << "/F_F_representor_norms.dat";
  std::ifstream F_F_file(F_F_filename.str().c_str());
  if(F_F_file.is_open()){
    for(unsigned int i = 1; i <= Q_f(); ++i){
      for(unsigned int j = 1; j <= Q_f(); ++j){
        F_F_file >> F_F_representor_norms(i,j);
      }
    }
    F_F_file.close();
    std::cout << " Read " << F_F_filename.str() << std::endl;
  }
  else{
    std::cerr << "Unable to open file " << F_F_filename.str() << " for reading!" << std::endl;
    exit(1);
  }
  
  // Read A_F_representor norms
  A_F_representor_norms.clear();
  std::stringstream A_F_filename;
  A_F_filename << directory_name << "/A_F_representor_norms.dat";
  std::ifstream A_F_file(A_F_filename.str().c_str());
  if(A_F_file.is_open()){
    for(unsigned int n = 1; n <= n_bf; ++n){
      FullColMatrixT A_F_matrix(Q_a(), Q_f());
      for(unsigned int i = 1; i <= Q_a(); ++i){
        for(unsigned int j = 1; j <= Q_f(); ++j){
          A_F_file >> A_F_matrix(i,j);
        }
      }
      A_F_representor_norms.push_back(A_F_matrix);
    }
    A_F_file.close();
    std::cout << " Read " << A_F_filename.str() << std::endl;
  }
  else{
    std::cerr << "Unable to open file " << A_F_filename.str() << " for reading!" << std::endl;
    exit(1);
  }
  
  // Read A_A_representor norms
  A_A_representor_norms.clear();
  std::stringstream A_A_filename;
  A_A_filename << directory_name << "/A_A_representor_norms.dat";
  std::ifstream A_A_file(A_A_filename.str().c_str());
  if(A_A_file.is_open()){
    for(unsigned int n1 = 1; n1 <= n_bf; ++n1){
      std::vector<FullColMatrixT> A_A_n_vec;
      for(unsigned int n2 = n1; n2 <= n_bf; ++n2){
        FullColMatrixT A_A_matrix(Q_a(), Q_a());
        for(unsigned int i = 1; i <= Q_a(); ++i){
          for(unsigned int j = 1; j <= Q_a(); ++j){
            A_A_file >> A_A_matrix(i,j);
          }
        }
        A_A_n_vec.push_back(A_A_matrix);
      }
      A_A_representor_norms.push_back(A_A_n_vec);
    }
    A_A_file.close();
    std::cout << " Read " << A_A_filename.str() << std::endl;
  }
  else{
    std::cerr << "Unable to open file " << A_A_filename.str() << " for reading!" << std::endl;
    exit(1);
  }
}


// -----------------------------------------------------------------------------------//

/* Protected methods */

template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::update_RB_A_matrices()
{
  typename CoeffVector::const_iterator it, it1, it2;
    if ((n_bf() == 1) && (RB_A_matrices.size() < Q_a())) {
        //RB_A_matrices.resize(Q_a());
        for (unsigned int q_a = 0; q_a < Q_a(); ++q_a) {
            FullColMatrixT A(1,1);
            RB_A_matrices.push_back(A);
        }
    }
    
    // Assemble new basis function vector
    DenseVectorT new_bf_dense;
    if(assembled_A_operator_matrices){
        // Build dense vector
        new_bf_dense.engine().resize((int)rb_basis_functions[n_bf()-1].size());
        int index_count = 1;
        for (it = rb_basis_functions[n_bf()-1].begin(); it != rb_basis_functions[n_bf()-1].end(); ++it, ++index_count) {
          new_bf_dense(index_count) = (*it).second;
        }
        
        //std::cout << "new_bf_dense: " << new_bf_dense << std::endl;
    }
    
    for (unsigned int q_a = 0; q_a < Q_a(); ++q_a) {
        if (n_bf() == 1) {
            RB_A_matrices[q_a].engine().resize(1, 1);
        }
        else {
            FullColMatrixT tmp(RB_A_matrices[q_a]);
            RB_A_matrices[q_a].engine().resize((int)n_bf(), (int)n_bf());
            RB_A_matrices[q_a](tmp.rows(), tmp.cols()) = tmp;
            for(unsigned int i = 1; i <= n_bf(); ++i) {
                RB_A_matrices[q_a](i, n_bf()) = 0.;
                RB_A_matrices[q_a](n_bf(), i) = 0.;
            }

        }
        
        DenseVectorT A_new_bf, A_new_bf_T;
        if(assembled_A_operator_matrices){
            A_new_bf = truth->A_operator_matrices[q_a] * new_bf_dense;
            A_new_bf_T = transpose(truth->A_operator_matrices[q_a]) * new_bf_dense;
            //std::cout << A_new_bf << std::endl;
            
            //std::cout << A_new_bf_T << std::endl;
        }
        
        for (unsigned int i = 1; i <= n_bf(); ++i) {
            if(assembled_A_operator_matrices){
                DenseVectorT bf_dense(rb_basis_functions[i-1].size());
                int index_count = 1;
                for (it = rb_basis_functions[i-1].begin(); it != rb_basis_functions[i-1].end(); ++it, ++index_count) {
                  bf_dense(index_count) = (*it).second;
                }
                RB_A_matrices[q_a](i, n_bf()) = bf_dense * A_new_bf;
                if(i != n_bf()){
                    RB_A_matrices[q_a](n_bf(), i) = A_new_bf_T * bf_dense;                    
                }
            }
            else{
                for (it1 = rb_basis_functions[n_bf()-1].begin(); it1 != rb_basis_functions[n_bf()-1].end(); ++it1) {
                    for (it2 = rb_basis_functions[i-1].begin(); it2 != rb_basis_functions[i-1].end(); ++it2) {
                        RB_A_matrices[q_a](n_bf(), i) += (*it1).second * (*it2).second 
                                                   * (*truth->A_operators[q_a])((*it1).first, (*it2).first);
                        if (i != n_bf()) {
                         RB_A_matrices[q_a](i, n_bf()) += (*it1).second * (*it2).second 
                                   * (*truth->A_operators[q_a])((*it2).first, (*it1).first);
                        }
                    }
                }                
            }
        }
        
        std::cout << "RB_A(" << q_a << ") = " << RB_A_matrices[q_a] << std::endl;
    }
    
}

template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::update_RB_F_vectors()
{
    typename CoeffVector::const_iterator it;
    if ((n_bf() == 1) && (RB_F_vectors.size() < Q_f())) {
        RB_F_vectors.resize(Q_f());
    }
    for (unsigned int q_f = 0; q_f < Q_f(); ++q_f) {
      if (n_bf() == 1) {
        RB_F_vectors[q_f].engine().resize(1);
      }
      else {
        DenseVectorT tmp(RB_F_vectors[q_f]);
        RB_F_vectors[q_f].engine().resize((int)n_bf());
        RB_F_vectors[q_f](tmp.range()) = tmp;
        RB_F_vectors[q_f](n_bf()) = 0.;            
      }

      for (it = rb_basis_functions[n_bf()-1].begin(); it != rb_basis_functions[n_bf()-1].end(); ++it) {
          RB_F_vectors[q_f](n_bf()) += (*it).second * (*truth->F_operators[q_f])((*it).first);
      }
      
      std::cout << "RB_F(" << q_f << ") = " << RB_F_vectors[q_f] << std::endl;
        
    }
}

template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::update_RB_inner_product()
{
  typename CoeffVector::const_iterator it, it1, it2;

    if (n_bf() == 1) {
        RB_inner_product.engine().resize(1, 1);
    }
    else {
        FullColMatrixT tmp(RB_inner_product);
        RB_inner_product.engine().resize((int)n_bf(), (int)n_bf());
        RB_inner_product(tmp.rows(), tmp.cols()) = tmp;
        for(unsigned int i = 1; i <= n_bf(); ++i) {
            RB_inner_product(i, n_bf()) = 0.;
            RB_inner_product(n_bf(), i) = 0.;
        }
    }
    
    for (unsigned int i = 1; i <= n_bf(); ++i) {
      RB_inner_product(n_bf(), i) = truth->inner_product(rb_basis_functions[n_bf()-1], rb_basis_functions[i-1]);
      if (i != n_bf()) {
        RB_inner_product(i, n_bf()) = RB_inner_product(n_bf(),i);
      }      
    }
}

template <typename T, typename TruthModel>
T
RBModel2D<T, TruthModel>::residual_dual_norm(const DenseVectorT& u_RB, const std::vector<T>& mu)
{
    T res_dual_norm = 0;
        
    int Qf = Q_f();
    int Qa = Q_a();
    DenseVectorT ThetaF(Qf);
    DenseVectorT ThetaA(Qa);
    for (int i = 1; i <= Qf; ++i) {
        ThetaF(i) = (*theta_f[i-1])(mu); 
    }
    for (int i = 1; i <= Qa; ++i) {
        ThetaA(i) =(*theta_a[i-1])(mu); 
    }

    DenseVectorT FF_T = F_F_representor_norms * ThetaF;
    
    res_dual_norm = ThetaF * FF_T;
        
    int N = u_RB.length();
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
    res_dual_norm += 2 * u_RB * T_AF_T;    
    
    DenseVectorT T_AA_T_u = T_AA_T * u_RB;
    res_dual_norm += u_RB * T_AA_T_u;    

  
    if(res_dual_norm < 0){
      std::cout << "Warning: Residual dual norm negative: " << std::setprecision(10) << res_dual_norm << std::endl;
      res_dual_norm = std::fabs(res_dual_norm);
    }
    
    return std::sqrt(res_dual_norm);
}

template <typename T, typename TruthModel>
T
RBModel2D<T, TruthModel>::alpha_LB(std::vector<T>& _param)
{
    T alpha_lb = DBL_MAX;
    for(unsigned int qa = 0; qa < Q_a(); ++qa){
        T reftheta = (*theta_a[qa])(ref_param);
        if (((*theta_a[qa])(_param) / reftheta) < alpha_lb) {
            alpha_lb = (*theta_a[qa])(_param) / reftheta;
        }
    }
    
    return alpha_lb;
}
} // namespace lawa

