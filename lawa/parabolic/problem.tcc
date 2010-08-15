namespace lawa{

template<typename T, typename Basis>
Problem<T, Basis>::Problem(const Basis& _basis)
	: basis(_basis)
{
}

template<typename T, typename Basis>
template <typename BilinearForm>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
Problem<T, Basis>::getStiffnessMatrix(BilinearForm& a, int J, T tol)
{   
    int j0 = basis.j0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-1;
    int offsetI = basis.mra.rangeI(j0).firstIndex()-1;
    
    int N = basis.mra.cardI(J);
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A(N,N);
    
    // SF * SF 
    for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
        for(int k2 = basis.mra.rangeI(j0).firstIndex(); k2 <= basis.mra.rangeI(j0).lastIndex(); ++k2){            
            T val = a(XBSpline, j0, k1, XBSpline, j0, k2);            
            if(fabs(val) > tol){
                A(k1-offsetI, k2-offsetI) = val;              
            }
        }
    }
           
    // SF * W
    for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
        for(int j = j0; j <= J-1; ++j){
            for(int k2 = basis.rangeJ(j).firstIndex(); k2 <= basis.rangeJ(j).lastIndex(); ++k2){
                T val = a(XBSpline, j0, k1, XWavelet, j, k2);            
                if(fabs(val) > tol){
                    A(k1-offsetI,  basis.cardJ(j) + k2 - offsetJ) = val;              
                }
            }
        }
    }
    
      // W * SF
    for(int k2 = basis.mra.rangeI(j0).firstIndex(); k2 <= basis.mra.rangeI(j0).lastIndex(); ++k2){
        for(int j = j0; j <= J-1; ++j){
            for(int k1 = basis.rangeJ(j).firstIndex(); k1 <= basis.rangeJ(j).lastIndex(); ++k1){
                T val = a(XWavelet, j, k1, XBSpline, j0, k2);            
                if(fabs(val) > tol){
                    A(basis.cardJ(j) + k1 - offsetJ, k2 - offsetI) = val;              
                }
            }
        }
    }
    
        // W * W
    for(int j = j0; j <= J-1; ++j){
        for(int k1 = basis.rangeJ(j).firstIndex(); k1 <= basis.rangeJ(j).lastIndex(); ++k1){
            for(int j_ = j0; j_ <= J-1; ++j_){
                for(int k2 = basis.rangeJ(j_).firstIndex(); k2 <= basis.rangeJ(j_).lastIndex(); ++k2){
                    T val = a(XWavelet, j, k1, XWavelet, j_, k2);            
                    if(fabs(val) > tol){
                        A(basis.cardJ(j) + k1 - offsetJ, basis.cardJ(j_) + k2 - offsetJ) = val;              
                    }
                }
            }
        }
    }
    
    A.finalize();
    
    return A;
}

template<typename T, typename Basis>
template <typename RHSIntegral>
flens::DenseVector<flens::Array<T> >
Problem<T, Basis>::getRHS(RHSIntegral& rhs, int J)
{
    int j0 = basis.j0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-1;
    int offsetI = basis.mra.rangeI(j0).firstIndex()-1;
    
    int N = basis.mra.cardI(J);
    flens::DenseVector<flens::Array<T> > f(N);
    
    // SF x F
    for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
        f(k1 - offsetI) = rhs(XBSpline, j0, k1);
    }
    
    // W x F
    for(int j = j0; j <= J-1; ++j){
        for(int k1 = basis.rangeJ(j).firstIndex(); k1 <= basis.rangeJ(j).lastIndex(); ++k1){
            f(basis.cardJ(j) + k1 - offsetJ) = rhs(XWavelet, j, k1);
        }
    }
    
    return f;
}

template<typename T, typename Basis>
template <typename Preconditioner>
flens::DiagonalMatrix<T>    
Problem<T, Basis>::getPreconditioner(Preconditioner& P, int J)
{
    int j0 = basis.j0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-1;
    int offsetI = basis.mra.rangeI(j0).firstIndex()-1;
    
    int N = basis.mra.cardI(J);
    flens::DenseVector<flens::Array<T> > D(N);
    
    // SF x SF
    for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
        D(k1 - offsetI) = P(XBSpline, j0, k1);
    }
    
    // W x F
    for(int j = j0; j <= J-1; ++j){
        for(int k1 = basis.rangeJ(j).firstIndex(); k1 <= basis.rangeJ(j).lastIndex(); ++k1){
            D(basis.cardJ(j) + k1 - offsetJ) = P(XWavelet, j, k1);
        }
    }
    
    return flens::DiagonalMatrix<T>(D);
}


}