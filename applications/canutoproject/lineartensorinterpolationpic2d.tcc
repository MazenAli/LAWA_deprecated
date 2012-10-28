namespace lawa {

template <typename T>
int
LinearTensorInterpolationPic2D<T>::N1;

template <typename T>
int
LinearTensorInterpolationPic2D<T>::N2;

template <typename T>
flens::DenseVector<Array<T> >
LinearTensorInterpolationPic2D<T>::sing_pts_x;

template <typename T>
flens::DenseVector<Array<T> >
LinearTensorInterpolationPic2D<T>::sing_pts_y;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
LinearTensorInterpolationPic2D<T>::coeffs;

template <typename T>
T
LinearTensorInterpolationPic2D<T>::diffusion_coeff;

template <typename T>
Basis<T,Primal,R,CDF>
LinearTensorInterpolationPic2D<T>::linearTensorInterpolBasis(2,2,0);

template <typename T>
void
LinearTensorInterpolationPic2D<T>::readPicture(const char* filename)
{
    std::ifstream infile (filename);
    int numOfRows = 0, numOfCols = 0;
    if (infile.is_open()) {
        std::string line;
        while (std::getline( infile, line, '\n' ) ) {
            ++numOfRows;
            std::istringstream line_ss(line);
            if (numOfRows==1) {
                std::string matrix_entry_string;
                while(std::getline( line_ss, matrix_entry_string, ',' )) {
                    ++numOfCols;
                }
            }
        }
        infile.close();
    }

    diffusion_coeff = 1.;
    N1 = numOfRows-1;
    N2 = numOfCols-1;
    std::cerr << "N1 = " << N1 << ", N2 = " << N2 << std::endl;
    coeffs.engine().resize(_(0,N1),_(0,N2));

    sing_pts_x.engine().resize(N1+1);
    sing_pts_y.engine().resize(N2+1);

    for (int i=1; i<=N1+1; ++i) {
        sing_pts_x(i) = (T)(i-1)/(T)N1;
    }
    for (int i=1; i<=N2+1; ++i) {
        sing_pts_y(i) = (T)(i-1)/(T)N2;
    }
    //std::cout << "sing_pts_x = " << sing_pts_x << std::endl;
    //std::cout << "sing_pts_y = " << sing_pts_y << std::endl;

    std::cout << "N1 = " << N1 << ", h1 = " << 1./N1 << std::endl;
    std::cout << "N2 = " << N2 << ", h2 = " << 1./N2 << std::endl;

    std::ifstream infile2 (filename);
    int k1=0, k2=0;
    if (infile2.is_open()) {
        std::cout << "File is open, ready to read..." << std::endl;
        std::string line;
        while (std::getline( infile2, line, '\n' ) ) {
            std::istringstream line_ss(line);
            std::string matrix_entry_string;
            k2 = 0;
            while(std::getline( line_ss, matrix_entry_string, ',' )) {
                T entry = atof(matrix_entry_string.c_str());
                coeffs(k1,k2) = entry;
                ++k2;
            }
            ++k1;
        }
        infile2.close();
    }
}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::evaluateInterpolation(T x1, T x2, int deriv_x1, int deriv_x2) {

    T y1 = N1*x1, y2 = N2*x2;
    int k11 = std::floor(y1), k12 = std::ceil(y1);
    int k21 = std::floor(y2), k22 = std::ceil(y2);

    T ret = 0.;
    if ((k11==k12) && (k21==k22)) {
        ret =  coeffs(k11,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2);
    }
    else if ((k11!=k12) && (k21==k22)) {
        ret =   coeffs(k11,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2)
              + coeffs(k12,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k12,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2);
    }
    else if ((k11==k12) && (k21!=k22)) {
        ret =    coeffs(k11,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2)
               + coeffs(k11,k22) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k22,deriv_x2);
    }
    else {
        ret =    coeffs(k11,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2)
               + coeffs(k11,k22) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k22,deriv_x2)
               + coeffs(k12,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k12,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2)
               + coeffs(k12,k22) * linearTensorInterpolBasis.mra.phi(y1,0,k12,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k22,deriv_x2);
    }

    T factor1 = deriv_x1==1 ? N1 : 1;
    T factor2 = deriv_x2==1 ? N2 : 1;
    return factor1 * factor2 * ret;

}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::evaluateInterpolation(T x1, T x2) {
    return LinearTensorInterpolationPic2D<T>::evaluateInterpolation(x1,x2,0,0);
}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::dx1_evaluateInterpolation(T x1, T x2) {
    return diffusion_coeff*LinearTensorInterpolationPic2D<T>::evaluateInterpolation(x1,x2,1,0);
}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::dx2_evaluateInterpolation(T x1, T x2) {
    return diffusion_coeff*LinearTensorInterpolationPic2D<T>::evaluateInterpolation(x1,x2,0,1);
}


template <typename T>
T
LinearTensorInterpolationPic2D<T>::evaluateLiftingFunction(T x1, T x2, int deriv_x1, int deriv_x2) {

    T y1 = N1*x1, y2 = N2*x2;

    int k1_firstIndex = 0, k1_lastIndex = N1;
    int k2_firstIndex = 0, k2_lastIndex = N2;
    int k11 = std::floor(y1), k12 = std::ceil(y1);
    int k21 = std::floor(y2), k22 = std::ceil(y2);

    int corner = 0, edge = 0;

    if (0<=x1 && x1<=sing_pts_x(1)) {
        if (0<=x2 && x2<=sing_pts_y(1))                         corner = 1;
        else if (sing_pts_y(1)<x2  && x2<=sing_pts_y(N2))       edge   = 4;
        else if (sing_pts_y(N2)<=x2 && x2<=sing_pts_y(N2+1))    corner = 4;
        else return 0;
    }
    else if (sing_pts_x(1)<x1 && x1<sing_pts_x(N1)) {
        if (0<=x2 && x2<=sing_pts_y(1))                         edge   = 1;
        else if (sing_pts_y(N2)<=x2 && x2<=sing_pts_y(N2+1))    edge   = 3;
        else return 0.;
    }
    else if (sing_pts_x(N1)<=x1 && x1<=sing_pts_x(N1+1)) {
        if (0<=x2 && x2<=sing_pts_y(1))                         corner = 2;
        else if (sing_pts_y(1)<=x2 && x2<sing_pts_y(N2))        edge   = 2;
        else if (sing_pts_y(N2)<=x2 && x2<=sing_pts_y(N2+1))    corner = 3;
        else return 0.;
    }
    else return 0.;

    //std::cerr << "corner = " << corner << ", edge = " << edge << std::endl;

    T ret = 0.;

    if (corner>0) {
        if (corner==1) {
            ret =    coeffs(k1_firstIndex,k2_firstIndex)   * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex,deriv_x1)   * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex,deriv_x2)
                   + coeffs(k1_firstIndex,k2_firstIndex+1) * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex,deriv_x1)   * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex+1,deriv_x2)
                   + coeffs(k1_firstIndex+1,k2_firstIndex) * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex+1,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex,deriv_x2);
        }
        else if (corner==2) {
            ret =    coeffs(k1_lastIndex-1,k2_firstIndex)  * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex-1,deriv_x1)  * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex,deriv_x2)
                   + coeffs(k1_lastIndex,k2_firstIndex)    * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex,deriv_x1)    * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex,deriv_x2)
                   + coeffs(k1_lastIndex,k2_firstIndex+1)  * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex,deriv_x1)    * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex+1,deriv_x2);
        }
        else if (corner==3) {
            ret =    coeffs(k1_lastIndex,k2_lastIndex-1)   * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex,deriv_x1)    * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex-1,deriv_x2)
                   + coeffs(k1_lastIndex,k2_lastIndex)     * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex,deriv_x1)    * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex,deriv_x2)
                   + coeffs(k1_lastIndex-1,k2_lastIndex)   * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex-1,deriv_x1)  * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex,deriv_x2);
        }
        else if (corner==4) {
            ret =    coeffs(k1_firstIndex+1,k2_lastIndex)  * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex+1,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex,deriv_x2)
                   + coeffs(k1_firstIndex,k2_lastIndex)    * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex,deriv_x1)   * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex,deriv_x2)
                   + coeffs(k1_firstIndex,k2_lastIndex-1)  * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex,deriv_x1)   * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex-1,deriv_x2);
        }
    }
    else {
        if (edge==1) {
            if (k11!=k12) ret =    coeffs(k11,k2_firstIndex) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex,deriv_x2)
                                 + coeffs(k12,k2_firstIndex) * linearTensorInterpolBasis.mra.phi(y1,0,k12,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex,deriv_x2);
            else          ret =    coeffs(k11,k2_firstIndex) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k2_firstIndex,deriv_x2);
        }
        else if (edge==2) {
            if (k21!=k22) ret =    coeffs(k1_lastIndex,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2)
                                 + coeffs(k1_lastIndex,k22) * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k22,deriv_x2);
            else          ret =    coeffs(k1_lastIndex,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k1_lastIndex,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2);
        }
        else if (edge==3) {
            if (k11!=k12) ret =    coeffs(k11,k2_lastIndex) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex,deriv_x2)
                                 + coeffs(k12,k2_lastIndex) * linearTensorInterpolBasis.mra.phi(y1,0,k12,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex,deriv_x2);
            else          ret =    coeffs(k11,k2_lastIndex) * linearTensorInterpolBasis.mra.phi(y1,0,k11,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k2_lastIndex,deriv_x2);
        }
        else if (edge==4) {
            if (k21!=k22) ret =    coeffs(k1_firstIndex,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2)
                                 + coeffs(k1_firstIndex,k22) * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k22,deriv_x2);
            else          ret =    coeffs(k1_firstIndex,k21) * linearTensorInterpolBasis.mra.phi(y1,0,k1_firstIndex,deriv_x1) * linearTensorInterpolBasis.mra.phi(y2,0,k21,deriv_x2);
        }
    }

    T factor1 = deriv_x1==1 ? N1 : 1;
    T factor2 = deriv_x2==1 ? N2 : 1;
    return factor1 * factor2 * ret;

}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::evaluateLiftingFunction(T x1, T x2) {
    return LinearTensorInterpolationPic2D<T>::evaluateLiftingFunction(x1, x2, 0, 0);
}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::dx1_evaluateLiftingFunction(T x1, T x2) {
    return diffusion_coeff*LinearTensorInterpolationPic2D<T>::evaluateLiftingFunction(x1, x2, 1, 0);
}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::dx2_evaluateLiftingFunction(T x1, T x2) {
    return diffusion_coeff*LinearTensorInterpolationPic2D<T>::evaluateLiftingFunction(x1, x2, 0, 1);
}


template <typename T>
T
LinearTensorInterpolationPic2D<T>::evaluateInterpolationMinusLiftingFunction(T x1, T x2) {
    return   LinearTensorInterpolationPic2D<T>::evaluateInterpolation(x1,x2,0,0)
           - LinearTensorInterpolationPic2D<T>::evaluateLiftingFunction(x1,x2,0,0);
}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::dx1_evaluateInterpolationMinusLiftingFunction(T x1, T x2) {
    return   diffusion_coeff*(LinearTensorInterpolationPic2D<T>::evaluateInterpolation(x1,x2,1,0)
           - LinearTensorInterpolationPic2D<T>::evaluateLiftingFunction(x1,x2,1,0) );
}

template <typename T>
T
LinearTensorInterpolationPic2D<T>::dx2_evaluateInterpolationMinusLiftingFunction(T x1, T x2) {
    return   diffusion_coeff*(LinearTensorInterpolationPic2D<T>::evaluateInterpolation(x1,x2,0,1)
           - LinearTensorInterpolationPic2D<T>::evaluateLiftingFunction(x1,x2,0,1) );
}

template <typename T>
void
LinearTensorInterpolationPic2D<T>::plotInterpolation(const char* filename, T h1, T h2) {
    std::ofstream plotfile(filename);
    for (T x1=0.; x1<=1.; x1+=h1) {
        for (T x2=0.; x2<=1.; x2+=h2) {
            plotfile << x1 << " " << x2 << " " << LinearTensorInterpolationPic2D<T>::evaluateInterpolation(x1, x2) << std::endl;
        }
        plotfile << std::endl;
    }
    plotfile.close();
}

}   // namespace lawa
