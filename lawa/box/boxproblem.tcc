/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

namespace lawa{
    
template<typename T, typename Basis>
BoxProblem<T, Basis>::
BoxProblem(Basis _basis)
    : basis(_basis)
{        
}


template<typename T, typename Basis>
template<typename BilinearForm>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
BoxProblem<T, Basis>::
getStiffnessMatrix(BilinearForm& a, int J_x, int J_y, T tol)
{   
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    BoxIndex<Basis> I(basis, J_x, J_y);
    
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A(basis.dim(J_x,J_y), basis.dim(J_x,J_y));
    
    bool spline = true;
    bool wavelet = false;
                                                 
     /* ============  v = Scaling Fct x Scaling Fct ==========================*/
     //std::cout << "===== v = SF * SF =======" << std::endl;
     Range<int> Rvx = b1.mra.rangeI(b1.j0);
     Range<int> Rvy = b2.mra.rangeI(b2.j0);
     Range<int> Rux = b1.mra.rangeI(b1.j0);
     Range<int> Ruy = b2.mra.rangeI(b2.j0);
     for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
       for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

         /* u = Scaling Fct x Scaling Fct */ 
         Rux = b1.mra.rangeI(b1.j0);
         Ruy = b2.mra.rangeI(b2.j0);    
         for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
           for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){    
                         
             T val = a(spline, b1.j0, kux, spline, b2.j0, kuy, 
                       spline, b1.j0, kvx, spline, b2.j0, kvy);
             if(fabs(val) > tol){
                 A(I(spline, b1.j0, kvx, spline, b2.j0, kvy), 
                   I(spline, b1.j0, kux, spline, b2.j0, kuy)) = val;
             }
             
           }
         }

         /* u = Scaling Fct x Wavelet */ 
         Rux = b1.mra.rangeI(b1.j0);
         for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++juy){
           Ruy = b2.rangeJ(juy); 
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
               T val = a(spline, b1.j0, kux, wavelet, juy, kuy, 
                         spline, b1.j0, kvx, spline,  b2.j0, kvy);
               if(fabs(val) > tol){
                   A(I(spline, b1.j0, kvx, spline, b2.j0, kvy),
                     I(spline, b1.j0, kux, wavelet, juy, kuy)) = val;
               }
               
             }
           }  
         }

         /* u = Wavelet x Scaling Function */ 
         Ruy = b2.mra.rangeI(b2.j0);
         for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jux){
           Rux = b1.rangeJ(jux); 
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(wavelet,  jux, kux, spline, b2.j0, kuy, 
                         spline, b1.j0, kvx, spline, b2.j0, kvy);
               if(fabs(val) > tol){
                   A(I(spline, b1.j0, kvx, spline, b2.j0, kvy),
                     I(wavelet, jux, kux, spline, b2.j0, kuy)) = val;
               }

             }
           }
         }

         /* u = Wavelet x Wavelet */ 
         for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0) - 1; ++jux){
           Rux = b1.rangeJ(jux); 
           for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, jux) - 1; ++juy){    
             Ruy = b2.rangeJ(juy); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(wavelet, jux, kux,  wavelet, juy, kuy, 
                           spline,  b1.j0, kvx, spline, b2.j0, kvy);
                 if(fabs(val) > tol){
                     A(I(spline, b1.j0, kvx, spline, b2.j0, kvy),
                       I(wavelet, jux, kux, wavelet, juy, kuy)) = val;
                 }

               }
             }
           }            
         }   
       }
     }    

     /* ============  v = Scaling Fct x Wavelet ==========================*/
    //std::cout << "===== v = SF * W =======" << std::endl;

     Rvx = b1.mra.rangeI(b1.j0);
     for(int jvy = b2.j0; jvy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; jvy++){
       Rvy = b2.rangeJ(jvy);
       for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
         for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

           /* u = Scaling Fct x Scaling Fct */ 
           Rux = b1.mra.rangeI(b1.j0);
           Ruy = b2.mra.rangeI(b2.j0);    
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(spline, b1.j0, kux, spline,  b2.j0, kuy, 
                         spline, b1.j0, kvx, wavelet, jvy, kvy);
               if(fabs(val) > tol){
                   A(I(spline, b1.j0, kvx, wavelet, jvy, kvy),
                     I(spline, b1.j0, kux, spline, b2.j0, kuy)) = val;
               }

             }
           }

           /* u = Scaling Fct x Wavelet */ 
           Rux = b1.mra.rangeI(b1.j0);
           for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++juy){
             Ruy = b2.rangeJ(juy); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                
                 T val = a(spline, b1.j0, kux, wavelet, juy, kuy, 
                           spline, b1.j0, kvx, wavelet, jvy, kvy);
                 if(fabs(val) > tol){  
                     A(I(spline, b1.j0, kvx, wavelet, jvy, kvy),
                       I(spline, b1.j0, kux, wavelet, juy, kuy)) = val;
                 }                
                
               }
             }  
           }

           /* u = Wavelet x Scaling Function */ 
           Ruy = b2.mra.rangeI(b2.j0);
           for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jux){
             Rux = b1.rangeJ(jux); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                 T val = a(wavelet, jux, kux,   spline, b2.j0, kuy, 
                           spline,  b1.j0, kvx, wavelet, jvy, kvy);
                 if(fabs(val) > tol){
                     A(I(spline, b1.j0, kvx, wavelet, jvy, kvy), 
                       I(wavelet, jux, kux, spline, b2.j0, kuy)) = val;
                 }
                 
               }
             }
           }

           /* u = Wavelet x Wavelet */ 
           for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jux){
             Rux = b1.rangeJ(jux);
             for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, jux) -1; ++juy){
               Ruy = b2.rangeJ(juy);
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                   T val = a(wavelet, jux, kux,   wavelet, juy, kuy, 
                             spline,  b1.j0, kvx, wavelet, jvy, kvy);
                   if(fabs(val) > tol){
                       A(I(spline, b1.j0, kvx, wavelet, jvy, kvy), 
                         I(wavelet, jux, kux, wavelet, juy, kuy)) = val;
                   }

                 }
               } 
             }
           }
         }
       }
     }

     /* ============  v = Wavelet x Scaling Fct ==========================*/
     //cout << "===== v = W * SF =======" << endl;

     Rvy = b2.mra.rangeI(b2.j0);
     for(int jvx = b1.j0; jvx <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; jvx++){
       Rvx = b1.rangeJ(jvx);
       for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
         for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

           /* u = Scaling Fct x Scaling Fct */ 
           Rux = b1.mra.rangeI(b1.j0);
           Ruy = b2.mra.rangeI(b2.j0);  
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(spline,  b1.j0, kux, spline, b2.j0, kuy, 
                         wavelet, jvx, kvx,   spline, b2.j0, kvy);
               if(fabs(val) > tol){
                   A(I(wavelet, jvx, kvx, spline, b2.j0, kvy),
                     I(spline, b1.j0, kux, spline, b2.j0, kuy)) = val;
               }
               
             }
           }

           /* u = Scaling Fct x Wavelet */ 
           Rux = b1.mra.rangeI(b1.j0);
           for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++juy){
             Ruy = b2.rangeJ(juy); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(spline,  b1.j0, kux, wavelet, juy, kuy, 
                           wavelet, jvx, kvx,   spline, b2.j0, kvy);
                 if(fabs(val) > tol){
                     A(I(wavelet, jvx, kvx, spline, b2.j0, kvy),
                       I(spline, b1.j0, kux, wavelet, juy, kuy)) = val;
                 }

               }
             }  
           }

           /* u = Wavelet x Scaling Function */ 
           Ruy = b2.mra.rangeI(b2.j0);
           for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jux){
             Rux = b1.rangeJ(jux); 
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(wavelet, jux, kux, spline, b2.j0, kuy, 
                           wavelet, jvx, kvx, spline, b2.j0, kvy);
                 if(fabs(val) > tol){
                     A(I(wavelet, jvx, kvx, spline, b2.j0, kvy),
                       I(wavelet, jux, kux, spline, b2.j0, kuy)) = val;
                 }
                 
               }
             }
           }

           /* u = Wavelet x Wavelet */ 
           for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jux){
             Rux = b1.rangeJ(jux);
             for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, jux) -1; ++juy){
               Ruy = b2.rangeJ(juy);
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                   T val = a(wavelet, jux, kux, wavelet, juy, kuy, 
                             wavelet, jvx, kvx, spline,  b2.j0, kvy);
                   if(fabs(val) > tol){
                       A(I(wavelet, jvx, kvx, spline, b2.j0, kvy),
                         I(wavelet, jux, kux, wavelet, juy, kuy)) = val;
                   }

                 }
               } 
             }
           }
         }
       }
     }   

     /* ============  v = Wavelet x Wavelet ==========================*/
     //cout << "===== v = W * W =======" << endl;  
     for(int jvx = b1.j0; jvx <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jvx){
       Rvx = b1.rangeJ(jvx);
       for(int jvy = b2.j0; jvy <= basis.J2_max(J_x, J_y, jvx) -1; ++jvy){
         Rvy = b2.rangeJ(jvy);
         for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
           for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

             /* u = Scaling Fct x Scaling Fct */  
             Rux = b1.mra.rangeI(b1.j0);
             Ruy = b2.mra.rangeI(b2.j0);
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                 T val = a(spline,  b1.j0, kux, spline,  b2.j0, kuy, 
                           wavelet, jvx, kvx,   wavelet, jvy, kvy);
                 if(fabs(val) > tol){
                     A(I(wavelet, jvx, kvx, wavelet, jvy, kvy),
                       I(spline,  b1.j0, kux, spline, b2.j0, kuy)) = val;
                 }

               }
             }

             /* u = Scaling Fct x Wavelet */ 
             Rux = b1.mra.rangeI(b1.j0);
             for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++juy){
               Ruy = b2.rangeJ(juy); 
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                   T val = a(spline, b1.j0, kux, wavelet, juy, kuy, 
                             wavelet, jvx, kvx,  wavelet, jvy, kvy);
                   if(fabs(val) > tol){
                       A(I(wavelet, jvx, kvx, wavelet, jvy, kvy),
                         I(spline, b1.j0, kux, wavelet, juy, kuy)) = val;
                   }

                 }
               }  
             }

             /* u = Wavelet x Scaling Function */ 
             Ruy = b2.mra.rangeI(b2.j0);
             for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jux){
               Rux = b1.rangeJ(jux); 
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                   T val = a(wavelet, jux, kux, spline,  b2.j0, kuy, 
                             wavelet, jvx, kvx, wavelet, jvy, kvy);
                   if(fabs(val) > tol){
                       A(I(wavelet, jvx, kvx, wavelet, jvy, kvy),
                         I(wavelet, jux, kux, spline, b2.j0, kuy)) = val;
                   }  

                 }
               }
             }

             /* u = Wavelet x Wavelet */ 
             for(int jux = b1.j0; jux <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jux){
               Rux = b1.rangeJ(jux);
               for(int juy = b2.j0; juy <= basis.J2_max(J_x, J_y, jux) -1; ++juy){
                 Ruy = b2.rangeJ(juy);
                 for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                   for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                     T val = a(wavelet, jux, kux, wavelet, juy, kuy, 
                               wavelet, jvx, kvx, wavelet, jvy, kvy);
                     if(fabs(val) > tol){
                         A(I(wavelet, jvx, kvx, wavelet, jvy, kvy),
                           I(wavelet, jux, kux, wavelet, juy, kuy)) = val;
                     }

                   }
                 } 
               }
             }           
           }
         } 
       }
     }
    
    A.finalize();
    return A;
 
}  


template<typename T, typename Basis>
template<typename RHSIntegral>
flens::DenseVector<flens::Array<T> >
BoxProblem<T, Basis>::getRHS(RHSIntegral& rhs, int J_x, int J_y)
{
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    BoxIndex<Basis> I(basis, J_x, J_y);
    
    bool spline = true;
    bool wavelet = false;
     
    flens::DenseVector<flens::Array<T> > f(basis.dim(J_x, J_y));

    /*  ============  v = Scaling Fct x Scaling Fct ==========================*/
    //std::cout << "SF x SF : " << std::endl;
    Range<int> Rvx = b1.mra.rangeI(b1.j0);
    Range<int> Rvy = b2.mra.rangeI(b2.j0);
    for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
      for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
          
          f(I(spline, b1.j0, kvx, spline, b2.j0, kvy)) 
              = rhs(spline, b1.j0, kvx, spline, b2.j0, kvy);

      }
    }

    /* ============  v = Scaling Fct x Wavelet ==========================*/
    //std::cout << "SF x W : " << std::endl;
    Rvx = b1.mra.rangeI(b1.j0);
    for(int jvy = b2.j0; jvy <= basis.J2_max(J_x, J_y, b1.j0-1) - 1; ++jvy){
      Rvy = b2.rangeJ(jvy);
      for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
        for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
        
          f(I(spline, b1.j0, kvx, wavelet, jvy, kvy)) 
              = rhs(spline, b1.j0, kvx, wavelet, jvy, kvy);
                     
        }
      }  
    }

    /* ============  v = Wavelet x Scaling Fct ==========================*/
    //std::cout << "W x SF : " << std::endl;
    Rvy = b2.mra.rangeI(b2.j0);
    for(int jvx = b1.j0; jvx <= basis.J1_max(J_x, J_y, b2.j0-1) - 1; ++jvx){
      Rvx = b1.rangeJ(jvx);
      for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
        for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
            
          f(I(wavelet, jvx, kvx, spline, b2.j0, kvy)) 
              = rhs(wavelet, jvx, kvx, spline, b2.j0, kvy);
          
        }
      }  
    }

    /*  ============  v = Wavelet x Wavelet ==========================*/
    //std::cout << "W x W : " << std::endl;
    for(int jvx = b1.j0; jvx <= basis.J1_max(J_x, J_y, b2.j0) -1; ++jvx){
      Rvx = b1.rangeJ(jvx);
      for(int jvy = b2.j0; jvy <= basis.J2_max(J_x, J_y, jvx) - 1; ++jvy){
        Rvy = b2.rangeJ(jvy);
        for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
          for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
            
            f(I(wavelet, jvx, kvx, wavelet, jvy, kvy)) 
                = rhs(wavelet, jvx, kvx, wavelet, jvy, kvy);

          }
        }
      }
    }
    
    return f;
}


template<typename T, typename Basis>
template<typename Preconditioner>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    
BoxProblem<T, Basis>::
getPreconditioner(Preconditioner& P, int J_x, int J_y)
{
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    BoxIndex<Basis> I(basis, J_x, J_y);
    
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > D(basis.dim(J_x, J_y), basis.dim(J_x, J_y));
    
    bool spline = true;
    bool wavelet = false;

    /* SF x SF */
    Range<int> Rx = b1.mra.rangeI(b1.j0);
    Range<int> Ry = b2.mra.rangeI(b2.j0);   
    for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
        for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
            
            int index = I(spline, b1.j0, kx, spline, b2.j0, ky);
            D(index, index) = P(spline, b1.j0, kx, spline, b2.j0, ky);
            
        }
    }
    
    /* SF x W */
    Rx = b1.mra.rangeI(b1.j0);
    for(int jy = b2.j0; jy <= basis.J2_max(J_x, J_y, b1.j0-1)-1; ++jy){
        Ry = b2.rangeJ(jy);
        for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
            for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                
                int index = I(spline, b1.j0, kx, wavelet, jy, ky);
                D(index, index) = P(spline, b1.j0, kx, wavelet, jy, ky);
                
            }
        }
    }
    
    /* W x SF */
    Ry = b2.mra.rangeI(b2.j0);
    for(int jx = b1.j0; jx <= basis.J1_max(J_x, J_y, b2.j0-1) -1; ++jx){
        Rx = b1.rangeJ(jx);
        for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
            for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                
                int index = I(wavelet, jx, kx, spline, b2.j0, ky);
                D(index, index) = P(wavelet, jx, kx, spline, b2.j0, ky);
            
            }
        }      
    }
    
    /* W x W */
    for(int jx = b1.j0; jx <= basis.J1_max(J_x, J_y, b2.j0) - 1; ++jx){
        Rx = b1.rangeJ(jx);
        for(int jy = b2.j0; jy <= basis.J2_max(J_x, J_y, jx) - 1; ++jy){
            Ry = b2.rangeJ(jy);
            for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
                for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                    
                    int index = I(wavelet, jx, kx, wavelet, jy, ky);
                    D(index, index) = P(wavelet, jx, kx, wavelet, jy, ky);
                
                }
            } 
        }
    }
    
    return D;
    
}
    
} // namespace lawa