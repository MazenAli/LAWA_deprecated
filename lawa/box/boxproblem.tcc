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
    
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral, typename Preconditioner>
BoxProblem<T, Basis, BilinearForm, RHSIntegral, Preconditioner>::BoxProblem(Basis _basis, BilinearForm _a, RHSIntegral _rhs, Preconditioner _P)
    : basis(_basis), a(_a), rhs(_rhs), P(_P)
{        
}


template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral, typename Preconditioner>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
BoxProblem<T, Basis, BilinearForm, RHSIntegral, Preconditioner>::getStiffnessMatrix(int J_x, int J_y, T tol)
{   
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    flens::DenseVector<flens::Array<T> >  N(2);
    N = b1.mra.cardI(J_x), b2.mra.cardI(J_y);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A(N(1)*N(2), N(1)*N(2));
    
    int osJx = b1.rangeJ(b1.j0).firstIndex() - 1;
    int osJy = b2.rangeJ(b2.j0).firstIndex() - 1;
    int osIx = b1.mra.rangeI(b1.j0).firstIndex() - 1;
    int osIy = b2.mra.rangeI(b2.j0).firstIndex() - 1;
    
    bool spline = true;
    bool wavelet = false;
                                                 
     /* ============  v = Scaling Fct x Scaling Fct ==========================*/
     //cout << "===== v = SF * SF =======" << endl;
     Range<int> Rvx = b1.mra.rangeI(b1.j0);
     Range<int> Rvy = b2.mra.rangeI(b2.j0);
     Range<int> Rux = b1.mra.rangeI(b1.j0);
     Range<int> Ruy = b2.mra.rangeI(b2.j0);
     int Cvx = b1.mra.cardI(b1.j0);
     int Cvy = b2.mra.cardI(b2.j0);
     int Cux, Cuy;       
     for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
       for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

         /* u = Scaling Fct x Scaling Fct */ 
         Rux = b1.mra.rangeI(b1.j0);
         Ruy = b2.mra.rangeI(b2.j0); 
         Cux = b1.mra.cardI(b1.j0);
         Cuy = b2.mra.cardI(b2.j0);    
         for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
           for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){    
                         
             T val = a(spline, spline, b1.j0, kux, b2.j0, kuy, spline, spline, b1.j0, kvx, b2.j0, kvy);
             if(fabs(val) > tol){
                 A((kvx-osIx-1)*Cvy + kvy-osIy, 
                   (kux-osIx-1)*Cuy + kuy-osIy) = val;
             }
             
           }
         }

         /* u = Scaling Fct x Wavelet */ 
         Rux = b1.mra.rangeI(b1.j0);
         Cux = b1.mra.cardI(b1.j0);
         for(int juy = b2.j0; juy <= J_y - 1; ++juy){
           Ruy = b2.rangeJ(juy); 
           Cuy = b2.mra.cardI(juy);    
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
               T val = a(spline, wavelet, b1.j0, kux, juy, kuy, spline, spline, b1.j0, kvx, b2.j0, kvy);
               if(fabs(val) > tol){
                   A((kvx-osIx-1)*Cvy + kvy-osIy, 
                     Cux*Cuy + (kux-osIx-1)*Cuy + kuy-osJy) = val;
               }
               
             }
           }  
         }

         /* u = Wavelet x Scaling Function */ 
         Ruy = b2.mra.rangeI(b2.j0);
         Cuy = b2.mra.cardI(b2.j0);
         for(int jux = b1.j0; jux <= J_x - 1; ++jux){
           Rux = b1.rangeJ(jux); 
           Cux = b1.mra.cardI(jux);    
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(wavelet, spline, jux, kux, b2.j0, kuy, spline, spline, b1.j0, kvx, b2.j0, kvy);
               if(fabs(val) > tol){
                   A((kvx-osIx-1)*Cvy + kvy-osIy, 
                     Cux*b2.mra.cardI(J_y) + (kux-osJx-1)*Cuy + kuy-osIy) = val;
               }

             }
           }
         }

         /* u = Wavelet x Wavelet */ 
         for(int jux = b1.j0; jux <= J_x - 1; ++jux){
           Rux = b1.rangeJ(jux); 
           Cux = b1.mra.cardI(jux);
           for(int juy = b2.j0; juy <= J_y - 1; ++juy){    
             Ruy = b2.rangeJ(juy); 
             Cuy = b2.mra.cardI(juy);
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(wavelet, wavelet, jux, kux, juy, kuy, spline, spline, b1.j0, kvx, b2.j0, kvy);
                 if(fabs(val) > tol){
                     A((kvx-osIx-1)*Cvy + kvy-osIy, 
                       Cux*b2.mra.cardI(J_y) + Cux*Cuy + (kux-osJx-1)*Cuy + kuy-osJy) = val;
                 }

               }
             }
           }            
         }   
       }
     }    

     /* ============  v = Scaling Fct x Wavelet ==========================*/
    // cout << "===== v = SF * W =======" << endl;

     Rvx = b1.mra.rangeI(b1.j0);
     Cvx = b1.mra.cardI(b1.j0);    
     for(int jvy = b2.j0; jvy <= J_y - 1; jvy++){
       Rvy = b2.rangeJ(jvy);
       Cvy = b2.mra.cardI(jvy);

       for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
         for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

           /* u = Scaling Fct x Scaling Fct */ 
           Rux = b1.mra.rangeI(b1.j0);
           Ruy = b2.mra.rangeI(b2.j0); 
           Cux = b1.mra.cardI(b1.j0);
           Cuy = b2.mra.cardI(b2.j0);    
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(spline, spline, b1.j0, kux, b2.j0, kuy, spline, wavelet, b1.j0, kvx, jvy, kvy);
               if(fabs(val) > tol){
                   A(Cvx*Cvy + (kvx-osIx-1)*Cvy + kvy-osJy, 
                     (kux - osIx - 1)*Cuy + kuy-osIy) = val;
               }

             }
           }

           /* u = Scaling Fct x Wavelet */ 
           Rux = b1.mra.rangeI(b1.j0);
           Cux = b1.mra.cardI(b1.j0);
           for(int juy = b2.j0; juy <= J_y - 1; ++juy){
             Ruy = b2.rangeJ(juy); 
             Cuy = b2.mra.cardI(juy);    
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                
                 T val = a(spline, wavelet, b1.j0, kux, juy, kuy, spline, wavelet, b1.j0, kvx, jvy, kvy);
                 if(fabs(val) > tol){
                     A(Cvx*Cvy + (kvx-osIx-1)*Cvy + kvy-osJy, 
                       Cux*Cuy + (kux-osIx-1)*Cuy + kuy-osJy) = val;
                 }                
                
               }
             }  
           }

           /* u = Wavelet x Scaling Function */ 
           Ruy = b2.mra.rangeI(b2.j0);
           Cuy = b2.mra.cardI(b2.j0);
           for(int jux = b1.j0; jux <= J_x - 1; ++jux){
             Rux = b1.rangeJ(jux); 
             Cux = b1.mra.cardI(jux);    
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                 T val = a(wavelet, spline, jux, kux, b2.j0, kuy, spline, wavelet, b1.j0, kvx, jvy, kvy);
                 if(fabs(val) > tol){
                     A(Cvx*Cvy + (kvx-osIx-1)*Cvy + kvy-osJy, 
                       Cux*b2.mra.cardI(J_y) + (kux-osJx-1)*Cuy + kuy-osIy) = val;
                 }
                 
               }
             }
           }

           /* u = Wavelet x Wavelet */ 
           for(int jux = b1.j0; jux <= J_x -1; ++jux){
             Rux = b1.rangeJ(jux);
             Cux = b1.mra.cardI(jux);
             for(int juy = b2.j0; juy <= J_y -1; ++juy){
               Ruy = b2.rangeJ(juy);
               Cuy = b2.mra.cardI(juy);
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                   T val = a(wavelet, wavelet, jux, kux, juy, kuy, spline, wavelet, b1.j0, kvx, jvy, kvy);
                   if(fabs(val) > tol){
                       A(Cvx*Cvy + (kvx-osIx-1)*Cvy + kvy-osJy, 
                         Cux*b2.mra.cardI(J_y) + Cux*Cuy + (kux-osJx-1)*Cuy + kuy-osIy) = val;
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
     Cvy = b2.mra.cardI(b2.j0);    
     for(int jvx = b1.j0; jvx <= J_x - 1; jvx++){
       Rvx = b1.rangeJ(jvx);
       Cvx = b1.mra.cardI(jvx);

       for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
         for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

           /* u = Scaling Fct x Scaling Fct */ 
           Rux = b1.mra.rangeI(b1.j0);
           Ruy = b2.mra.rangeI(b2.j0); 
           Cux = b1.mra.cardI(b1.j0);
           Cuy = b2.mra.cardI(b2.j0);    
           for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
             for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
               
               T val = a(spline, spline, b1.j0, kux, b2.j0, kuy, wavelet, spline, jvx, kvx, b2.j0, kvy);
               if(fabs(val) > tol){
                   A(Cvx*b2.mra.cardI(J_y) + (kvx-osJx-1)*Cvy + kvy-osIy,
                     (kux-osJx-1)*Cuy + kuy-osIy) = val;
               }
               
             }
           }

           /* u = Scaling Fct x Wavelet */ 
           Rux = b1.mra.rangeI(b1.j0);
           Cux = b1.mra.cardI(b1.j0);
           for(int juy = b2.j0; juy <= J_y - 1; ++juy){
             Ruy = b2.rangeJ(juy); 
             Cuy = b2.mra.cardI(juy);    
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(spline, wavelet, b1.j0, kux, juy, kuy, wavelet, spline, jvx, kvx, b2.j0, kvy);
                 if(fabs(val) > tol){
                     A(Cvx*b2.mra.cardI(J_y) + (kvx-osJx-1)*Cvy + kvy-osIy,
                       Cux*Cuy + (kux-osIx-1)*Cuy + kuy-osJy) = val;
                 }

               }
             }  
           }

           /* u = Wavelet x Scaling Function */ 
           Ruy = b2.mra.rangeI(b2.j0);
           Cuy = b2.mra.cardI(b2.j0);
           for(int jux = b1.j0; jux <= J_x - 1; ++jux){
             Rux = b1.rangeJ(jux); 
             Cux = b1.mra.cardI(jux);    
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                 T val = a(wavelet, spline, jux, kux, b2.j0, kuy, wavelet, spline, jvx, kvx, b2.j0, kvy);
                 if(fabs(val) > tol){
                     A(Cvx*b2.mra.cardI(J_y) + (kvx-osJx-1)*Cvy + kvy-osIy,
                       Cux*b2.mra.cardI(J_y) + (kux-osJx-1)*Cuy + kuy-osIy) = val;
                 }
                 
               }
             }
           }

           /* u = Wavelet x Wavelet */ 
           for(int jux = b1.j0; jux <= J_x -1; ++jux){
             Rux = b1.rangeJ(jux);
             Cux = b1.mra.cardI(jux);
             for(int juy = b2.j0; juy <= J_y -1; ++juy){
               Ruy = b2.rangeJ(juy);
               Cuy = b2.mra.cardI(juy);
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                   T val = a(wavelet, wavelet, jux, kux, juy, kuy, wavelet, spline, jvx, kvx, b2.j0, kvy);
                   if(fabs(val) > tol){
                       A(Cvx*b2.mra.cardI(J_y) + (kvx-osJx-1)*Cvy + kvy-osIy,
                         Cux*b2.mra.cardI(J_y) + Cux*Cuy + (kux-osJx-1)*Cuy + kuy-osJy) = val;
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
     for(int jvx = b1.j0; jvx <= J_x -1; ++jvx){
       Rvx = b1.rangeJ(jvx);
       Cvx = b1.mra.cardI(jvx);
       for(int jvy = b2.j0; jvy <= J_y -1; ++jvy){
         Rvy = b2.rangeJ(jvy);
         Cvy = b2.mra.cardI(jvy);
         for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
           for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){

             /* u = Scaling Fct x Scaling Fct */  
             Rux = b1.mra.rangeI(b1.j0);
             Ruy = b2.mra.rangeI(b2.j0);
             Cux = b1.mra.cardI(b1.j0);
             Cuy = b2.mra.cardI(b2.j0);  
             for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
               for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                 T val = a(spline, spline, b1.j0, kux, b2.j0, kuy, wavelet, wavelet, jvx, kvx, jvy, kvy);
                 if(fabs(val) > tol){
                     A(Cvx*b2.mra.cardI(J_y) + Cvx*Cvy + (kvx-osJx-1)*Cvy + kvy-osJy,
                       (kux-osJx-1)*Cuy + kuy-osIy) = val;
                 }

               }
             }

             /* u = Scaling Fct x Wavelet */ 
             Rux = b1.mra.rangeI(b1.j0);
             Cux = b1.mra.cardI(b1.j0);
             for(int juy = b2.j0; juy <= J_y - 1; ++juy){
               Ruy = b2.rangeJ(juy); 
               Cuy = b2.mra.cardI(juy);    
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                   
                   T val = a(spline, wavelet, b1.j0, kux, juy, kuy, wavelet, wavelet, jvx, kvx, jvy, kvy);
                   if(fabs(val) > tol){
                       A(Cvx*b2.mra.cardI(J_y) + Cvx*Cvy + (kvx-osJx-1)*Cvy + kvy-osJy,
                         Cux*Cuy + (kux-osIx-1)*Cuy + kuy-osJy) = val;
                   }

                 }
               }  
             }

             /* u = Wavelet x Scaling Function */ 
             Ruy = b2.mra.rangeI(b2.j0);
             Cuy = b2.mra.cardI(b2.j0);
             for(int jux = b1.j0; jux <= J_x - 1; ++jux){
               Rux = b1.rangeJ(jux); 
               Cux = b1.mra.cardI(jux);    
               for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                 for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                 
                   T val = a(wavelet, spline, jux, kux, b2.j0, kuy, wavelet, wavelet, jvx, kvx, jvy, kvy);
                   if(fabs(val) > tol){
                       A(Cvx*b2.mra.cardI(J_y) + Cvx*Cvy + (kvx-osJx-1)*Cvy + kvy-osJy,
                         Cux*b2.mra.cardI(J_y) + (kux-osJx-1)*Cuy + kuy-osIy) = val;
                   }  

                 }
               }
             }

             /* u = Wavelet x Wavelet */ 
             for(int jux = b1.j0; jux <= J_x -1; ++jux){
               Rux = b1.rangeJ(jux);
               Cux = b1.mra.cardI(jux);
               for(int juy = b2.j0; juy <= J_y -1; ++juy){
                 Ruy = b2.rangeJ(juy);
                 Cuy = b2.mra.cardI(juy);
                 for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
                   for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){

                     T val = a(wavelet, wavelet, jux, kux, juy, kuy, wavelet, wavelet, jvx, kvx, jvy, kvy);
                     if(fabs(val) > tol){
                         A(Cvx*b2.mra.cardI(J_y) + Cvx*Cvy + (kvx-osJx-1)*Cvy + kvy-osJy,
                           Cux*b2.mra.cardI(J_y) + Cux*Cuy + (kux-osJx-1)*Cuy + kuy-osJy) = val;
                     }

                   }
                 } 
               }
             }           
           }
         } 
       }
     }
    
    return A;
 
}  


template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral, typename Preconditioner>
flens::DenseVector<flens::Array<T> >
BoxProblem<T, Basis, BilinearForm, RHSIntegral, Preconditioner>::getRHS(int J_x, int J_y)
{
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    flens::DenseVector<flens::Array<T> >  N(2);
    N = b1.mra.cardI(J_x), b2.mra.cardI(J_y);
    
    int osJx = b1.rangeJ(b1.j0).firstIndex() - 1;
    int osJy = b2.rangeJ(b2.j0).firstIndex() - 1;
    int osIx = b1.mra.rangeI(b1.j0).firstIndex() - 1;
    int osIy = b2.mra.rangeI(b2.j0).firstIndex() - 1;
    
    bool spline = true;
    bool wavelet = false;
     
    flens::DenseVector<flens::Array<T> > f(b1.mra.cardI(J_x)* b2.mra.cardI(J_y));

    /*  ============  v = Scaling Fct x Scaling Fct ==========================*/
    // cout << "SF x SF : " << endl;
    Range<int> Rvx = b1.mra.rangeI(b1.j0);
    Range<int> Rvy = b2.mra.rangeI(b2.j0);
    int Cvx = b1.mra.cardI(b1.j0);
    int Cvy = b2.mra.cardI(b2.j0);
    for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
      for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
          
          f((kvx-osIx-1)*Cvy + kvy-osIy) = rhs(true, b1.j0, kvx, true, b2.j0, kvy);

      }
    }

    /* ============  v = Scaling Fct x Wavelet ==========================*/
    //cout << "SF x W : " << endl;
    Rvx = b1.mra.rangeI(b1.j0);
    Cvx = b1.mra.cardI(b1.j0);
    for(int jvy = b2.j0; jvy <= J_y - 1; ++jvy){
      Rvy = b2.rangeJ(jvy);
      Cvy = b2.mra.cardI(jvy);
      for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
        for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
          
          f(Cvx*Cvy + (kvx-osIx-1)*Cvy + kvy-osJy) = rhs(true, b1.j0, kvx, false, jvy, kvy);
           
        }
      }  
    }

    /* ============  v = Wavelet x Scaling Fct ==========================*/
    // cout << "W x SF : " << endl;
    Rvy = b2.mra.rangeI(b2.j0);
    Cvy = b2.mra.cardI(b2.j0);
    for(int jvx = b1.j0; jvx <= J_x - 1; ++jvx){
      Rvx = b1.rangeJ(jvx);
      Cvx = b1.mra.cardI(jvx);
      for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
        for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
          
          f(Cvx*b2.mra.cardI(J_y) + (kvx-osJx-1)*Cvy + kvy-osIy) 
          = rhs(false, jvx, kvx, true, b2.j0, kvy);
          
        }
      }  
    }

    /*  ============  v = Wavelet x Wavelet ==========================*/
    //cout << "W x W : " << endl;
    for(int jvx = b1.j0; jvx <= J_x -1; ++jvx){
      Rvx = b1.rangeJ(jvx);
      Cvx = b1.mra.cardI(jvx);
      for(int jvy = b2.j0; jvy <= J_y - 1; ++jvy){
        Rvy = b2.rangeJ(jvy);
        Cvy = b2.mra.cardI(jvy);
        for(int kvx = Rvx.firstIndex(); kvx <= Rvx.lastIndex(); ++kvx){
          for(int kvy = Rvy.firstIndex(); kvy <= Rvy.lastIndex(); ++kvy){
             
            f(Cvx*b2.mra.cardI(J_y) + Cvx*Cvy + (kvx-osJx-1)*Cvy + kvy-osJy) 
            = rhs(false, jvx, kvx, false, jvy, kvy);

          }
        }
      }
    }
    
    return f;
}


template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral, typename Preconditioner>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    
BoxProblem<T, Basis, BilinearForm, RHSIntegral, Preconditioner>::getPreconditioner(int J_x, int J_y)
{
    typedef typename Basis::FirstBasisType FirstBasis;
    typedef typename Basis::SecondBasisType SecondBasis;
    FirstBasis b1 = basis.first;
    SecondBasis b2 = basis.second;
    
    flens::DenseVector<flens::Array<T> >  N(2);
    N = b1.mra.cardI(J_x), b2.mra.cardI(J_y);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > D(N(1)*N(2), N(1)*N(2));
    
    int osJx = b1.rangeJ(b1.j0).firstIndex() - 1;
    int osJy = b2.rangeJ(b2.j0).firstIndex() - 1;
    int osIx = b1.mra.rangeI(b1.j0).firstIndex() - 1;
    int osIy = b2.mra.rangeI(b2.j0).firstIndex() - 1;
    
    bool spline = true;
    bool wavelet = false;

    /* SF x SF */
    Range<int> Rx = b1.mra.rangeI(b1.j0);
    Range<int> Ry = b2.mra.rangeI(b2.j0);
    int Cx = b1.mra.cardI(b1.j0);
    int Cy = b2.mra.cardI(b2.j0);    
    for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
        for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
            
            int index = (kx-osIx-1)*Cy + ky-osIy;
            D(index, index) = P(spline, b1.j0, kx, spline, b2.j0, ky);
            
        }
    }
    
    /* SF x W */
    for(int jy = b2.j0; jy <= J_y-1; ++jy){
        Ry = b2.rangeJ(jy);
        Cy = b2.mra.cardI(jy);
        for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
            for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                
                int index = Cx*Cy + (kx-osIx-1)*Cy + ky-osJy;
                D(index, index) = P(spline, b1.j0, kx, wavelet, jy, ky);
                
            }
        }
    }
    
    /* W x SF */
    Ry = b2.mra.rangeI(b2.j0);
    Cy = b2.mra.cardI(b2.j0);
    for(int jx = b1.j0; jx <= J_x -1; ++jx){
        Rx = b1.rangeJ(jx);
        Cx = b1.mra.cardI(jx);
        for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
            for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                
                int index = Cx*b2.mra.cardI(J_y) + (kx-osJx-1)*Cy + ky-osIy;
                D(index, index) = P(wavelet, jx, kx, spline, b2.j0, ky);
            
            }
        }      
    }
    
    /* W x W */
    for(int jx = b1.j0; jx <= J_x - 1; ++jx){
        Rx = b1.rangeJ(jx);
        Cx = b1.mra.cardI(jx);
        for(int jy = b2.j0; jy <= J_y - 1; ++jy){
            Ry = b2.rangeJ(jy);
            Cy = b2.mra.cardI(jy);
            for(int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx){
                for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                    
                    int index = Cx*b2.mra.cardI(J_y) + Cx*Cy + (kx-osJx-1)*Cy + ky-osJy;
                    D(index, index) = P(wavelet, jx, kx, wavelet, jy, ky);
                
                }
            } 
        }
    }
    
    return D;
    
}
    
} // namespace lawa