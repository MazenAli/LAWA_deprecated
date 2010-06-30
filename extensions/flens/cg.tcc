/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#include <cassert>
#include <limits>

namespace lawa {

template <typename MA, typename VX, typename VB>
int
cg(const MA &A, VX &x, const VB &b, typename _cg<VB>::T tol,
   long maxIterations)
{
    typename _cg<VB>::T alpha, beta, rNormSquare, rNormSquarePrev;
    typename _cg<VB>::AuxVector Ap, r, p;

    if (x.length()!=A.numCols()) {
        x.engine().resize(A.numCols());
    }

    r = A*x - b;
    p = -1*r;
    rNormSquare = r*r;
    for (long k=1; k<=maxIterations; k++) {
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << sqrt(rNormSquare)
                << std::endl;
        #endif
        if (sqrt(rNormSquare)<=tol) {
            return k-1;
        }
        Ap = A*p;
        alpha = rNormSquare/(p * Ap);
        x += alpha*p;
        r += alpha*Ap;

        rNormSquarePrev = rNormSquare;
        rNormSquare = r*r;
        beta = rNormSquare/rNormSquarePrev;
        p = beta*p - r;
    }
    return maxIterations;
}

// Algorithm 9.2, Y. Saad: Iterative Methods for Sparse Linear Systems
// for solving Ax=b with P^T A P u = P^T b, u=P^{-1} x
// Note the role of P and P^T is switched.
// This algorithm uses split preconditioning defined as in
// (6.2), K. Urban: Wavelet Methods for Elliptic PDEs
// Also non-symmetric preconditioning is possible as mentioned in remark 6.2
template <typename Prec, typename MA, typename VX, typename VB>
int
pcg(const Prec &P, const MA &A, VX &x, const VB &b,
    typename _cg<VB>::T tol, long maxIterations)
{
    typename _cg<VB>::T pNormSquare, alpha, beta, rHatq, rHatqPrev;
    typename _cg<VB>::AuxVector r, rHat, p, Ap;

    if (x.length()!=A.numCols()) {
        x.engine().resize(A.numCols());
    }
    r = b - A*x;
    rHat = transpose(P)*r;
    p = P*rHat;
    rHatq = rHat*rHat;
    for (long k=1; k<=maxIterations; k++) {
        Ap = A*p;
        alpha = rHatq/(Ap*p);
        x += alpha*p;
        rHat = rHat - alpha * transpose(P)*Ap;
        rHatqPrev = rHatq;
        rHatq =rHat*rHat;
        beta = rHatq/rHatqPrev;
        p *= beta;
        p += P*rHat;
        pNormSquare = p*p;
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << sqrt(pNormSquare)
                << std::endl;
        #endif
        if (sqrt(pNormSquare)<=tol) {
            return k-1;
        }
    }
    return maxIterations;
}

/*
template <typename Prec, typename MA, typename VX, typename VB>
int
pcg(const Prec &B, const MA &A, VX &x, const VB &b,
    typename _cg<VB>::T tol, long maxIterations)
{
    typename _cg<VB>::T pNormSquare, alpha, beta, rq, rqPrev;
    typename _cg<VB>::AuxVector r, q, p, Ap;

    r = A*x - b;
    q = B*r;

    p = q;
    // TODO: next line results in an error with T = long double. WHY???
    // p = -q;
    p *= -1;
    rq = r*q;

    for (long k=1; k<=maxIterations; k++) {
        pNormSquare = p*p;
        if (sqrt(pNormSquare)<tol) {
            return k-1;
        }
        Ap = A*p;
        alpha = rq/(p*Ap);
        x += alpha*p;

        r += alpha*Ap;
        q = B*r;

        rqPrev = rq;
        rq = r*q;
        beta = rq/rqPrev;
        p = beta*p - q;
    }
    return maxIterations;
}
*/
} // namespace lawa
