/*
 * eqnsys.h - equations system solver class definitions
 *
 * Copyright (C) 2004, 2005, 2006 Stefan Jahn <stefan@lkcc.org>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __EQNSYS_H__
#define __EQNSYS_H__

enum algo_type {
  ALGO_INVERSE = 0x0001,
  ALGO_GAUSS = 0x0002,
  ALGO_GAUSS_JORDAN = 0x0004,
  ALGO_LU_FACTORIZATION_CROUT = 0x0008,
  ALGO_LU_FACTORIZATION_DOOLITTLE = 0x0010,
  ALGO_LU_SUBSTITUTION_CROUT = 0x0020,
  ALGO_LU_SUBSTITUTION_DOOLITTLE = 0x0040,
  ALGO_LU_DECOMPOSITION = 0x0028,
  ALGO_LU_DECOMPOSITION_CROUT = 0x0028,
  ALGO_LU_DECOMPOSITION_DOOLITTLE = 0x0050,
  ALGO_JACOBI = 0x0080,
  ALGO_GAUSS_SEIDEL = 0x0100,
  ALGO_SOR = 0x0200,
  ALGO_QR_DECOMPOSITION = 0x0400,
  ALGO_QR_DECOMPOSITION_LS = 0x0800,
  ALGO_SV_DECOMPOSITION = 0x1000,
  // testing
  ALGO_QR_DECOMPOSITION_2 = 0x2000,
};

enum pivot_type {
  PIVOT_NONE = 0x01,
  PIVOT_PARTIAL = 0x02,
  PIVOT_FULL = 0x04,
};

#include "tmatrix.h"
#include "tvector.h"

namespace qucs {

template <class nr_type_t> class eqnsys {
public:
  eqnsys();
  eqnsys(eqnsys &) = delete;
  ~eqnsys();
  void setAlgo(int a) { algo = a; }
  int getAlgo() { return algo; }
  void passEquationSys(tmatrix<nr_type_t> *, tvector<nr_type_t> *, tvector<nr_type_t> *);
  void solve();

private:
  int update;
  int algo;
  int pivoting;
  int *rMap;
  int *cMap;
  int N;
  double *nPvt;

  tmatrix<nr_type_t> *A;
  tmatrix<nr_type_t> *V;
  tvector<nr_type_t> *B;
  tvector<nr_type_t> *X;
  tvector<nr_type_t> *R;
  tvector<nr_type_t> *T;
  tvector<double> *S;
  tvector<double> *E;

  void solve_inverse();
  void solve_gauss();
  void solve_gauss_jordan();
  void solve_lu_crout();
  void solve_lu_doolittle();
  void factorize_lu_crout();
  void factorize_lu_doolittle();
  void substitute_lu_crout();
  void substitute_lu_doolittle();
  void solve_qr();
  void solve_qr_ls();
  void solve_qrh();
  void factorize_qrh();
  void substitute_qrh();
  void factorize_qr_householder();
  void substitute_qr_householder();
  void substitute_qr_householder_ls();
  nr_type_t householder_create_left(int);
  void householder_apply_left(int, nr_type_t);
  nr_type_t householder_left(int);
  nr_type_t householder_create_right(int);
  void householder_apply_right(int, nr_type_t);
  void householder_apply_right_extern(int, nr_type_t);
  nr_type_t householder_right(int);
  double euclidian_c(int, int r = 1);
  double euclidian_r(int, int c = 1);
  void givens_apply_u(int, int, double, double);
  void givens_apply_v(int, int, double, double);
  void solve_svd();
  void chop_svd();
  void factorize_svd();
  void substitute_svd();
  void diagonalize_svd();
  void solve_iterative();
  void solve_sor();
  double convergence_criteria();
  void ensure_diagonal();
  void ensure_diagonal_MNA();
  int countPairs(int, int &, int &);
  void preconditioner();
};

} // namespace qucs

#include "eqnsys.cpp"

#endif /* __EQNSYS_H__ */
