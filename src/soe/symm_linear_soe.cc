/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2011 F.E.Karaoulanis [http://www.nemesis-project.org]     *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "soe/symm_linear_soe.h"
#include <stdio.h>

SymmLinearSOE::SymmLinearSOE()
  :SOE() {
  myTag = TAG_SOE_LINEAR_SYMM;
}
SymmLinearSOE::~SymmLinearSOE() {
}
int SymmLinearSOE::
insertMatrixIntoA(const Matrix& Ke, const IDContainer& EFTable, double factor) {
  isLUFactored = false;
  for (unsigned i = 0; i < EFTable.size(); i++)
    for (unsigned j = 0; j <=i ; j++) {
      if (EFTable[i] < 0) continue;
      if (EFTable[j] < 0) continue;
      int ii = EFTable[i];
      int jj = EFTable[j];
      if (EFTable[i] < EFTable[j]) {
        ii = EFTable[j];
        jj = EFTable[i];
      }
      A[jj*size_-static_cast<int>((jj+1)*jj*0.5)+ii]+=factor*Ke(i, j);
    }
  return 0;
}
void SymmLinearSOE::set_size() {
  // If the size has not changed do not resize arrays
  if (size_ == pA->get_model()->get_num_eqns()) {
    return;
  } else {
    size_ = pA->get_model()->get_num_eqns();
  }
  A.resize(static_cast<int>(0.5*size_*(size_+1)));
  B.resize(size_);
  X.resize(size_);
  IPIV.resize(size_);

  size_t d = (static_cast<int>(0.5*size_*(size_+1))+size_+size_)*sizeof(double);
  size_t i = size_*sizeof(int);
  printf("soe: Allocated %6.2fmb of memory for %d dofs.\n",
    static_cast<double>(d+i)/(1024*1024), size_);
}
void SymmLinearSOE::print() {
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      if (j <= i) {
        printf("%12.4f ", A[j*size_-static_cast<int>((j+1)*j*0.5)+i]);
      } else {
        printf("%12.4f ", A[i*size_-static_cast<int>((i+1)*i*0.5)+j]);
      }
      printf("| %12.4f\n", B[i]);
    }
  }
}
int SymmLinearSOE::get_eigen_sign() {
/// @todo Unfinished function.
//  for (int i = 0;i < size_;i++)
//    if (fabs(A(i, i))<1e-12)  return 0;
//    else if (A(i, i)<0)   return -1;
  return 1;
}
int SymmLinearSOE::solve() {
  int N = size_;
  int NRHS = 1;
  int INFO;
  int LDB = N;
  char c = 'L';
  X = B;
  if (!isLUFactored) {
    dsptrf(&c, &N, &A[0], &IPIV[0], &INFO, 1);
    if (INFO != 0)
      throw SException("[nemesis:%d] %s", 1102, "SOE: lapack DSPTRF failed.");
    isLUFactored = true;
  }
  // Solve the system A*X = B, overwriting B with X.
  dsptrs(&c, &N, &NRHS, &A[0], &IPIV[0], &X[0], &LDB, &INFO, 1);
  return 0;
}
