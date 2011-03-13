/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
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

#include "soe/full_linear_soe.h"

FullLinearSOE::FullLinearSOE()
  :SOE() {
  myTag = TAG_SOE_LINEAR_FULL;
}
FullLinearSOE::~FullLinearSOE() {
}
int FullLinearSOE::
insertMatrixIntoA(const Matrix& Ke, const IDContainer& EFTable, double factor) {
  isLUFactored = false;
  for (unsigned i = 0; i < EFTable.size(); i++)
    for (unsigned j = 0; j < EFTable.size(); j++) {
      if ((EFTable[i]) < 0) continue;
      if ((EFTable[j]) < 0) continue;
      A[EFTable[j]*size_+EFTable[i]]+=factor*Ke(i, j);
    }
  return 0;
}
void FullLinearSOE::set_size() {
  // If the size has not changed do not resize arrays
  if (size_ == pA->get_model()->get_num_eqns()) {
    return;
  } else {
    size_ = pA->get_model()->get_num_eqns();
  }
  A.resize(size_*size_);
  B.resize(size_);
  X.resize(size_);
  IPIV.resize(size_);

  size_t d = (size_*size_+2*size_)*sizeof(double);
  size_t i = size_*sizeof(int);
  printf("soe: Allocated %6.2fmb of memory for %d dofs.\n",
    static_cast<double>(d+i)/(1024*1024), size_);
}
void FullLinearSOE::print() {
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++)  cout << A[j*size_+i] << ' ';
        cout << " | " << B[i] << endl;
  }
}
int FullLinearSOE::get_eigen_sign() {
  for (int i = 0; i < size_; i++)
    if (fabs(A[i*(size_+1)]) < 1e-12)  return 0;
    else if (A[i*(size_+1)]  < 0)      return -1;
  return 1;
}
int FullLinearSOE::solve() {
  this->print();
  int N = size_;
  int NRHS = 1;
  int LDA = N;
  int LDB = N;
  int INFO;
  char c='N';
//  this->print();
//  for (int i = 0;i < size_;i++)
//  {
//        cout << B[i]<<endl;
//  }
//  cout << endl;
  X = B;
  // Compute the LU factorization of the band matrix A.
  if (!isLUFactored) {
    dgetrf(&N, &N, &A[0], &LDA, &IPIV[0], &INFO);
    if (INFO != 0)
      throw SException("[nemesis:%d] %s", 1101, "SOE: lapack DGETRF failed.\n");
    isLUFactored = true;
  }
  // Solve the system A*X = B, overwriting B with X.
  dgetrs(&c, &N, &NRHS, &A[0], &LDA, &IPIV[0], &X[0], &LDB, &INFO, 1);
  // for (int i = 0;i < size_;i++) cout << B[i]<<endl; cout << endl;
    return 0;
}
