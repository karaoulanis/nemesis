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

#include "soe/eigen_soe.h"
#include <stdio.h>

EigenSOE::EigenSOE()
  :SOE() {
  myTag = TAG_NONE;
}

EigenSOE::~EigenSOE() {
}

int EigenSOE::
insertMatrixIntoA(const Matrix& Ke, const IDContainer& EFTable, double factor) {
  isLUFactored = false;
  for (unsigned i = 0;i < EFTable.size();i++)
    for (unsigned j = 0; j < EFTable.size(); j++) {
      if (EFTable[i] < 0) continue;
      if (EFTable[j] < 0) continue;
      A[EFTable[j]*size_+EFTable[i]]+=factor*Ke(i, j);
    }
  return 0;
}

int EigenSOE::
insertMatrixIntoM(const Matrix& Ke, const IDContainer& EFTable, double factor) {
  isLUFactored = false;
  for (unsigned i = 0;i < EFTable.size();i++)
    for (unsigned j = 0; j < EFTable.size(); j++) {
      if (EFTable[i] < 0) continue;
      if (EFTable[j] < 0) continue;
      M[EFTable[j]*size_+EFTable[i]]+=factor*Ke(i, j);
    }
  return 0;
}

void EigenSOE::set_size() {
  // If the size has not changed do not resize arrays
  if (size_ == pA->get_model()->get_num_eqns()) {
    return;
  } else {
    size_ = pA->get_model()->get_num_eqns();
  }
  A.resize(size_*size_);
  M.resize(size_*size_);
  X.resize(size_);
  ALPHAR.resize(size_);
  ALPHAI.resize(size_);
  BETA.resize(size_);
  VL.resize(size_*size_);
  VR.resize(size_*size_);
  WORK.resize(8*size_);

  size_t d = (4*size_*size_+9*size_)*sizeof(double);
  size_t i = size_*sizeof(int);
  printf("soe: Allocated %6.2fmb of memory.\n",
    static_cast<double>(d+i)/(1024*1024));
}

void EigenSOE::print() {
  printf("A\n");
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      printf("%12.4f ", A[j*size_+i]);
    }
    printf("\n");
  }
  printf("M\n");
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      printf("%12.4f ", M[j*size_+i]);
    }
    printf("\n");
  }
}

int EigenSOE::get_eigen_sign() {
  return 1;
}

void EigenSOE::zeroM() {
  for (unsigned i = 0;i < M.size(); i++) M[i]=0;
}

int EigenSOE::solve() {
  char JOBVL='N';
  char JOBVR='V';
  int N = size_;
  int LDA = N;
  int LDB = N;
  int LDVL = N;
  int LDVR = N;
  int LWORK = 8*N;
  int INFO;
  dggev(&JOBVL, &JOBVR, &N, &A[0], &LDA, &M[0], &LDB, &ALPHAR[0], &ALPHAI[0],
    &BETA[0], &VL[0], &LDVL, &VR[0], &LDVR, &WORK[0], &LWORK, &INFO, 1, 1);
  if (INFO != 0)
    throw SException("[nemesis:%d] %s", 1110, "SOE: lapack DSYGV failed.");
  for (int i = 0;i < size_;i++) X[i]=ALPHAR[i]/BETA[i];
    return 0;
}
