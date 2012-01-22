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

#include "soe/band_linear_soe.h"
#include <stdio.h>
#include "model/model.h"
#include "numeric/lapack.h"
#include "numeric/matrix.h"

BandLinearSOE::BandLinearSOE()
    : SOE(),
      lowerBandwidth(0),
      upperBandwidth(0),
      nRows(0) {
  myTag = TAG_SOE_LINEAR_BAND;
}

BandLinearSOE::BandLinearSOE(Model* model)
    : SOE(model),
      lowerBandwidth(0),
      upperBandwidth(0),
      nRows(0) {
  myTag = TAG_SOE_LINEAR_BAND;
}

BandLinearSOE::~BandLinearSOE() {
}

int BandLinearSOE::
InsertMatrixIntoA(const Matrix& Ke, const IDContainer& EFTable, double factor) {
  isLUFactored = false;
  for (unsigned i = 0;i < EFTable.size();i++)
    for (unsigned j = 0; j < EFTable.size(); j++) {
      if (EFTable[i] < 0) continue;
      if (EFTable[j] < 0) continue;
      unsigned index = EFTable[j]*nRows + lowerBandwidth + upperBandwidth +
                       EFTable[i]-EFTable[j];
      A[index]+=factor*Ke(i, j);
    }
  return 0;
}

void BandLinearSOE::set_size() {
  UndirectedGraph G;
  model_->get_undirected_graph(G);
  lowerBandwidth = bandwidth(G);
  upperBandwidth = lowerBandwidth;

  // If the size has not changed do not resize arrays
  if (size_ == model_->get_num_eqns()) {
    return;
  }

  size_ = model_->get_num_eqns();
  nRows = 2*lowerBandwidth+upperBandwidth+1;
  A.resize(nRows*size_);
  B.Resize(size_);
  X.Resize(size_);
  IPIV.resize(size_);

  size_t d = (nRows+2)*size_*sizeof(double);
  size_t i = size_*sizeof(int);
  printf("soe: Allocated %6.2fmb of memory for %d dofs.\n",
    static_cast<double>(d+i)/(1024*1024), size_);
}

void BandLinearSOE::print() {
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      printf("%12.4f ", A[j*nRows+(lowerBandwidth+upperBandwidth+i-j)]);
    }
    printf("| %12.4f\n", B[i]);
  }
}

int BandLinearSOE::get_eigen_sign() {
  for (int i = 0;i < size_;i++)
    if (fabs(A[lowerBandwidth+upperBandwidth+i*size_]) < 1e-12)  return  0;
    else if (A[lowerBandwidth+upperBandwidth+i*size_]  < 0)      return -1;
  return 1;
}

int BandLinearSOE::solve() {
  // Input data
  int N = size_;
  int KL = lowerBandwidth;
  int KU = upperBandwidth;
  int NRHS = 1;
  int INFO;
  int LDAB = 2*KL+KU+1;
  int LDB = N;
  char c='N';
  X = B;

  // Compute the LU factorization of the band matrix A.
  if (!isLUFactored) {
    dgbtrf(&N, &N, &KL, &KU, &A[0], &LDAB, &IPIV[0], &INFO);
    if (INFO != 0)
      throw SException("[nemesis:%d] %s", 1103, "SOE: lapack DGBTRF failed.\n");
    isLUFactored = true;
  }
  // Solve the system A*X = B, overwriting B with X.
  dgbtrs(&c, &N, &KL, &KU, &NRHS, &A[0], &LDAB, &IPIV[0], &X[0], &LDB, &INFO, 1);
  return 0;
}
