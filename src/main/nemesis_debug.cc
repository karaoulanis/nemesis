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

#include "main/nemesis_debug.h"
#include <stdio.h>
#include "numeric/lapack.h"

void report(const double d, const char* name, int total, int decimal) {
  printf("%s =", name);
  num::print_d(d, total, decimal);
  printf("\n");
}

void report(const Matrix& m, const char* name, int total, int decimal) {
  printf("%s [%d, %d] =\n", name, m.get_rows(), m.get_cols());
  for (int i = 0; i < m.get_rows(); i++) {
    for (int j = 0; j < m.get_cols(); j++) {
      num::print_d(m(i, j), total, decimal);
    }
    printf("\n");
  }
}

void report(const Vector& v, const char* name, int total, int decimal) {
  printf("%s [%d] =", name, v.get_size());
  for (int i = 0; i < v.get_size(); i++) {
    num::print_d(v[i], total, decimal);
  }
  printf("\n");
}

void add(Matrix* K, int row, int col, const Matrix& B1, const Matrix& C,
         const Matrix B2, double c1, double c0) {
  int m = B1.get_rows();
  int n = B2.get_cols();
  int pos     = row*K->get_cols()+col;
  int colK    = K->get_cols();
  double* pK  = K->get_data();
  double* pB1 = B1.get_data();
  double* pB2 = B2.get_data();
  double* pC = C.get_data();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      pK[pos+i*colK+j]*=c0;
  for (int i = 0; i < n; i++)
    for (int k = 0; k < m; k++)
      for (int l = 0; l < m; l++)
        for (int j = 0; j < n; j++)
          pK[pos+i*colK+j]+=c1*pB1[k*n+i]*pC[k*m+l]*pB2[l*n+j];
}

void add(Vector* R, int row, const Matrix& BT, const Vector& V,
         double c1, double c0) {
  int m = BT.get_rows();
  int n = BT.get_cols();
  double* pV = V.get_data();
  double* pR = R->get_data();
  double* pBT = BT.get_data();
  for (int i = 0; i < n; i++)
    pR[i]*=c0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      pR[row+i]+=c1*pBT[j*n+i]*pV[j];
}

void add2(Vector* R, int row, const Matrix& BT, const Vector& V,
          double c1, double c0) {
  int m = BT.get_rows();
  int n = BT.get_cols();
  double* pV = V.get_data();
  double* pR = R->get_data();
  double* pBT = BT.get_data();
  for (int i = 0; i < m; i++)
    pR[i]*=c0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      pR[i]+=c1*pBT[i*n+j]*pV[row+j];
}

// ****************************************************************************
//
// ****************************************************************************
void add_BTCB(Matrix* K, int row, int col, const int* perm,
              const Matrix& B1, const Matrix& C, const Matrix B2,
              double c1, double c0) {
  int m = B1.get_rows();
  int n = B2.get_cols();
  int pos = row*K->get_cols()+col;
  int colK = K->get_cols();
  int colC = C.get_cols();
  double* pK = K->get_data();
  double* pB1 = B1.get_data();
  double* pB2 = B2.get_data();
  double* pC = C.get_data();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      pK[pos+i*colK+j]*=c0;
  for (int i = 0; i < n; i++)
    for (int k = 0; k < m; k++)
      for (int l = 0; l < m; l++)
        for (int j = 0; j < n; j++)
          pK[pos+i*colK+j]+=c1*pB1[k*n+i]*pC[perm[k]*colC+perm[l]]*pB2[l*n+j];
}

void add_BTv(Vector* R, int row, const int* perm,
             const Matrix& B, const Vector& v, double c1, double c0) {
  int m = B.get_rows();
  int n = B.get_cols();
  double* pV = v.get_data();
  double* pR = R->get_data();
  double* pB = B.get_data();
  for (int i = 0; i < n; i++)
    pR[i]*=c0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      pR[row+i]+=c1*pB[j*n+i]*pV[perm[j]];
}

void add_Bv(Vector* R, int row, const int* perm,
            const Matrix& B, const Vector& v, double c1, double c0) {
  int m = B.get_rows();
  int n = B.get_cols();
  double* pV = v.get_data();
  double* pR = R->get_data();
  double* pB = B.get_data();
  for (int i = 0; i < m; i++)
    pR[i]*=c0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      pR[perm[i]]+=c1*pB[i*n+j]*pV[row+j];
}

void spectralDecomposition(const Vector& s, Vector& sP, Matrix& sV) {
  int N = 3;
  char JOBZ='V';
  char UPLO='L';
  int LDA = 3;
  int LWORK = 102;
  int INFO;
  static Vector WORK(LWORK);

  sV(0, 0) = s[0];
  sV(0, 1) = s[3];
  sV(0, 2) = s[5];
  sV(1, 0) = 0.;
  sV(1, 1) = s[1];
  sV(1, 2) = s[4];
  sV(2, 0) = 0.;
  sV(2, 1) = 0.;
  sV(2, 2) = s[2];

  dsyev(&JOBZ, &UPLO, &N, sV.get_data(), &LDA, sP.get_data(), WORK.get_data(),
        &LWORK, &INFO);
  // cout << WORK[0]<<endl;
  double d, d0, d1, d2;
  d  = sP[0];
  d0 = sV(0, 0);
  d1 = sV(0, 1);
  d2 = sV(0, 2);
  sP[0]    = sP[2];
  sV(0, 0) = sV(2, 0);
  sV(0, 1) = sV(2, 1);
  sV(0, 2) = sV(2, 2);
  sP[2]    = d;
  sV(2, 0) = d0;
  sV(2, 1) = d1;
  sV(2, 2) = d2;
}
