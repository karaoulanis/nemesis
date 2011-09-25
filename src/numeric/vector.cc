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

#include "numeric/vector.h"
#include "numeric/lapack.h"
#include "numeric/numeric.h"

Vector::Vector()
:size_(0), data_(0) {
}

Vector::Vector(int size)
:size_(size), data_(0) {
  try {
    data_ = new double[size_];
  } catch(std::bad_alloc) {
    throw SException("[nemesis:%d] %s", 1001, "Run out of memory.\n");
  }
}

Vector::Vector(int size, double value)
:size_(size), data_(0) {
  try {
    data_ = new double[size_];
  } catch(std::bad_alloc) {
    throw SException("[nemesis:%d] %s", 1001, "Run out of memory.\n");
  }
  for (int i = 0; i < size_; i++) {
    data_[i] = value;
  }
}

Vector::Vector(const Vector& vector)
:size_(vector.size_), data_(0) {
  if (size_ != 0) {
    try {
      data_ = new double[size_];
    } catch(std::bad_alloc) {
      throw SException("[nemesis:%d] %s", 1001, "Run out of memory.\n");
    }
    for (int i = 0; i < size_; i++) {
      data_[i] = vector.data_[i];
    }
  }
}

Vector::~Vector() {
  if (data_ != 0) delete[] data_;
}

Vector& Vector::operator=(const Vector& v) {
  #ifdef _DEBUG
  num::check::array_size(v.size_, size_);
  #endif
  if (this != &v)
    for (int i = 0;i < size_;i++) data_[i]=v.data_[i];
  return *this;
}

const Vector& Vector::Eigenvalues() {
  #ifdef _DEBUG
  num::check::array_size(size_, 6);
  #endif
  int N = 3;
  char JOBZ='N';
  char UPLO='L';
  int LDA = 3;
  int LWORK = 102;
  int INFO;
  static Vector res(3);
static Vector A(9);
  static Vector WORK(LWORK);
  res.Clear();

  A[0]=data_[0];
  A[1]=data_[3];
  A[2]=data_[5];

  A[3]=0.;
  A[4]=data_[1];
  A[5]=data_[4];

  A[6]=0.;
  A[7]=0.;
  A[8]=data_[2];

  dsyev(&JOBZ, &UPLO, &N, A.get_data(), &LDA, res.get_data(), WORK.get_data(),
        &LWORK, &INFO, 1, 1);
  //  dsyev(&JOBZ, &UPLO, &N, A.data(), &LDA, res.data(), WORK.data(), &LWORK,
  //        &INFO);
  //  cout << "Optimal LWORK : "<<WORK[0]<<endl;
  double d = res[0];
  res[0]=res[2];
  res[2]=d;
  return res;
}
