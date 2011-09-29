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

// Included files
#include <gtest/gtest.h>
#include "numeric/matrix.h"

class MatrixTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

// test: Matrix::Matrix();
TEST_F(MatrixTest, DefaultConstructor) {
  Matrix m;
  EXPECT_EQ(m.get_rows(),    0) << "Matrix rows are not 0 initialized";
  EXPECT_EQ(m.get_cols(),    0) << "Matrix cols are not 0 initialized";
  EXPECT_EQ(m.get_data(),
    static_cast<double*>(NULL)) << "Matrix data is  not 0 initialized";
  EXPECT_EQ(m.get_size(),    0) << "Matrix size is  not 0 initialized";
}

// test: Matrix::Matrix(int rows, int cols);
TEST_F(MatrixTest, ConstructorNoValues) {
  Matrix m(3, 2);
  EXPECT_EQ(m.get_rows(), 3) << "Matrix rows are not set correctly";
  EXPECT_EQ(m.get_cols(), 2) << "Matrix cols are not set correctly";
  EXPECT_EQ(m.get_size(), 6) << "Matrix size is  not set correctly";
}

// test: Matrix::Matrix(int rows, int cols, double value);
TEST_F(MatrixTest, ConstructorWithValues) {
  Matrix m(3, 2, 2.);
  EXPECT_EQ(m.get_rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m.get_cols(), 2) << "Matrix cols not set correctly";
  EXPECT_EQ(m.get_size(), 6) << "Matrix size is  not set correctly";
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m(i, j), 2.) << "Matrix value error";
    }
  }
}

// test: Matrix::Matrix(const Matrix& matrix);
TEST_F(MatrixTest, CopyConstructor) {
  Matrix m1(3, 2, 2.);
  Matrix m2(m1);
  EXPECT_EQ(m2.get_rows(), 3) << "Matrix rows are not set correctly";
  EXPECT_EQ(m2.get_cols(), 2) << "Matrix cols are not set correctly";
  EXPECT_EQ(m2.get_size(), 6) << "Matrix size is  not set correctly";
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m2(i, j), 2.) << "Matrix value error";
    }
  }
}

// test: ~Matrix();

// test: Matrix::Matrix& operator=(const Matrix& m);
TEST_F(MatrixTest, CopyAssignment) {
  Matrix m1(3, 2, 2.);
  Matrix m2 = m1;
  EXPECT_EQ(m2.get_rows(), 3) << "Matrix rows are not set correctly";
  EXPECT_EQ(m2.get_cols(), 2) << "Matrix cols are not set correctly";
  EXPECT_EQ(m2.get_size(), 6) << "Matrix size is  not set correctly";
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m2(i, j), 2.) << "Matrix value error";
    }
  }
}

// inline double& operator()(int i, int j) {
TEST_F(MatrixTest, AccessOperator) {
  Matrix m(3, 2);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      m(i, j) = j*3. + i;
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m(i, j), j*3. + i) << "Matrix value error";
    }
  }
}

/// @todo:
// inline double operator()(int i, int j) const {

// inline int get_rows() const {
TEST_F(MatrixTest, GetRows) {
  Matrix m(3, 2, 2.);
  EXPECT_EQ(m.get_rows(), 3) << "Matrix rows are not set correctly";
}

// inline int get_cols() const {
TEST_F(MatrixTest, GetCols) {
  Matrix m(3, 2, 2.);
  EXPECT_EQ(m.get_cols(), 2) << "Matrix cols are not set correctly";
}

/// @todo
// inline double* get_data() {
// inline double* get_data() const {
// inline double* get_data() const {

// test: void Matirx::resize(int m, int n)
TEST_F(MatrixTest, MatrixResize) {
  Matrix m(3, 2);
  m.Resize(2, 3);
  EXPECT_EQ(m.get_rows(), 2) << "Matrix rows are not set correctly";
  EXPECT_EQ(m.get_cols(), 3) << "Matrix cols are not set correctly";
  EXPECT_EQ(m.get_size(), 6) << "Matrix size is  not set correctly";
}

// test: void Matirx::resize(int m, int n double c)
TEST_F(MatrixTest, MatrixResizeWithValues) {
  Matrix m(2, 3);
  m.Resize(3, 2, 2.);
  EXPECT_EQ(m.get_rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m.get_cols(), 2) << "Matrix cols not set correctly";
  EXPECT_EQ(m.get_size(), 6) << "Matrix size is  not set correctly";
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m(i, j), 2.) << "Matrix value error";
    }
  }
}

// test: void Matrix::Clear()
TEST_F(MatrixTest, MatrixClear) {
  Matrix m(2, 3, 2.);
  m.Clear();
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m(i, j), 0.) << "Matrix value error";
    }
  }
}

// inline Matrix& operator+=(double value) {
TEST_F(MatrixTest, MatrixOperatorUnaryPlus) {
  Matrix m(2, 3, 2.);
  m += 1.;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m(i, j), 3.) << "Matrix value error";
    }
  }
}

// inline Matrix& operator-=(double value) {
TEST_F(MatrixTest, MatrixOperatorUnaryMinus) {
  Matrix m(2, 3, 2.);
  m -= 1.;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m(i, j), 1.) << "Matrix value error";
    }
  }
}

// inline Matrix& operator*=(double value) {
TEST_F(MatrixTest, MatrixOperatorUnaryProduct) {
  Matrix m(2, 3, 2.);
  m *= 2.;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m(i, j), 4.) << "Matrix value error";
    }
  }
}

// inline Matrix& operator/=(double value) {
TEST_F(MatrixTest, MatrixOperatorUnaryDivision) {
  Matrix m(2, 3, 2.);
  m /= 2.;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m(i, j), 1.) << "Matrix value error";
    }
  }
}

// inline Matrix operator*(const double value) const {
TEST_F(MatrixTest, MatrixScalarProduct1) {
  Matrix m1(2, 3, 2.);
  Matrix m2(2, 3);
  m2 = m1*2.0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m2(i, j), 4.) << "Matrix value error";
    }
  }
}

// inline friend Matrix operator*(const double c, const Matrix& m) {
TEST_F(MatrixTest, MatrixScalarProduct2) {
  Matrix m1(2, 3, 2.);
  Matrix m2(2, 3);
  m2 = 2.0*m1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m2(i, j), 4.) << "Matrix value error";
    }
  }
}

// inline Matrix operator/(const double c) const {
TEST_F(MatrixTest, MatrixScalarDivision) {
  Matrix m1(2, 3, 2.);
  Matrix m2(2, 3);
  m2 = m1/2.0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m2(i, j), 1.) << "Matrix value error";
    }
  }
}

// inline void add_cM(double c, const Matrix& m, double c0 = 1.0) {
// inline void add_cMM(double c, const Matrix& m1, const Matrix& m2,

// inline Matrix& operator+=(const Matrix& m) {
TEST_F(MatrixTest, MatrixPlusEqualMatrix) {
  Matrix m1(2, 3, 2.);
  Matrix m2(2, 3, 1.);
  m1 += m2;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m1(i, j), 3.) << "Matrix value error";
    }
  }
}

// inline Matrix& operator-=(const Matrix& m) {
TEST_F(MatrixTest, MatrixMinusEqualMatrix) {
  Matrix m1(2, 3, 2.);
  Matrix m2(2, 3, 1.);
  m1 -= m2;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m1(i, j), 1.) << "Matrix value error";
    }
  }
}

// inline Matrix operator+(const Matrix& m) const {
TEST_F(MatrixTest, MatrixPlusMatrix) {
  Matrix m0(2, 3);
  Matrix m1(2, 3, 2.);
  Matrix m2(2, 3, 1.);
  m0 = m1 + m2;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m0(i, j), 3.) << "Matrix value error";
    }
  }
}

// inline Matrix operator-(const Matrix& m) const {
TEST_F(MatrixTest, MatrixMinusMatrix) {
  Matrix m0(2, 3);
  Matrix m1(2, 3, 2.);
  Matrix m2(2, 3, 1.);
  m0 = m1 - m2;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m0(i, j), 1.) << "Matrix value error";
    }
  }
}

// inline Matrix operator*(const Matrix& m) const {
TEST_F(MatrixTest, MatrixDotMatrix) {
  Matrix m0(2, 2);
  Matrix m1(2, 2, 2.);
  Matrix m2(2, 2, 2.);
  m0 = m1 * m2;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m0(i, j), 8.) << "Matrix value error";
    }
  }
}

// inline Matrix& operator+() {
TEST_F(MatrixTest, PlusMatrix) {
  Matrix m(2, 3, 2.);
  m =  +m;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m(i, j), 2.) << "Matrix value error";
    }
  }
}

// inline Matrix operator-() {
TEST_F(MatrixTest, MinusMatrix) {
  Matrix m(2, 3, 2.);
  m =  -m;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m(i, j), -2.) << "Matrix value error";
    }
  }
}

// inline friend Matrix Transpose(const Matrix& m) {
TEST_F(MatrixTest, Transpose) {
  Matrix m(2, 3);
  m(0, 0) = 11;
  m(0, 1) = 12;
  m(0, 2) = 13;
  m(1, 0) = 21;
  m(1, 1) = 22;
  m(1, 2) = 23;
  Matrix mT = Transpose(m);
  EXPECT_EQ(mT(0, 0), m(0, 0)) << "Matrix value error";
  EXPECT_EQ(mT(0, 1), m(1, 0)) << "Matrix value error";
  EXPECT_EQ(mT(1, 0), m(0, 1)) << "Matrix value error";
  EXPECT_EQ(mT(1, 1), m(1, 1)) << "Matrix value error";
  EXPECT_EQ(mT(2, 0), m(0, 2)) << "Matrix value error";
  EXPECT_EQ(mT(2, 1), m(1, 2)) << "Matrix value error";
}

// inline Vector operator*(const Vector& v) const {
TEST_F(MatrixTest, MatrixVector) {
  Matrix m(2, 3);
  m(0, 0) = 11.;
  m(0, 1) = 12.;
  m(0, 2) = 13.;
  m(1, 0) = 21.;
  m(1, 1) = 22.;
  m(1, 2) = 23.;
  Vector v(3);
  v[0] = 1.;
  v[1] = 2.;
  v[2] = 3.;
  Vector res = m*v;
  EXPECT_EQ(res[0],  74.) << "Vector value error";
  EXPECT_EQ(res[1], 134.) << "Vector value error";
  
}
// inline Matrix& Append(const Matrix& m, int row, int col, double c = 1.0,
//                         double c0 = 0.) {
TEST_F(MatrixTest, Append) {
  Matrix m1;
  Matrix m2;
  // append at 1,1 with c = 1.0 (default) and c0 = 0. (default).
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  2  2]
  //     [ 1  1  1]     [ 2  2]         [ 1  2  2]
  m1.Resize(3, 3, 1.);
  m2.Resize(2, 2., 2.);
  m1.Append(m2, 1, 1);
  EXPECT_EQ(m1(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m1(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m1(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(1, 1), 2.) << "Matrix value error";
  EXPECT_EQ(m1(1, 2), 2.) << "Matrix value error";
  EXPECT_EQ(m1(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(2, 1), 2.) << "Matrix value error";
  EXPECT_EQ(m1(2, 2), 2.) << "Matrix value error";
  // append at 1,1 with c = 2.0 and c0 = 0. (default).
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  4  4]
  //     [ 1  1  1]     [ 2  2]         [ 1  4  4]
  m1.Resize(3, 3, 1.);
  m2.Resize(2, 2., 2.);
  m1.Append(m2, 1, 1, 2.);
  EXPECT_EQ(m1(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m1(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m1(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(1, 1), 4.) << "Matrix value error";
  EXPECT_EQ(m1(1, 2), 4.) << "Matrix value error";
  EXPECT_EQ(m1(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(2, 1), 4.) << "Matrix value error";
  EXPECT_EQ(m1(2, 2), 4.) << "Matrix value error";
  // append at 1,1 with c = 2.0 and c0 = 2.
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  6  6]
  //     [ 1  1  1]     [ 2  2]         [ 1  6  6]
  m1.Resize(3, 3, 1.);
  m2.Resize(2, 2., 2.);
  m1.Append(m2, 1, 1, 2., 2.);
  EXPECT_EQ(m1(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m1(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m1(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(1, 1), 6.) << "Matrix value error";
  EXPECT_EQ(m1(1, 2), 6.) << "Matrix value error";
  EXPECT_EQ(m1(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m1(2, 1), 6.) << "Matrix value error";
  EXPECT_EQ(m1(2, 2), 6.) << "Matrix value error";
}

// inline Matrix& AppendRow(const Vector& v, int row, int col, double c = 1.0,
//                            double c0 = 0.) {
TEST_F(MatrixTest, AppendRow) {
  Matrix m;
  Vector v;
  // append at 1,1 with c = 1.0 (default) and c0 = 0. (default).
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  2  2]
  //     [ 1  1  1]                     [ 1  1  1]
  m.Resize(3, 3, 1.);
  v.Resize(2, 2.);
  m.AppendRow(v, 1, 1);
  EXPECT_EQ(m(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 1), 2.) << "Matrix value error";
  EXPECT_EQ(m(1, 2), 2.) << "Matrix value error";
  EXPECT_EQ(m(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 2), 1.) << "Matrix value error";
  // append at 1,1 with c = 2.0 and c0 = 0. (default).
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  4  4]
  //     [ 1  1  1]                     [ 1  1  1]
  m.Resize(3, 3, 1.);
  v.Resize(2, 2.);
  m.AppendRow(v, 1, 1, 2.);
  EXPECT_EQ(m(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 1), 4.) << "Matrix value error";
  EXPECT_EQ(m(1, 2), 4.) << "Matrix value error";
  EXPECT_EQ(m(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 2), 1.) << "Matrix value error";
  // append at 1,1 with c = 2.0 and c0 = 2.
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  6  6]
  //     [ 1  1  1]                     [ 1  1  1]
  m.Resize(3, 3, 1.);
  v.Resize(2, 2.);
  m.AppendRow(v, 1, 1, 2., 2.);
  EXPECT_EQ(m(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 1), 6.) << "Matrix value error";
  EXPECT_EQ(m(1, 2), 6.) << "Matrix value error";
  EXPECT_EQ(m(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 2), 1.) << "Matrix value error";
}

// inline Matrix& AppendCol(const Vector& v, int row, int col, double c = 1.0,
//                            double c0 = 0.) {
TEST_F(MatrixTest, AppendCol) {
  Matrix m;
  Vector v;
  // append at 1,1 with c = 1.0 (default) and c0 = 0. (default).
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  2  1]
  //     [ 1  1  1]                     [ 1  2  1]
  m.Resize(3, 3, 1.);
  v.Resize(2, 2.);
  m.AppendCol(v, 1, 1);
  EXPECT_EQ(m(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 1), 2.) << "Matrix value error";
  EXPECT_EQ(m(1, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 1), 2.) << "Matrix value error";
  EXPECT_EQ(m(2, 2), 1.) << "Matrix value error";
  // append at 1,1 with c = 2.0 and c0 = 0. (default).
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  4  1]
  //     [ 1  1  1]                     [ 1  4  1]
  m.Resize(3, 3, 1.);
  v.Resize(2, 2.);
  m.AppendCol(v, 1, 1, 2.);
  EXPECT_EQ(m(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 1), 4.) << "Matrix value error";
  EXPECT_EQ(m(1, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 1), 4.) << "Matrix value error";
  EXPECT_EQ(m(2, 2), 1.) << "Matrix value error";
  // append at 1,1 with c = 2.0 and c0 = 2.
  // Matrix should be:
  //     [ 1  1  1]                     [ 1  1  1]          
  // m = [ 1  1  1] v = [ 2  2] --> m = [ 1  6  1]
  //     [ 1  1  1]                     [ 1  6  1]
  m.Resize(3, 3, 1.);
  v.Resize(2, 2.);
  m.AppendCol(v, 1, 1, 2., 2.);
  EXPECT_EQ(m(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 1), 6.) << "Matrix value error";
  EXPECT_EQ(m(1, 2), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(2, 1), 6.) << "Matrix value error";
  EXPECT_EQ(m(2, 2), 1.) << "Matrix value error";
}

// void Solve(Vector& x, const Vector& b);
// friend Matrix Inverse(const Matrix& m);
// friend double Det(const Matrix& m);
// friend std::ostream& operator<<(std::ostream& s, const Matrix& m) {
// void add_BTCB(int row, int col, const int* perm, const Matrix& B1,

// inline Matrix Identity(int n) {
TEST_F(MatrixTest, Identity) {
  Matrix m = Identity(3);
  EXPECT_EQ(m(0, 0), 1.) << "Matrix value error";
  EXPECT_EQ(m(0, 1), 0.) << "Matrix value error";
  EXPECT_EQ(m(0, 2), 0.) << "Matrix value error";
  EXPECT_EQ(m(1, 0), 0.) << "Matrix value error";
  EXPECT_EQ(m(1, 1), 1.) << "Matrix value error";
  EXPECT_EQ(m(1, 2), 0.) << "Matrix value error";
  EXPECT_EQ(m(2, 0), 0.) << "Matrix value error";
  EXPECT_EQ(m(2, 1), 0.) << "Matrix value error";
  EXPECT_EQ(m(2, 2), 1.) << "Matrix value error";
}

// inline Matrix VVT(const Vector& v1, const Vector& v2) {


