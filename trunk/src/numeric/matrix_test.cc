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
  EXPECT_EQ(m.get_rows(), 0) << "Matrix rows is not null initialized";
  EXPECT_EQ(m.get_cols(), 0) << "Matrix cols is not null initialized";
  // ASSERT_EQ(m.get_size(),0)    /// @todo: implement
  // EXPECT_EQ(m.data(),NULL)     /// @todo: implement;
}

// test: Matrix::Matrix(int rows, int cols);
TEST_F(MatrixTest, SimpleConstructor) {
  Matrix m(3, 2);
  EXPECT_EQ(m.get_rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m.get_cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)  /// @todo: implement;
}

// test: Matrix::Matrix(int rows, int cols, double value);
TEST_F(MatrixTest, ConstructorWithDefaultValues) {
  Matrix m(3, 2, 2.);
  EXPECT_EQ(m.get_rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m.get_cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
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
  EXPECT_EQ(m2.get_rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m2.get_cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)  /// @todo: implement;
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
  EXPECT_EQ(m2.get_rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m2.get_cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m2(i, j), 2.) << "Matrix value error";
    }
  }
  Matrix m3(4, 4, 2.);
  m3 = m1;
  EXPECT_EQ(m3.get_rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m3.get_cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m3(i, j), 2.) << "Matrix value error";
    }
  }
}

// inline double& operator()(int i, int j) {
// inline double operator()(int i, int j) const {
// inline int get_rows() const {
// inline int get_cols() const {
// inline double* get_data() {
// inline double* get_data() const {
// inline double* get_data() const {
// void Resize(int rows, int cols);
// void Resize(int cols, int rows, double value);
// void Clear();
// inline Matrix& operator+=(double value) {
// inline Matrix& operator-=(double value) {
// inline Matrix& operator*=(double value) {
// inline Matrix& operator/=(double value) {
// inline Matrix operator*(const double value) const {
// inline friend Matrix operator*(const double c, const Matrix& m) {
// inline Matrix operator/(const double c) const {
// inline void add_cM(double c, const Matrix& m, double c0 = 1.0) {
// inline void add_cMM(double c, const Matrix& m1, const Matrix& m2,
// inline Matrix& operator+=(const Matrix& m) {
// inline Matrix& operator-=(const Matrix& m) {
// inline Matrix operator+(const Matrix& m) const {
// inline Matrix operator-(const Matrix& m) const {
// inline Matrix operator*(const Matrix& m) const {
// inline Matrix& operator+() {
// inline Matrix operator-() {
// inline friend Matrix Transpose(const Matrix& m) {
// inline Vector operator*(const Vector& v) const {
// inline Matrix& Append(const Matrix& m, int row, int col, double c = 1.0,
//                         double c0 = 0.) {
// inline Matrix& AppendRow(const Vector& v, int row, int col, double c = 1.0,
//                            double c0 = 0.) {
// inline Matrix& AppendCol(const Vector& v, int row, int col, double c = 1.0,
//                            double c0 = 0.) {
// void Solve(Vector& x, const Vector& b);
// friend Matrix Inverse(const Matrix& m);
// friend double det(const Matrix& m);
// friend std::ostream& operator<<(std::ostream& s, const Matrix& m) {
// void add_BTCB(int row, int col, const int* perm, const Matrix& B1,
// inline Matrix Identity(int n) {
// inline Matrix VVT(const Vector& v1, const Vector& v2) {

// test: void Matirx::resize(int m, int n)
TEST_F(MatrixTest, MatrixResize) {
  Matrix m(3, 2);
  m.resize(2, 3);
  EXPECT_EQ(m.get_rows(), 2) << "Matrix rows not set correctly";
  EXPECT_EQ(m.get_cols(), 3) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
}

// test: void Matirx::resize(int m, int n double c)
TEST_F(MatrixTest, MatrixResizeWithDefaultValues) {
  Matrix m(2, 3);
  m.resize(3, 2, 2.);
  EXPECT_EQ(m.get_rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m.get_cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m(i, j), 2.) << "Matrix value error";
    }
  }
}

// test: void Matrix::clear()
TEST_F(MatrixTest, MatrixClear) {
  Matrix m(2, 3, 2.);
  m.clear();
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(m(i, j), 0.) << "Matrix value error";
    }
  }
}

