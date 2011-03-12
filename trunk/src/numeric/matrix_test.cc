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
  EXPECT_EQ(m.rows(), 0) << "Matrix rows is not null initialized";
  EXPECT_EQ(m.cols(), 0) << "Matrix cols is not null initialized";
  // ASSERT_EQ(m.get_size(),0) /// @todo: implement
  // EXPECT_EQ(m.data(),NULL)    /// @todo: implement;
}

// test: Matrix::Matrix(int m, int n);
TEST_F(MatrixTest, SimpleConstructor) {
  Matrix m(3, 2);
  EXPECT_EQ(m.rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m.cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
}

// test: Matrix::Matrix(int m, int n, double c);
TEST_F(MatrixTest, ConstructorWithDefaultValues) {
  Matrix m(3, 2, 2.);
  EXPECT_EQ(m.rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m.cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m(i, j), 2.) << "Matrix value error";
    }
  }
}

// test: Matrix::Matrix(const Matrix& m);
TEST_F(MatrixTest, CopyConstructor) {
  Matrix m1(3, 2, 2.);
  Matrix m2(m1);
  EXPECT_EQ(m2.rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m2.cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
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
  EXPECT_EQ(m2.rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m2.cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m2(i, j), 2.) << "Matrix value error";
    }
  }
  Matrix m3(4, 4, 2.);
  m3 = m1;
  EXPECT_EQ(m3.rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m3.cols(), 2) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m3(i, j), 2.) << "Matrix value error";
    }
  }
}

// test: void Matirx::resize(int m, int n)
TEST_F(MatrixTest, MatrixResize) {
  Matrix m(3, 2);
  m.resize(2, 3);
  EXPECT_EQ(m.rows(), 2) << "Matrix rows not set correctly";
  EXPECT_EQ(m.cols(), 3) << "Matrix cols not set correctly";
  // EXPECT_EQ(m.get_size(),6)
}

// test: void Matirx::resize(int m, int n double c)
TEST_F(MatrixTest, MatrixResizeWithDefaultValues) {
  Matrix m(2, 3);
  m.resize(3, 2, 2.);
  EXPECT_EQ(m.rows(), 3) << "Matrix rows not set correctly";
  EXPECT_EQ(m.cols(), 2) << "Matrix cols not set correctly";
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
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(m(i, j), 0.) << "Matrix value error";
    }
  }
}

