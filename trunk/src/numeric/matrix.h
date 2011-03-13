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

#ifndef SRC_NUMERIC_MATRIX_H_
#define SRC_NUMERIC_MATRIX_H_

// Included files
#include "exception/sexception.h"
#include "numeric/array_check.h"
#include "numeric/vector.h"

// Forward declarations
class Vector;

class Matrix {
 private:
  int rows_, cols_, size_;
  double* data_;
 public:
  Matrix();
  Matrix(int m, int n);
  Matrix(int m, int n, double c);
  Matrix(const Matrix& m);
  ~Matrix();
  Matrix& operator=(const Matrix& m);
  /**
   * Implements (i, j) operator: m(i, j)
   * Range checking is provided for the debug version.
   * @param i The number of row.
   * @param j The number of column.
   */
  inline double& operator()(int i, int j) {
    #ifdef _DEBUG
    array_range_check(i, j, rows_, cols_);
    #endif
    return data_[i*cols_+j];
  }
  /**
   * Implements (i, j) operator: m(i, j)
   * Range checking is provided for the debug version.
   * @param i The number of row.
   * @param j The number of column.
   */
  inline double operator()(int i, int j) const {
    #ifdef _DEBUG
    array_range_check(i, j, rows_, cols_);
    #endif
    return data_[i*cols_+j];
  }
  /**
   * Returns the number of rows of the Matrix.
   */
  inline int rows() const {
    return rows_;
  }
  /**
   * Returns the number of columns of the Matrix.
   */
  inline int cols() const {
    return cols_;
  }
  /**
   * Returns a pointer to Matrix data.
   */
  inline double* data() {
    return data_;
  }
  /**
   * Returns a const pointer to Matrix data.
   */
  inline double* data() const {
    return data_;
  }
  void resize(int m, int n);
  void resize(int m, int n, double c);
  void clear();
  /**
   * Implements += operator: this+=c.
   * @param c All elements of the Matrix are added by this factor.
   * @return  A reference to the Matrix.
   */
  inline Matrix& operator+=(double c) {
    for (int i = 0;i < size_;i++) data_[i]+=c;
    return *this;
  }
  /**
   * Implements += operator: this-=c.
   * @param c This factor is substracted by all elements of the Matrix.
   * @return  A reference to the Matrix.
   */
  inline Matrix& operator-=(double c) {
    for (int i = 0;i < size_;i++) data_[i]-=c;
    return *this;
  }
  /**
   * Implements *= operator: this*=c.
   * @param c All elements of the Matrix are multiplied with this factor.
   * @return  A reference to the Matrix.
   */
  inline Matrix& operator*=(double c) {
    for (int i = 0;i < size_;i++) data_[i]*=c;
    return *this;
  }
    /**
   * Implements /= operator: this/=c.
   * @param c All elements of this Matrix are divided with this factor.
   * @return  A reference to the Matrix.
   */
  inline Matrix& operator/=(double c) {
    for (int i = 0;i < size_;i++) data_[i]/=c;
    return *this;
  }
  /**
   * Implements * operator: m = this*c.
   * @param c All elements of the Matrix are multiplied with this factor.
   * @return  A reference to the Matrix.
   */
  inline Matrix operator*(const double c) const {
    Matrix res(*this);
    res*=c;
    return res;
  }
  /**
   * Implements * operator: m = c*this, using member * operator.
   * @param c All elements of the Matrix are multiplied with this factor.
   * @param m A reference to the matrix that will be multiplied.
   * @return  Created Matrix.
   */
    inline friend Matrix operator*(const double c, const Matrix& m) {
    return m*c;
  }
  /**
   * Implements / operator: m = this/c.
   * @param c All elements of the Matrix are divided with this factor.
   * @return  Created Matrix.
   */
  inline Matrix operator/(const double c) const {
    Matrix res(*this);
    res/=c;
    return res;
  }
  /**
   * Implements c0*this+=c*m.
   * No size checking is provided to avoid double checking with +/- operators.
   * @param c  Factor for the matrix to be added.
   * @param m Matrix to be added.
   * @param c0 Factor for the existing Matrix (this).
   */
  inline void add_cM(double c, const Matrix& m, double c0 = 1.0) {
    if (c0 == 0.0)     for (int i = 0;i < size_;i++) data_[i]=0;
    else if (c0 != 1.0)  for (int i = 0;i < size_;i++) data_[i]*=c0;
    for (int i = 0;i < size_;i++) data_[i]+=c*m.data_[i];
  }
  /**
   * Implements c0*this+=c*m1*m2.
   * No size checking is provided to avoid double checking with * operator.
   * @param c  Factor for the matrix multiplication.
   * @param m1 First Matrix to be multiplied.
   * @param m2 Second Matrix to be multiplied.
   * @param c0 Factor for the existing Matrix (this).
   */
  inline void add_cMM(double c, const Matrix& m1, const Matrix& m2,
                      double c0 = 0.) {
    if (c0 == 0.0)     for (int i = 0;i < size_;i++) data_[i]=0;
    else if (c0 != 1.0)  for (int i = 0;i < size_;i++) data_[i]*=c0;
    int bCols = m1.cols();
    for (int i = 0;i < rows_;i++)
      for (int j = 0;j < cols_;j++)
        for (int k = 0;k < bCols;k++)
          data_[i*cols_+j]+=c*m1.data_[i*bCols+k]*m2.data_[k*cols_+j];
  }
  /**
   * Implements += operator: this+=m.
   * Size checking is provided for the debug version.
   * @param m The given Matrix.
   * @return  A reference to the Matrix.
   */
  inline Matrix& operator+=(const Matrix& m) {
    #ifdef _DEBUG
    array_size_check(rows_, cols_, m.rows_, m.cols_);
    #endif
    for (int i = 0;i < size_;i++) data_[i]+=m.data_[i];
    return *this;
  }
  /**
   * Implements -= operator: this-=m.
   * Size checking is provided for the debug version.
   * @param m The given Matrix.
   * @return  A reference to the Matrix.
   */
  inline Matrix& operator-=(const Matrix& m) {
    #ifdef _DEBUG
    array_size_check(rows_, cols_, m.rows_, m.cols_);
    #endif
    for (int i = 0;i < size_;i++) data_[i]-=m.data_[i];
    return *this;
  }
  /**
   * Implements + operator: this+m.
   * Uses c*this+=c1*m, with c0 = 1.0 and c = 1.0.
   * Size checking is provided for the debug version.
   * @param m The given Matrix.
   * @return  Created Matrix.
   */
  inline Matrix operator+(const Matrix& m) const {
    #ifdef _DEBUG
    array_size_check(rows_, cols_, m.rows_, m.cols_);
    #endif
    Matrix res(*this);
    res.add_cM(1.0, m, 1.0);
    return res;
  }
  /**
   * Implements - operator: this-m.
   * Uses c*this+=c1*m with c0 = 1.0 and c=-1.0.
   * Size checking is provided for the debug version.
   * @param m The given Matrix.
   * @return  Created Matrix.
   */
  inline Matrix operator-(const Matrix& m) const {
    #ifdef _DEBUG
    array_size_check(rows_, cols_, m.rows_, m.cols_);
    #endif
    Matrix res(*this);
    res.add_cM(-1.0, m, 1.0);
    return res;
  }
  /**
   * Implements * operator: m1*m2.
   * Uses c0*this+=c*m1*m2, with c0 = 0 c = 1.0 and m1 = this.
   * Size checking is provided for the debug version.
   * @param m The given Matrix.
   * @return  Created Matrix.
   */
  inline Matrix operator*(const Matrix& m) const {
    #ifdef _DEBUG
    array_size_check(cols_, m.rows_);
    #endif
    Matrix res(rows_, m.cols_);
    res.add_cMM(1.0, *this, m, 0.0);
    return res;
  }
  /**
   * Implements + operator: +this
   * @return  A reference to the Matrix.
   */
  inline Matrix& operator+() {
    return *this;
  }
  /**
   * Implements - operator: -this.
   * @return  Created Matrix.
   */
  inline Matrix operator-() {
    Matrix res(cols_, rows_);
    for (int i = 0;i < size_;i++) res.data_[i]=-data_[i];
    return res;
  }
  /**
   * Implements the Transpose of a Matrix: Transpose(m)
   * @param m The given Matrix.
   * @return  Created Matrix.
   */
  inline friend Matrix Transpose(const Matrix& m) {
    Matrix res(m.cols_, m.rows_);
    for (int i = 0;i < m.rows_;i++)
      for (int j = 0;j < m.cols_;j++)
        res.data_[j*m.rows_+i]=m.data_[i*m.cols_+j];
    return res;
  }
  /**
   * Implements Matrix*Vector.
   * Size checking is provided for the debug version.
   * @param v The given Vector.
   * @return  Created Vector.
   */
  inline Vector operator*(const Vector& v) const {
    #ifdef _DEBUG
    array_size_check(cols_, v.size());
    #endif
    Vector res(rows_, 0.);
    for (int i = 0;i < rows_;i++)
      for (int j = 0;j < cols_;j++) res[i]+=data_[i*cols_+j]*v[j];
    return res;
  }
  /**
   * Appends the entries of a Matrix m.
   * Range checking is provided for the debug version.
   * @param m   The Matrix to be copied.
   * @param row The row from which copy starts.
   * @param col The column from which copy starts.
   * @param c   A factor to be multiplied with the appended entries.
   * @param c0  A factor to be multiplied with the existing entries.
   */
  inline Matrix& append(const Matrix& m, int row, int col, double c = 1.0,
                        double c0 = 0.) {
    #ifdef _DEBUG
    array_range_check(row+m.rows_, col+m.cols_, rows_, cols_);
    #endif
    if (c0 == 0.0)
      for (int i = 0;i < m.rows_;i++)
        for (int j = 0;j < m.cols_;j++) data_[(row+i)*cols_+(col+j)]=0.;
    else if (c0 != 1.0)
      for (int i = 0;i < m.rows_;i++)
        for (int j = 0;j < m.cols_;j++) data_[(row+i)*cols_+(col+j)]*=c0;
    for (int i = 0;i < m.rows_;i++)
      for (int j = 0;j < m.cols_;j++)
        data_[(row+i)*cols_+(col+j)]+=c*m.data_[i*m.cols_+j];
    return *this;
  }
  /**
   * Appends the entries of a Vector vin a given row.
   * Range checking is provided for the debug version.
   * @param v   The Vector to be copied.
   * @param row The row from which copy starts.
   * @param col The column from which copy starts.
   * @param c   A factor to be multiplied with the appended entries.
   * @param c0  A factor to be multiplied with the existing entries.
   */
  inline Matrix& appendRow(const Vector& v, int row, int col, double c = 1.0,
                           double c0 = 0.) {
    #ifdef _DEBUG
    array_range_check(row, col+v.size(), rows_, cols_);
    #endif
    if (c0 == 0.0)
      for (int i = 0;i < v.size();i++) data_[row*cols_+col+i]=0.;
    else if (c0 != 1.0)
      for (int i = 0;i < v.size();i++) data_[row*cols_+col+i]*=c0;
    for (int i = 0;i < v.size();i++)     data_[row*cols_+col+i]+=c*v[i];
    return *this;
  }
  /**
   * Appends the entries of a Vector vin a given column.
   * Range checking is provided for the debug version.
   * @param v   The Vector to be copied.
   * @param row The row from which copy starts.
   * @param col The column from which copy starts.
   * @param c   A factor to be multiplied with the appended entries.
   * @param c0  A factor to be multiplied with the existing entries.
   */
  inline Matrix& appendCol(const Vector& v, int row, int col, double c = 1.0,
                           double c0 = 0.) {
    #ifdef _DEBUG
    array_range_check(row+v.size(), col, rows_, cols_);
    #endif
    if (c0 == 0.0)
      for (int i = 0;i < v.size();i++) data_[(row+i)*cols_+col]=0.;
    else if (c0 != 1.0)
      for (int i = 0;i < v.size();i++) data_[(row+i)*cols_+col]*=c0;
    for (int i = 0;i < v.size();i++)     data_[(row+i)*cols_+col]+=c*v[i];
    return *this;
  }

  void solve(Vector& x, const Vector& b);
  friend Matrix Inverse(const Matrix& m);
  friend double det(const Matrix& m);

  friend std::ostream& operator<<(std::ostream& s, const Matrix& m) {
    s << 1200 << ' ' << m.rows_ << ' ' << m.cols_ << ' ';
    for (int i = 0;i < m.size_;i++) s << m.data_[i] << ' ';
    return s;
  }
  void add_BTCB(int row, int col, const int* perm, const Matrix& B1,
    const Matrix& C, const Matrix B2, double c1, double c0 = 0.) {
    int m = B1.rows();
    int n = B2.cols();
    int pos = row*cols_+col;
    int colC = C.cols();
    double* pB1 = B1.data();
    double* pB2 = B2.data();
    double* pC = C.data();
    for (int i = 0;i < n;i++)
      for (int j = 0;j < n;j++)
        data_[pos+i*cols_+j]*=c0;
    for (int i = 0;i < n;i++)
      for (int k = 0;k < m;k++)
        for (int l = 0;l < m;l++)
          for (int j = 0;j < n;j++)
            data_[pos+i*cols_+j]+=c1*pB1[k*n+i]*pC[perm[k]*colC+perm[l]]*pB2[l*n+j];
  }
};
/**
* Returns an identity matrix of size n.
*/
inline Matrix Identity(int n) {
  Matrix res(n, n, 0.);
  for (int i = 0;i < n;i++) res(i, i)=1.;
  return res;
}
/**
 * Implements Vector1*Transpose(Vector).
 * Size checking is provided for the debug version.
 * It is implemented here to avoid conflicts.
 * A temporary object is created.
 * @param v1 The column vector.
 * @param v2 The row (transposed) vector.
 * @return   Created Matrix.
 */
inline Matrix VVT(const Vector& v1, const Vector& v2) {
  #ifdef _DEBUG
  array_size_check(v1.size(), v2.size());
  #endif
  Matrix res(v1.size(), v2.size());
  for (int i = 0;i < v1.size();i++)
    for (int j = 0;j < v2.size();j++)
      res(i, j)=v1[i]*v2[j];
  return res;
}
#endif  // SRC_NUMERIC_MATRIX_H_
