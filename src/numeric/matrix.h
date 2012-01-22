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

#ifndef SRC_NUMERIC_MATRIX_H_
#define SRC_NUMERIC_MATRIX_H_

// Included files
#include "exception/sexception.h"
#include "numeric/numeric.h"
#include "numeric/vector.h"

// Forward declarations
class Vector;

class Matrix {
 public:
  /**
   * Default constructor
   */
  Matrix();

  /**
   * Constructor.
   * Matrix is not initialized.
   * @param rows Number of rows.
   * @param cols Number of columns.
   */
  Matrix(int rows, int cols);

  /**
   * Constructor
   * @param rows Number of rows.
   * @param cols Number of columns.
   * @param value Value to be set to all entries.
   */
  Matrix(int rows, int cols, double value);

  /**
   * Copy constructor
   * @param matrix The matrix to be copied from.
   */
  Matrix(const Matrix& matrix);

  /**
   * Destructor
   */
  ~Matrix();

  /**
   * Copy assignment.
   * @param matrix The matrix to be assigned.
   */
  Matrix& operator=(const Matrix& m);

  /**
   * Implements (i, j) operator: m(i, j)
   * Range checking is provided for the debug version.
   * @param i The number of row.
   * @param j The number of column.
   */
  inline double& operator()(int i, int j) {
    #ifdef _DEBUG
    num::check::array_range(i, j, rows_, cols_);
    #endif
    return data_[i*cols_+j];
  }

  /**
   * Implements constant (i, j) operator: m(i, j)
   * Range checking is provided for the debug version.
   * @param i The number of row.
   * @param j The number of column.
   */
  inline double operator()(int i, int j) const {
    #ifdef _DEBUG
    num::check::array_range(i, j, rows_, cols_);
    #endif
    return data_[i*cols_+j];
  }

  /**
   * Returns the number of rows of the Matrix.
   */
  inline int get_rows() const {
    return rows_;
  }

  /**
   * Returns the number of columns of the Matrix.
   */
  inline int get_cols() const {
    return cols_;
  }

  /**
   * Returns the array size.
   */
  inline int get_size() const {
    return size_;
  }

  /**
   * Returns a pointer to Matrix data.
   */
  inline double* get_data() {
    return data_;
  }

  /**
   * Returns a const pointer to Matrix data.
   */
  inline double* get_data() const {
    return data_;
  }

  /**
  * Resizes the Matrix.
  * Does not preserve data and does not initialize entries.
  * If new Matrix of the same size no allocation takes place.
  * Throws an exception for bad allocation.
  * @param rows The number of rows.
  * @param cols The number of columns.
  */
  void Resize(int rows, int cols);

  /**
  * Resizes the Matrix and initializes all entries to value.
  * Does not preserve data and does not initialize entries.
  * If new Matrix of the same size no allocation takes place.
  * Throws an exception for bad allocation.
  * @param rows The number of rows.
  * @param cols The number of columns.
  * @param value Value to be set to all entries.
  */
  void Resize(int cols, int rows, double value);

  /**
   * Set all entries of the Matrix equal to zero.
   */
  void Clear();

  /**
   * Implements += operator: this+=value.
   * @param value Added to all elements of the Matrix.
   * @return A reference to the Matrix.
   */
  inline Matrix& operator+=(double value) {
    for (int i = 0; i < size_; i++) {
      data_[i] += value;
    }
    return *this;
  }

  /**
   * Implements -= operator: this-=value.
   * @param value Subtracted from all elements of the Matrix.
   * @return A reference to the Matrix.
   */
  inline Matrix& operator-=(double value) {
    for (int i = 0; i < size_; i++) {
      data_[i] -= value;
    }
    return *this;
  }

  /**
   * Implements *= operator: this*=value.
   * @param value Multiplied to all elements of the Matrix.
   * @return A reference to the Matrix.
   */
  inline Matrix& operator*=(double value) {
    for (int i = 0; i < size_; i++) {
      data_[i] *= value;
    }
    return *this;
  }

  /**
   * Implements /= operator: this/=value.
   * @param value All elements of the Matrix are devided by this value.
   * @return  A reference to the Matrix.
   */
  inline Matrix& operator/=(double value) {
    for (int i = 0; i < size_; i++) {
      data_[i] /= value;
    }
    return *this;
  }

  /**
   * Implements * operator: m = this*value.
   * @param value Multiplied to all elements of the Matrix.
   * @return A newly created Matrix.
   * @todo Highly inefficient, do not create object.
   */
  inline Matrix operator*(const double value) const {
    Matrix res(*this);
    res *= value;
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
  inline void Add_cM(double c, const Matrix& m, double c0 = 1.0) {
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
  inline void Add_cMM(double c, const Matrix& m1, const Matrix& m2,
                      double c0 = 0.) {
    if (c0 == 0.0)     for (int i = 0;i < size_;i++) data_[i]=0;
    else if (c0 != 1.0)  for (int i = 0;i < size_;i++) data_[i]*=c0;
    int bCols = m1.get_cols();
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
    num::check::array_size(rows_, cols_, m.rows_, m.cols_);
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
    num::check::array_size(rows_, cols_, m.rows_, m.cols_);
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
    num::check::array_size(rows_, cols_, m.rows_, m.cols_);
    #endif
    Matrix res(*this);
    res.Add_cM(1.0, m, 1.0);
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
    num::check::array_size(rows_, cols_, m.rows_, m.cols_);
    #endif
    Matrix res(*this);
    res.Add_cM(-1.0, m, 1.0);
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
    num::check::array_size(cols_, m.rows_);
    #endif
    Matrix res(rows_, m.cols_);
    res.Add_cMM(1.0, *this, m, 0.0);
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
    num::check::array_size(cols_, v.size());
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
  inline Matrix& Append(const Matrix& m, int row, int col, double c = 1.0,
                        double c0 = 0.) {
    #ifdef _DEBUG
    num::check::array_range(row+m.rows_, col+m.cols_, rows_, cols_);
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
  inline Matrix& AppendRow(const Vector& v, int row, int col, double c = 1.0,
                           double c0 = 0.) {
    #ifdef _DEBUG
    num::check::array_range(row, col+v.size(), rows_, cols_);
    #endif
    if (c0 == 0.0)
      for (int i = 0; i < v.get_size(); i++)  data_[row*cols_+col+i]=0.;
    else if (c0 != 1.0)
      for (int i = 0; i < v.get_size(); i++)  data_[row*cols_+col+i]*=c0;
    for (int i = 0; i < v.get_size(); i++)    data_[row*cols_+col+i]+=c*v[i];
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
  inline Matrix& AppendCol(const Vector& v, int row, int col, double c = 1.0,
                           double c0 = 0.) {
    #ifdef _DEBUG
    num::check::array_range(row+v.size(), col, rows_, cols_);
    #endif
    if (c0 == 0.0)
      for (int i = 0; i < v.get_size(); i++) data_[(row+i)*cols_+col]=0.;
    else if (c0 != 1.0)
      for (int i = 0; i < v.get_size(); i++) data_[(row+i)*cols_+col]*=c0;
    for (int i = 0; i < v.get_size(); i++)   data_[(row+i)*cols_+col]+=c*v[i];
    return *this;
  }

  /**
   *
   */
  void Solve(Vector& x, const Vector& b);

  /**
   *
   */
  friend Matrix Inverse(const Matrix& m);

  /**
   *
   */
  friend double Det(const Matrix& m);

  /**
   *
   */
  void Add_BTCB(int row, int col, const int* perm, const Matrix& B1,
    const Matrix& C, const Matrix B2, double c1, double c0 = 0.) {
    int m = B1.get_rows();
    int n = B2.get_cols();
    int pos = row*cols_+col;
    int colC = C.get_cols();
    double* pB1 = B1.get_data();
    double* pB2 = B2.get_data();
    double* pC = C.get_data();
    for (int i = 0;i < n;i++)
      for (int j = 0;j < n;j++)
        data_[pos+i*cols_+j]*=c0;
    for (int i = 0;i < n;i++)
      for (int k = 0;k < m;k++)
        for (int l = 0;l < m;l++)
          for (int j = 0;j < n;j++)
            data_[pos+i*cols_+j]+=
              c1*pB1[k*n+i]*pC[perm[k]*colC+perm[l]]*pB2[l*n+j];
  }
 private:
  int rows_;
  int cols_;
  int size_;
  double* data_;
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
  num::check::array_size(v1.get_size(), v2.get_size());
  #endif
  Matrix res(v1.get_size(), v2.get_size());
  for (int i = 0;i < v1.get_size();i++)
    for (int j = 0;j < v2.get_size();j++)
      res(i, j)=v1[i]*v2[j];
  return res;
}
#endif  // SRC_NUMERIC_MATRIX_H_
