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
#include "numeric/matrix.h"
#include "exception/sexception.h"
#include "numeric/array_check.h"
#include "numeric/lu.h"
#include "numeric/vector.h"

/**
 * Default constructor.
 * Initializes everything to zero.
 */
Matrix::Matrix()
    : rows_(0),
      cols_(0),
      size_(0),
      data_(0) {
}

/**
 * Constructor.
 * Allocates a size m*n double.
 * Exception handling for bad allocation is provided.
 * @param n The number of rows.
 * @param m The number of columns.
 */
Matrix::Matrix(int m, int n)
    : rows_(m),
      cols_(n),
      size_(m*n) {
  try {
    data_ = new double[size_];
  } catch(std::bad_alloc) {
    throw SException("[nemesis:%d] %s", 1001, "Run out of memory.\n");
  }
}

/**
* Constructor.
* Allocates a size m*n double and initializes entries to c.
* Exception handling for bad allocation is provided.
* @param n The number of rows.
* @param m The number of columns.
* @param c Initial value for all entries.
*/
Matrix::Matrix(int m, int n, double c)
    : rows_(m),
      cols_(n),
      size_(m*n) {
  try {
      data_ = new double[size_];
    } catch(std::bad_alloc) {
      throw SException("[nemesis:%d] %s", 1001, "Run out of memory.\n");
    }
    for (int i = 0; i < size_; i++) data_[i]=c;
  }

/**
 * Copy constructor.
 * Allocates a size m.rows_*m.cols_ double and copies entries from Matrix m.
 * Exception handling for bad allocation is provided.
 * @param m The Matrix that is copied.
 */
Matrix::Matrix(const Matrix& m)
    : rows_(m.rows_),
      cols_(m.cols_),
      size_(m.size_) {
  if (size_ > 0) {
    try {
      data_ = new double[size_];
    } catch(std::bad_alloc) {
      throw SException("[nemesis:%d] %s", 1001, "Run out of memory.\n");
    }
    for (int i = 0; i < size_; i++) data_[i]=m.data_[i];
  } else {
      data_ = 0;
  }
}

/**
 * Destructor.
 */
Matrix::~Matrix() {
  if (data_ != 0) delete[] data_;
}

/**
 * Copy assignment.
 * Size checking is provided for the debug version.
 * @param m The Matrix that is copied.
 * @return  A reference to the Matrix.
 */
Matrix& Matrix::operator=(const Matrix& m) {
  #ifdef _DEBUG
  array_size_check(rows_, cols_, m.rows_, m.cols_);
  #endif
  if (this != &m)
    for (int i = 0; i < size_; i++) data_[i]=m.data_[i];
  return *this;
}

/**
* Resizes the Matrix.
* Does not preserves data and does not initialize entries.
* If new Matrix of the same size no allocation takes place.
* Exception handling for bad allocation is provided.
* @param m The number of rows.
* @param n The number of columns.
*/
void Matrix::resize(int m, int n) {
  rows_ = m;
  cols_ = n;
  if (size_ != m*n) {
    size_ = m*n;
    if (data_ != NULL) delete[] data_;
    try {
      data_ = new double[size_];
    } catch(std::bad_alloc) {
      throw SException("[nemesis:%d] %s", 1001, "Run out of memory.\n");
    }
  }
}

/**
* Resizes the Matrix.
* Does not preserves data but initializes all entries to c.
* If new Matrix of the same size no allocation takes place.
* Exception handling for bad allocation is provided.
* @param m The number of rows.
* @param n The number of columns.
* @param c Initial value for all entries.
*/
void Matrix::resize(int m, int n, double c) {
  rows_ = m;
  cols_ = n;
  if (size_ != m*n) {
    size_ = m*n;
    if (data_ != NULL) delete[] data_;
    try {
      data_ = new double[size_];
    } catch(std::bad_alloc) {
      throw SException("[nemesis:%d] %s", 1001, "Run out of memory.\n");
    }
  }
  for (int i = 0; i < size_; i++) data_[i]=c;
}

/**
* Clears the contents of a Matrix.
*/
void Matrix::clear() {
  for (int i = 0; i < size_; i++) data_[i]=0.;
}

