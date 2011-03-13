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

#ifndef SRC_SOE_SOE_H_
#define SRC_SOE_SOE_H_

// C++ system files
#include <cmath>
#include <iostream>
#include <valarray>

// Project files (alphabetically)
#include "analysis/analysis.h"
#include "analysis/analysis_object.h"
#include "containers/containers.h"
#include "numeric/lapack.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

// Forward declarations
class Analysis;
class Model;

class SOE: public AnalysisObject {
 protected:
  int size_;
  std::valarray<double> A;
  Vector X;
  Vector B;
  bool isLUFactored;
  std::valarray<int> IPIV;
  public:
  SOE();
  virtual ~SOE();

  virtual int insertMatrixIntoA(const Matrix& Ke, const IDContainer& EFTable,
                  double factor = 1.0)=0;
  virtual int insertVectorIntoB(const Vector& Ve, const IDContainer& EFTable,
                  double factor = 1.0);
  virtual int insertMatrixIntoA(const Matrix& Be, const IDContainer& EFTable,
                  const IDContainer& SFTable, double factor = 1.0);

  virtual void zeroA();
  virtual void zeroB();
  virtual void zeroX();
  virtual void zero();

  const Vector& get_X();
  const Vector& get_B();
  void addB(const Vector& v);
  void set_B(const Vector& v);
  void set_X(const Vector& v);

  virtual void set_size()=0;
  virtual int get_eigen_sign()=0;
  virtual void print()=0;
  virtual void printSolution();

  virtual int solve()=0;
  int plotGraph(const char* s);
};
#endif  // SRC_SOE_SOE_H_
