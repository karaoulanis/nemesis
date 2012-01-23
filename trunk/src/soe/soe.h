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

#ifndef SRC_SOE_SOE_H_
#define SRC_SOE_SOE_H_

// C++ system files
#include <cmath>
#include <valarray>

// Project files (alphabetically)
#include "analysis/analysis_object.h"
#include "containers/containers.h"
#include "numeric/vector.h"

// Forward declarations
class Model;
class Matrix;

class SOE: public AnalysisObject {
 public:
  SOE();
  explicit SOE(Model* model);
  virtual ~SOE();

  virtual int InsertMatrixIntoA(const Matrix& Ke, const IDContainer& EFTable,
                  double factor = 1.0)=0;
  virtual int insertVectorIntoB(const Vector& Ve, const IDContainer& EFTable,
                  double factor = 1.0);
  virtual int InsertMatrixIntoA(const Matrix& Be, const IDContainer& EFTable,
                  const IDContainer& SFTable, double factor = 1.0);

  virtual void zeroA();
  virtual void zeroB();
  virtual void zeroX();
  virtual void zero();

  const Vector& get_X();
  const Vector& get_B();
  void addB(const Vector& v);
  void set_model(Model* model);
  void set_B(const Vector& v);
  void set_X(const Vector& v);

  virtual void set_size()=0;
  virtual int get_eigen_sign()=0;
  virtual void print()=0;

  virtual int solve()=0;

 protected:
  Model* model_;
  int size_;
  std::valarray<double> A;
  Vector X;
  Vector B;
  std::valarray<int> IPIV;
  bool isLUFactored;
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  SOE(const SOE&);
  void operator=(const SOE&);
};
#endif  // SRC_SOE_SOE_H_
