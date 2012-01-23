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

#ifndef SRC_ELEMENTS_QUAD4E_H_
#define SRC_ELEMENTS_QUAD4E_H_

#include <vector>
#include "elements/quad4.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

class Quad4e: public Quad4 {
 public:
  // Constructors and Destructor
  Quad4e();
  Quad4e(int id,
         std::vector<Node*> nodes,
         MultiaxialMaterial* material,
         double thickness);
  ~Quad4e();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();

  void Update();
  void formKR();
  double get_J(double xi, double eta);
  void formTo(double xi, double eta);
  void formNu(double xi, double eta);
  void formBu(double xi, double eta);
  void formBe(double xi, double eta);
  void formT0(double xi, double eta);

 private:
  static Matrix Bu;
  static Matrix Be;
  static Matrix T0;
  static Vector Nu;
  Vector alpha;
};
#endif  // SRC_ELEMENTS_QUAD4E_H_
