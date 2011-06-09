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

#ifndef SRC_ELEMENTS_TETRAHEDRON4_DISP_H_
#define SRC_ELEMENTS_TETRAHEDRON4_DISP_H_

#include <vector>
#include "elements/element.h"
#include "numeric/matrix.h"

class MatPoint;
class MultiaxialMaterial;

class Tetrahedron4Disp: public Element {
 public:
  // Constructors and Destructor
  Tetrahedron4Disp();
  Tetrahedron4Disp(int id,
        std::vector<Node*> nodes,
        MultiaxialMaterial* material);
  ~Tetrahedron4Disp();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();

  void update();
  void commit();

  void findShapeFunctions();
  bool checkIfAllows(FEObject* f);
  void recoverStresses();
 protected:
  std::vector<MatPoint*> myMatPoints;
  static Matrix N;
  static double V;
};
#endif  // SRC_ELEMENTS_TETRAHEDRON4_DISP_H_
