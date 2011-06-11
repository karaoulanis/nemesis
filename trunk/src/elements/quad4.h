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

#ifndef SRC_ELEMENTS_QUAD4_H_
#define SRC_ELEMENTS_QUAD4_H_

#include <vector>
#include "elements/element.h"

class MatPoint;
class MultiaxialMaterial;

class Quad4: public Element {
 public:
  // Constructors and Destructor
  Quad4();
  Quad4(int id,
        std::vector<Node*> nodes,
        MultiaxialMaterial* material,
        double thickness);
  ~Quad4();


  void findShapeFunctionsAt(MatPoint* pMatPoint);
  void AddInitialStresses(int direction, double h1, double s1,
                          double h2, double s2, double K0);
  void recoverStresses();
  void commit();
  int get_num_plastic_points();

 protected:
  int p1;  // Integration rule xi
  int p2;  // Integration rule eta
  double thickness_;
  std::vector<MatPoint*> myMatPoints;
  std::vector<MultiaxialMaterial*> materials;
  static Matrix N;
  static double detJ;
};
#endif  // SRC_ELEMENTS_QUAD4_H_
