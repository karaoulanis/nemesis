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

#ifndef SRC_ELEMENTS_BRICK8_H_
#define SRC_ELEMENTS_BRICK8_H_

#include <vector>
#include "elements/element.h"

class MatPoint;
class MultiaxialMaterial;

class Brick8: public Element {
 public:
  // Constructors and Destructor
  Brick8();
  Brick8(int id,
        std::vector<Node*> nodes,
        MultiaxialMaterial* material);
  virtual ~Brick8();

  virtual const Matrix& get_K();
  virtual const Matrix& get_M();
  virtual const Vector& get_R();

  virtual void Update();
  virtual void Commit();

  void AddInitialStresses(int direction, double h1, double s1,
                          double h2, double s2, double K0);
  void recoverStresses();
  int get_num_plastic_points();

  void shapeFunctions();
  virtual void get_B(Matrix* B, int node, int gPoint)=0;

 protected:
  std::vector<MatPoint*> myMatPoints;
  static double shp[8][4][8];
  static double detJ[8];
  static std::vector<int> perm;
};
#endif  // SRC_ELEMENTS_BRICK8_H_
