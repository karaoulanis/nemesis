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

#ifndef SRC_ELEMENTS_TIMOSHENKO2D_H_
#define SRC_ELEMENTS_TIMOSHENKO2D_H_

#include <vector>
#include "elements/element.h"

class Node;
class UniaxialMaterial;
class CrossSection;

class Timoshenko2d: public Element {
 public:
  Timoshenko2d();
  Timoshenko2d(int id,
         std::vector<Node*> nodes,
         UniaxialMaterial* material,
         CrossSection* section,
         int rule);
  ~Timoshenko2d();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();
  const Vector& get_Rgrad();

  void update()                         {return;}
  void commit()                         {return;}
  void recoverStresses();
  void shapeFunctions(int n, double xi, double* N, double* dN);

 protected:
  UniaxialMaterial* myUniMaterial;
  CrossSection* mySection;
  double L;
  double cosX[2];
  int gPoints;
  static const double GaussCoords[4][4];
  static const double GaussWeights[4][4];

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Timoshenko2d(const Timoshenko2d&);
  void operator=(const Timoshenko2d&);
};
#endif  // SRC_ELEMENTS_TIMOSHENKO2D_H_
