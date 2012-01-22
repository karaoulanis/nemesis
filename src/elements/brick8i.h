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

#ifndef SRC_ELEMENTS_BRICK8I_H_
#define SRC_ELEMENTS_BRICK8I_H_

#include <vector>
#include "elements/brick8.h"

class Brick8i: public Brick8 {
 public:
  // Constructors and Destructor
  Brick8i();
  Brick8i(int id,
        std::vector<Node*> nodes,
        MultiaxialMaterial* material);
  ~Brick8i();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();
  void Update();
  void Commit();
  void get_B(Matrix* /*B*/, int /*node*/, int /*gPoint*/) {}
 private:
  static double shpStd[8][4][8];
  static double shpInc[3][4][8];
  static double detJ[8];

  Vector aTrial;
  Vector aConvg;
  void shapeFunctions();
  void get_Bstd(Matrix* B, int node, int gPoint);
  void get_BInc(Matrix* B, int node, int gPoint);
  void get_Kdd(Matrix* K);
  void get_Kda(Matrix* K);
  void get_Kaa(Matrix* K);
};

#endif  // SRC_ELEMENTS_BRICK8I_H_
