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

#ifndef SRC_ELEMENTS_BRICK8B_H_
#define SRC_ELEMENTS_BRICK8B_H_

#include <vector>
#include "elements/brick8.h"

class Brick8b: public Brick8 {
 public:
  // Constructors and Destructor
  Brick8b();
  Brick8b(int id,
        std::vector<Node*> nodes,
        MultiaxialMaterial* material);
  ~Brick8b();

  void get_B(Matrix* B, int node, int gPoint);
};
#endif  // SRC_ELEMENTS_BRICK8B_H_
