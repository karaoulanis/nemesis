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

#ifndef SRC_ELEMENTS_QUAD4_DISP_AXISYMMETRIC_H_
#define SRC_ELEMENTS_QUAD4_DISP_AXISYMMETRIC_H_

#include "elements/quad4.h"

class Quad4DispAxisymmetric: public Quad4 {
  public:
  // Constructors and Destructor
  Quad4DispAxisymmetric();
  Quad4DispAxisymmetric(int ID, int Node_1, int Node_2, int Node_3, int Node_4,
                        int MatID,
                        int integrationRuleXi, int integrationRuleEta);
  ~Quad4DispAxisymmetric();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();

  void update();
};

#endif  // SRC_ELEMENTS_QUAD4_DISP_AXISYMMETRIC_H_
