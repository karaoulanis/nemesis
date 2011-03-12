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

#ifndef SRC_MODEL_ELIMINATION_MODEL_ELEMENT_H_
#define SRC_MODEL_ELIMINATION_MODEL_ELEMENT_H_

#include "model/model_element.h"

class EliminationModelElement : public ModelElement {
 private:
 public:
  EliminationModelElement();
  EliminationModelElement(const IDContainer& FTable, Element* pElement);
  ~EliminationModelElement();

  void add_K(double factor = 1.0);
  void add_M(double factor = 1.0);
  void add_C(double factor = 1.0);
  void add_R(double factor = 1.0);
  void add_Reff(double factor = 1.0);
  void add_Rgrad(double factor = 1.0);
};
#endif  // SRC_MODEL_ELIMINATION_MODEL_ELEMENT_H_
