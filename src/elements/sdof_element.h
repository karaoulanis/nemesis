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

#ifndef SRC_ELEMENTS_SDOF_ELEMENT_H_
#define SRC_ELEMENTS_SDOF_ELEMENT_H_

#include <vector>
#include "elements/element.h"
#include "material/sdof_material.h"

class SDofElement: public Element {
 public:
  // Constructors and Destructor
  SDofElement();
  SDofElement(int id,
      std::vector<Node*> nodes,
      int dof,
      SDofMaterial* material);
  ~SDofElement();

  const Matrix& get_K();
  const Matrix& get_M();
  const Vector& get_R();

  void Update(const double Dt=0.) {return;}
  void Commit() {return;}
 private:
  SDofMaterial* mySDofMaterial;
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  SDofElement(const SDofElement&);
  void operator=(const SDofElement&);
};

#endif  // SRC_ELEMENTS_SDOF_ELEMENT_H_
