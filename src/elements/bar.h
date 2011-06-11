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

#ifndef SRC_ELEMENTS_BAR_H_
#define SRC_ELEMENTS_BAR_H_

#include <vector>
#include "elements/element.h"

class CrossSection;
class UniaxialMaterial;

class Bar: public Element {
 public:
  // Constructors and Destructor
  Bar();
  Bar(int id,
      std::vector<Node*> nodes,
      UniaxialMaterial* material,
      CrossSection* iSec,
      CrossSection* jSec,
      int dim);
  ~Bar();

  void update();
  void commit();
  const Matrix& get_M();
  const Vector& get_Reff();
  void recoverStresses();
 protected:
  int dim_;
  CrossSection* iSection;
  CrossSection* jSection;
  double L0;
  double A0;
  Vector cosX;
  UniaxialMaterial* myUniMaterial;
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Bar(const Bar&);
  void operator=(const Bar&);
};
#endif  // SRC_ELEMENTS_BAR_H_
