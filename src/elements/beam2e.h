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

#ifndef SRC_ELEMENTS_BEAM2E_H_
#define SRC_ELEMENTS_BEAM2E_H_

#include <vector>
#include "elements/element.h"

class CrossSection;
class UniaxialMaterial;

class Beam2e: public Element {
 public:
  Beam2e();
  Beam2e(int id,
         std::vector<Node*> nodes,
         UniaxialMaterial* material,
         CrossSection* section);
  ~Beam2e();

  const Matrix& get_K();
    const Matrix& get_M();
  const Vector& get_R();
  const Vector& get_Rgrad();

  bool checkIfAllows(FEObject* /*f*/)   {return true;}
  void update()           {return;}
  void commit()           {return;}
  void recoverStresses();

 protected:
  UniaxialMaterial* myUniMaterial;
  CrossSection* mySection;
  double L;
  double cosX[2];

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  Beam2e(const Beam2e&);
  void operator=(const Beam2e&);
};
#endif  // SRC_ELEMENTS_BEAM2E_H_
