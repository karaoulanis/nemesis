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

#ifndef SRC_LOADCASE_ELEMENTAL_LOAD_H_
#define SRC_LOADCASE_ELEMENTAL_LOAD_H_

#include "elements/element.h"
#include "loadcase/load.h"
#include "numeric/vector.h"

class Element;

class ElementalLoad: public Load {
 protected:
  Element* myElement;   /// A pointer to the corresponding Element.
  static Vector P;
  public:
  // Constructor and destructors
  ElementalLoad();
  ElementalLoad(int elemID);
  virtual ~ElementalLoad();
/*
  // Access to member data
  int set_load_direction(LoadDirection direction);
  int set_user_direction(Vector* direction);
  void set_A(Vector& aValues);
  void set_P(Vector& pValues);

  LoadDirection get_load_direction();
  Vector* get_user_direction();
  const Vector& get_A();
  const Vector& get_P();
*/  
  // Apply load
  virtual const Vector& get_P()=0;
  void apply(double fact, double t);
};


#endif  // SRC_LOADCASE_ELEMENTAL_LOAD_H_
