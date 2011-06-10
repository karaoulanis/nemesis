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

#ifndef SRC_MODEL_MODEL_ELEMENT_H_
#define SRC_MODEL_MODEL_ELEMENT_H_

#include "model/model_object.h"

class Constraint;
class ModelElement : public ModelObject {
 public:
  // Constructors
  ModelElement();
  ModelElement(const IDContainer& FTable, Element* pElement,
               Constraint* pConstraint);

  virtual ~ModelElement();

  // Access to data members
  Element* get_element()   {return myElement;}
  Constraint* get_constraint() {return myConstraint;}

  virtual void update();
  virtual void commit();

  virtual void add_K(double factor = 1.0)=0;
  virtual void add_M(double factor = 1.0)=0;
  virtual void add_C(double factor = 1.0)=0;
  virtual void add_R(double factor = 1.0)=0;
  virtual void add_Reff(double factor = 1.0)=0;

  virtual void add_Kgrad(double /*factor = 1.0*/) {}

 protected:
  Element* myElement;
  Constraint* myConstraint;

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  ModelElement(const ModelElement&);
  void operator=(const ModelElement&);
};

#endif  // SRC_MODEL_MODEL_ELEMENT_H_
