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

#include "model/model_element.h"

/**
 * Default constructor.
 */
ModelElement::ModelElement() {
}
/**
 * Simple constructor.
 * Initializes the ModelObject, which in turn initializes the FEObject and
 * passes the FTable to the ModelObject. This constructor will be used when
 * other objects than elements (such as constraints) will need to occupy a
 * ModelNode.
 */
ModelElement::ModelElement(const IDContainer& FTable, Element* pElement,
                           Constraint* pConstraint)
  :ModelObject(FTable), myElement(pElement), myConstraint(pConstraint) {
}
/**
 * Destructor.
 */
ModelElement::~ModelElement() {
}
void ModelElement::update() {
  if (myElement != 0) myElement->update();
}
void ModelElement::commit() {
  if (myElement != 0) myElement->commit();
}
