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


#include "model/standard_model_element.h"

/**
 * Default constructor.
 * Initializes only the ModelObject, which in turn initializes the FEObject.
 */
StandardModelElement::StandardModelElement()
  :ModelElement() {
}
/**
 * Constructor.
 * Initializes the ModelObject, which in turn initializes the FEObject, passes
 * the FTable to the ModelObject and copies the address of it's element.
 */
StandardModelElement::StandardModelElement(const IDContainer& FTable,
                                           Element* pElement)
  :ModelElement(FTable, pElement, 0) {
  myMatrix = theStaticMatrices[FTable.size()];
  myVector = theStaticVectors[FTable.size()];
}
StandardModelElement::~StandardModelElement() {
}
void StandardModelElement::add_K(double factor) {
  myMatrix->Add_cM(factor, myElement->get_K());
}
void StandardModelElement::add_M(double factor) {
  myMatrix->Add_cM(factor, myElement->get_M());
}
void StandardModelElement::add_C(double factor) {
  myMatrix->Add_cM(factor, myElement->get_C());
}
void StandardModelElement::add_R(double factor) {
  myVector->Add_cV(factor, myElement->get_R());
}
void StandardModelElement::add_Reff(double factor) {
  myVector->Add_cV(factor, myElement->get_Reff());
}
void StandardModelElement::add_Rgrad(double factor) {
  myVector->Add_cV(factor, myElement->get_Rgrad());
}
