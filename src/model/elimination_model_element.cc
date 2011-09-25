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


#include "model/elimination_model_element.h"

/**
 * Default constructor
 */
EliminationModelElement::EliminationModelElement()
  :ModelElement() {
}
/**
 * Constructor
 */
EliminationModelElement::EliminationModelElement(const IDContainer& FTable,
                       Element* pElement)
  :ModelElement(FTable, pElement, 0) {
  myMatrix = theStaticMatrices[FTable.size()];
  myVector = theStaticVectors[FTable.size()];
}
/**
 * Destructor
 */
EliminationModelElement::~EliminationModelElement() {
}

void EliminationModelElement::add_K(double factor) {
  myMatrix->Add_cM(factor, myElement->get_K());
}
void EliminationModelElement::add_M(double factor) {
  myMatrix->Add_cM(factor, myElement->get_M());
}
void EliminationModelElement::add_C(double factor) {
  myMatrix->Add_cM(factor, myElement->get_C());
}
void EliminationModelElement::add_R(double factor) {
  myVector->add_cV(factor, myElement->get_R());
}
void EliminationModelElement::add_Reff(double factor) {
  myVector->add_cV(factor, myElement->get_Reff());
}
void EliminationModelElement::add_Rgrad(double factor) {
  myVector->add_cV(factor, myElement->get_Rgrad());
}
